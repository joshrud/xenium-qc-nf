"""Transcript-level QC on transcripts.parquet: quality, background, z, FOVs, density, seams, hotspots, assignment."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

from .common import Bundle, StepOutput, setup_logging, show_grid, timed, to_grid

CONTROL_CATS = ["negative_control_probe", "negative_control_codeword"]


def _mad_low(s, k=3):
    return s < s.median() - k * stats.median_abs_deviation(s, scale="normal")


def _mad_high(s, k=3):
    return s > s.median() + k * stats.median_abs_deviation(s, scale="normal")


def run(bundle: str, outdir: str, bin_um: int = 20, ctrl_bin_um: int = 200, strip_um: int = 4,
        seam_window: int = 31, z_sample: int = 2_000_000):
    log = setup_logging("transcripts", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "transcripts", log)
    M = out.metrics
    lf = b.transcripts()
    q20 = lf.filter(pl.col("qv") >= 20)
    genes_q20 = q20.filter(pl.col("is_gene"))

    # ---- 1. codeword categories and QV distribution
    with timed(log, "codeword categories"):
        by_cat = (lf.group_by("codeword_category")
                    .agg(n=pl.len(), frac_q20=(pl.col("qv") >= 20).mean(), n_q20=(pl.col("qv") >= 20).sum(),
                         frac_in_cells=(pl.col("cell_id") != "UNASSIGNED").mean(),
                         frac_in_nucleus=(pl.col("overlaps_nucleus") == 1).mean())
                    .sort("n", descending=True).collect().to_pandas().set_index("codeword_category"))
        by_cat["frac_of_all"] = by_cat["n"] / by_cat["n"].sum()
        out.table(by_cat.reset_index(), "qc_codeword_categories")
        gene_cats = [c for c in by_cat.index if c.endswith("_gene")]
        M["n_transcripts"] = int(by_cat.n.sum())
        M["gene_frac_q20"] = float((by_cat.loc[gene_cats, "n_q20"].sum()) / by_cat.loc[gene_cats, "n"].sum())
        ctrl = [c for c in CONTROL_CATS if c in by_cat.index]
        M["control_frac_q20"] = float(by_cat.loc[ctrl, "n_q20"].sum() / by_cat.loc[ctrl, "n"].sum()) if ctrl else None
        M["n_control_q20"] = int(by_cat.loc[ctrl, "n_q20"].sum()) if ctrl else 0
        log.info("gene transcripts >=Q20: %.3f; control codewords/probes >=Q20: %s", M["gene_frac_q20"], M["control_frac_q20"])

        qv = (lf.select(pl.col("qv").round(0).clip(0, 40).alias("qv"), "is_gene")
                .group_by("qv", "is_gene").agg(n=pl.len()).collect().to_pandas())
        fig, ax = plt.subplots(figsize=(6, 3))
        for is_gene, label in [(True, "gene codewords"), (False, "control codewords")]:
            h = qv[qv.is_gene == is_gene].sort_values("qv")
            ax.plot(h.qv, h.n / h.n.sum(), label=label)
        ax.axvline(20, color="k", ls="--", lw=1); ax.set_xlabel("QV"); ax.set_ylabel("fraction"); ax.legend()
        ax.set_title("QV distribution (dashed line: Q20 cutoff)")
        out.figure(fig, "qc_qv_distribution")

    # ---- 2. genes against negative-control background
    with timed(log, "genes vs negative controls"):
        per_feature = (q20.group_by("feature_name", "is_gene")
                          .agg(n=pl.len(), frac_in_cells=(pl.col("cell_id") != "UNASSIGNED").mean(),
                               codeword_category=pl.col("codeword_category").unique().sort().str.join("+"))
                          .collect().to_pandas())
        controls = per_feature[per_feature.codeword_category.isin(CONTROL_CATS)]
        genes = per_feature[per_feature.is_gene].copy()
        bg_q99 = controls["n"].quantile(0.99) if len(controls) else 0
        genes["fold_over_control_max"] = genes["n"] / max(controls["n"].max() if len(controls) else 1, 1)
        genes["above_background"] = genes["n"] > bg_q99
        out.table(genes.sort_values("n"), "qc_gene_detection")
        low_bg = sorted(genes.loc[~genes.above_background, "feature_name"])
        M.update({"n_genes": int(len(genes)), "n_genes_above_background": int(genes.above_background.sum()),
                  "genes_below_background": low_bg, "control_q99_counts": float(bg_q99)})
        log.info("genes above 99th pct of controls: %d / %d; below: %s", genes.above_background.sum(), len(genes), low_bg)
        fig, ax = plt.subplots(figsize=(7, 3))
        bins = np.logspace(0, np.log10(max(per_feature.n.max(), 10)) + 0.1, 60)
        ax.hist(genes.n, bins=bins, alpha=0.7, label="genes")
        ax.hist(controls.n, bins=bins, alpha=0.7, label="negative controls")
        ax.axvline(max(bg_q99, 1), color="k", ls="--", lw=1)
        ax.set_xscale("log"); ax.set_xlabel("total Q20 counts per feature"); ax.legend()
        ax.set_title("Per-feature counts: genes vs negative controls")
        out.figure(fig, "qc_gene_vs_negative_controls")

    # ---- 3. z distribution
    with timed(log, "z distribution"):
        zs = genes_q20.select("z_location").collect().to_series()
        zs = zs.sample(min(z_sample, len(zs)), seed=0).to_numpy()
        lo, hi = np.percentile(zs, [5, 95])
        M["z_span_5_95_um"] = float(hi - lo)
        thick = b.metrics.get("thickness_transcripts_high_quality", np.nan)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.hist(zs, bins=150, density=True)
        ax.axvline(lo, color="k", ls="--", lw=1); ax.axvline(hi, color="k", ls="--", lw=1)
        ax.set_xlabel("z (µm)")
        ax.set_title(f"Q20 gene transcripts in z: 5-95% span {hi - lo:.1f} µm (10x thickness {thick:.2f} µm)")
        out.figure(fig, "qc_z_distribution")

    # ---- 4. per-FOV consistency
    with timed(log, "per-FOV QC"):
        fov = (lf.group_by("fov_name")
                 .agg(n_all=pl.len(),
                      n_gene_q20=((pl.col("qv") >= 20) & pl.col("is_gene")).sum(),
                      frac_q20=(pl.col("qv") >= 20).filter(pl.col("is_gene")).mean(),
                      frac_control_q20=((pl.col("qv") >= 20) & ~pl.col("is_gene")).sum() / (pl.col("qv") >= 20).sum(),
                      x=pl.col("x_location").median(), y=pl.col("y_location").median(),
                      x_min=pl.col("x_location").min(), x_max=pl.col("x_location").max(),
                      y_min=pl.col("y_location").min(), y_max=pl.col("y_location").max())
                 .collect().to_pandas().set_index("fov_name"))
        fov["flag_low_q20"] = _mad_low(fov.frac_q20)
        fov["flag_low_density"] = _mad_low(np.log10(fov.n_gene_q20 + 1))
        fov["flag_high_control"] = _mad_high(fov.frac_control_q20)
        out.table(fov.reset_index(), "qc_fov")
        M.update({"n_fov": int(len(fov)), "fov_low_q20": int(fov.flag_low_q20.sum()),
                  "fov_low_density": int(fov.flag_low_density.sum()), "fov_high_control": int(fov.flag_high_control.sum())})
        log.info("%d FOVs; flagged low Q20 %d, low density %d, high control %d", len(fov),
                 M["fov_low_q20"], M["fov_low_density"], M["fov_high_control"])
        fig, axs = plt.subplots(1, 3, figsize=(16, 3.5))
        for ax, col, title in zip(axs, ["frac_q20", "n_gene_q20", "frac_control_q20"],
                                  ["Q20 fraction", "Q20 gene transcripts", "control fraction (Q20)"]):
            sca = ax.scatter(fov.x, fov.y, c=fov[col], s=25, marker="s", cmap="viridis")
            fl = fov[fov.flag_low_q20 | fov.flag_high_control]
            ax.scatter(fl.x, fl.y, s=60, facecolors="none", edgecolors="red", marker="s")
            ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(title); fig.colorbar(sca, ax=ax, shrink=0.7)
        fig.suptitle("Per-FOV QC (red: flagged for Q20 or control fraction)")
        out.figure(fig, "qc_fov_maps")

    # ---- 5. density, z and assignment maps on a 20 µm grid
    with timed(log, f"{bin_um} µm bins"):
        bins = (genes_q20.with_columns(bx=(pl.col("x_location") // bin_um).cast(pl.Int32),
                                       by=(pl.col("y_location") // bin_um).cast(pl.Int32))
                  .group_by("bx", "by")
                  .agg(n=pl.len(), assigned=(pl.col("cell_id") != "UNASSIGNED").mean(),
                       z_med=pl.col("z_location").median(),
                       z_iqr=pl.col("z_location").quantile(0.75) - pl.col("z_location").quantile(0.25))
                  .collect().to_pandas())
        in_tissue = bins.n >= bin_um * bin_um / 20  # >= 1 transcript per 20 µm² on average
        bins["log_density"] = np.log10(bins.n / (bin_um * bin_um))
        M["tissue_area_mm2"] = float(in_tissue.sum() * bin_um * bin_um / 1e6)
        log.info("tissue bins: %d (%.1f mm²)", in_tissue.sum(), M["tissue_area_mm2"])
        fig, ax = plt.subplots(figsize=(14, 7))
        show_grid(ax, to_grid(bins, "log_density"), bin_um,
                  f"log10 Q20 gene transcripts per µm² ({bin_um} µm bins), FOV outlines in white",
                  vmin=bins.log_density[in_tissue].quantile(0.01), vmax=bins.log_density.quantile(0.995))
        for _, r in fov.iterrows():
            ax.add_patch(plt.Rectangle((r.x_min, r.y_min), r.x_max - r.x_min, r.y_max - r.y_min,
                                       fill=False, ec="white", lw=0.4, alpha=0.6))
        out.figure(fig, "qc_transcript_density")

    # ---- 6. stitching seams: thin-strip density, detrended by a running median, vs distance to FOV borders
    with timed(log, "seam profile"):
        def edges_of(lo, hi):
            e = np.unique(np.round(np.concatenate([fov[lo], fov[hi]]) / strip_um) * strip_um)
            return e[(e > e.min()) & (e < e.max())]

        def strip_profile(coord, edges):
            c = (genes_q20.group_by((pl.col(coord) // strip_um).cast(pl.Int32).alias("s")).agg(n=pl.len())
                   .collect().to_pandas().set_index("s").n.sort_index())
            c = c.reindex(range(c.index.min(), c.index.max() + 1), fill_value=0)
            rel = c / c.rolling(seam_window, center=True, min_periods=10).median()
            pos = (c.index.to_numpy() + 0.5) * strip_um
            if len(edges) == 0:
                return pd.Series(dtype=float)
            d = pos[:, None] - edges[None, :]
            offset = d[np.arange(len(d)), np.abs(d).argmin(1)]
            prof = pd.DataFrame({"offset": np.round(offset / strip_um) * strip_um, "rel": rel.to_numpy()})
            prof = prof[np.isfinite(prof.rel) & (prof.offset.abs() <= 60)]
            return prof.groupby("offset").rel.mean()

        px = strip_profile("x_location", edges_of("x_min", "x_max"))
        py = strip_profile("y_location", edges_of("y_min", "y_max"))
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(px.index, px, label="x borders"); ax.plot(py.index, py, label="y borders")
        ax.axhline(1, color="k", lw=0.5); ax.axvline(0, color="k", ls="--", lw=0.5); ax.set_ylim(0.8, 1.2)
        ax.set_xlabel("distance to nearest FOV border (µm)"); ax.set_ylabel("density / local median"); ax.legend()
        ax.set_title("Seam check in slide coordinates")
        out.figure(fig, "qc_fov_seam_profile")
        near = lambda p: p[p.index.to_series().abs() <= 8] if len(p) else pd.Series([1.0])  # noqa: E731
        M["seam_max_deviation"] = float(max((near(px) - 1).abs().max(), (near(py) - 1).abs().max()))
        log.info("largest density deviation within 8 µm of a FOV border: %.1f%%", 100 * M["seam_max_deviation"])

    # ---- 7. negative-control hotspots
    with timed(log, "negative-control hotspots"):
        cb = (lf.with_columns(bx=(pl.col("x_location") // ctrl_bin_um).cast(pl.Int32),
                              by=(pl.col("y_location") // ctrl_bin_um).cast(pl.Int32))
                .group_by("bx", "by").agg(n=pl.len(), n_ctrl=pl.col("codeword_category").is_in(CONTROL_CATS).sum())
                .collect().to_pandas())
        rate = cb.n_ctrl.sum() / cb.n.sum()
        cb["expected"] = cb.n * rate
        cb["obs_over_exp"] = cb.n_ctrl / cb.expected
        cb["p"] = stats.poisson.sf(cb.n_ctrl - 1, cb.expected)
        cb["hotspot"] = (cb.p < 0.001 / len(cb)) & (cb.obs_over_exp > 2) & (cb.n > 1000)
        out.table(cb.assign(x_um=cb.bx * ctrl_bin_um, y_um=cb.by * ctrl_bin_um)
                    .sort_values("obs_over_exp", ascending=False), "qc_negative_control_bins")
        M.update({"control_hotspots": int(cb.hotspot.sum()), "control_bins": int(len(cb)), "control_rate": float(rate)})
        log.info("global control fraction %.2e; hotspot bins %d of %d", rate, cb.hotspot.sum(), len(cb))
        fig, ax = plt.subplots(figsize=(12, 5.5))
        g = np.full((cb.by.max() + 1, cb.bx.max() + 1), np.nan)
        ok = cb.n > 1000
        g[cb.by[ok], cb.bx[ok]] = np.log2((cb.n_ctrl[ok] + 0.5) / cb.expected[ok])
        im = ax.imshow(g, cmap="RdBu_r", vmin=-2, vmax=2, extent=(0, g.shape[1] * ctrl_bin_um, g.shape[0] * ctrl_bin_um, 0))
        hs = cb[cb.hotspot]
        ax.scatter((hs.bx + 0.5) * ctrl_bin_um, (hs.by + 0.5) * ctrl_bin_um, s=80, facecolors="none", edgecolors="k")
        plt.colorbar(im, ax=ax, shrink=0.7, label="log2 observed / expected")
        ax.set_title(f"Negative-control enrichment per {ctrl_bin_um} µm bin (circles: significant hotspots)")
        out.figure(fig, "qc_negative_control_hotspots")

    # ---- 8. section tilt and local thickness
    with timed(log, "z maps"):
        tb = bins[in_tissue].copy()
        med, mad = tb.z_iqr.median(), stats.median_abs_deviation(tb.z_iqr, scale="normal")
        tb["thick"] = tb.z_iqr > med + 5 * mad
        M.update({"z_iqr_median_um": float(med), "thick_bin_fraction": float(tb.thick.mean())})
        log.info("z IQR median %.2f µm; thick bins %.2f%% of tissue", med, 100 * tb.thick.mean())
        fig, axs = plt.subplots(1, 2, figsize=(18, 4.5))
        show_grid(axs[0], to_grid(tb, "z_med"), bin_um, "median z (µm): tilt / warping", cmap="coolwarm")
        show_grid(axs[1], to_grid(tb, "z_iqr"), bin_um, "z IQR (µm): local thickness",
                  vmin=tb.z_iqr.quantile(0.01), vmax=tb.z_iqr.quantile(0.995))
        out.figure(fig, "qc_z_maps")

    # ---- 9. transcript assignment to cells
    with timed(log, "assignment"):
        fig, axs = plt.subplots(1, 2, figsize=(18, 4.5), gridspec_kw={"width_ratios": [2.2, 1]})
        show_grid(axs[0], to_grid(bins[in_tissue], "assigned"), bin_um, "fraction of Q20 transcripts assigned to a cell",
                  vmin=0, vmax=1)
        nd = (genes_q20.select((pl.col("cell_id") != "UNASSIGNED").alias("assigned"),
                               pl.col("nucleus_distance").clip(0, 30))
                .filter(pl.int_range(pl.len()) % 50 == 0).collect().to_pandas())
        for a, lab in [(True, "assigned"), (False, "unassigned")]:
            axs[1].hist(nd.nucleus_distance[nd.assigned == a], bins=60, density=True, histtype="step", label=lab)
        axs[1].set_xlabel("distance to nearest nucleus (µm, clipped at 30)"); axs[1].legend()
        axs[1].set_title("Distance to nearest nucleus")
        out.figure(fig, "qc_assignment")
        M["assigned_fraction_q20_genes"] = float(nd.assigned.mean())
        M["unassigned_within_5um_of_nucleus"] = float((nd.nucleus_distance[~nd.assigned] <= 5).mean())
        M["lowest_assigned_genes"] = (genes.set_index("feature_name").frac_in_cells.sort_values()
                                      .head(10).round(3).to_dict())
    out.write_metrics()
