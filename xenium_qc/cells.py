"""Cell-level QC: segmentation, per-cell distributions, filters, gene filter, spillover, normalization."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

from .common import Bundle, StepOutput, load_markers, sample_rows, setup_logging, timed


def _mad_outlier(x, k=5):
    x = np.log10(x + 1)
    med, mad = np.median(x), stats.median_abs_deviation(x, scale="normal")
    return (x < med - k * mad) | (x > med + k * mad)


def _leak(X, a, b):
    """Mean count of B in A-high cells, relative to B in B-high cells."""
    ref = X.loc[X[b] >= 3, b].mean()
    return X.loc[X[a] >= 3, b].mean() / ref if ref else np.nan


def run(bundle: str, outdir: str, markers: str | None = None, min_counts: int = 10, min_cells_per_gene: int = 5,
        mad_k: float = 5.0, target_sum: float = 100.0, write_h5ad: bool = True):
    log = setup_logging("cells", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "cells", log)
    M = out.metrics
    cfg = load_markers(markers)

    with timed(log, "loading cell_feature_matrix.h5 and cells.parquet"):
        raw = sc.read_10x_h5(b.file("cell_feature_matrix.h5"), gex_only=False)
        cells = b.cells()
        raw = raw[raw.obs_names.isin(cells.index)].copy()
        ftype = raw.var["feature_types"]
        rna = raw[:, ftype == "Gene Expression"].copy()
        rna.var_names_make_unique()
        rna.obs = cells.loc[rna.obs_names].copy()
        rna.obsm["spatial"] = rna.obs[["x_centroid", "y_centroid"]].to_numpy()
        log.info("%d cells x %d genes; feature types: %s", rna.n_obs, rna.n_vars, ftype.value_counts().to_dict())

    sc.pp.calculate_qc_metrics(rna, percent_top=None, log1p=False, inplace=True)
    obs = rna.obs
    obs["nucleus_frac"] = obs.nucleus_area / obs.cell_area
    ctrl_cols = [c for c in ["control_probe_counts", "control_codeword_counts", "unassigned_codeword_counts",
                             "deprecated_codeword_counts", "genomic_control_counts"] if c in obs]
    obs["control_counts"] = obs[ctrl_cols].sum(axis=1)

    # ---- segmentation method
    if "segmentation_method" in obs:
        seg = (obs.groupby("segmentation_method", observed=True)
                  .agg(n_cells=("total_counts", "size"), median_counts=("total_counts", "median"),
                       median_genes=("n_genes_by_counts", "median"), median_area=("cell_area", "median")))
        seg["fraction"] = seg.n_cells / seg.n_cells.sum()
        out.table(seg.reset_index(), "qc_segmentation_methods")
        M["segmentation_fractions"] = seg.fraction.round(4).to_dict()
        exp = [k for k in seg.index if "expansion" in str(k).lower()]
        M["nuclear_expansion_fraction"] = float(seg.loc[exp, "fraction"].sum()) if exp else 0.0
        log.info("segmentation methods: %s", M["segmentation_fractions"])

    # ---- distributions
    with timed(log, "distributions"):
        fig, axs = plt.subplots(1, 4, figsize=(17, 3))
        for ax, (col, logx) in zip(axs, [("total_counts", True), ("n_genes_by_counts", True),
                                         ("cell_area", True), ("nucleus_frac", False)]):
            v = obs[col].dropna()
            ax.hist(np.log10(v + 1) if logx else v, bins=100)
            ax.set_xlabel(f"log10({col} + 1)" if logx else col)
        axs[0].axvline(np.log10(min_counts + 1), color="r", ls="--", lw=1)
        fig.suptitle(f"Per-cell QC distributions (red: {min_counts}-count cutoff)")
        out.figure(fig, "qc_cell_distributions")

        sub = sample_rows(obs, 100_000)
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.hexbin(np.log10(sub.cell_area), np.log10(sub.total_counts + 1), gridsize=80, bins="log", cmap="viridis")
        ax.set_xlabel("log10 cell area (µm²)"); ax.set_ylabel("log10(counts + 1)"); ax.set_title("Counts vs area")
        out.figure(fig, "qc_counts_vs_area")

        sub = sample_rows(obs, 150_000)
        fig, axs = plt.subplots(1, 2, figsize=(16, 4))
        for ax, col, title in [(axs[0], "total_counts", "log10 counts per cell"), (axs[1], "cell_area", "log10 cell area")]:
            v = np.log10(sub[col] + 1)
            sca = ax.scatter(sub.x_centroid, sub.y_centroid, c=v, s=0.2, cmap="magma", rasterized=True,
                             vmin=v.quantile(0.02), vmax=v.quantile(0.98))
            ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(title); fig.colorbar(sca, ax=ax, shrink=0.7)
        out.figure(fig, "qc_spatial_maps")
        M.update({"median_counts": float(obs.total_counts.median()), "median_genes": float(obs.n_genes_by_counts.median()),
                  "median_cell_area": float(obs.cell_area.median())})

    # ---- filters
    obs["fail_low_counts"] = obs.total_counts < min_counts
    obs["fail_no_nucleus"] = obs.nucleus_area.isna() | (obs.nucleus_area <= 0)
    obs["flag_area_outlier"] = _mad_outlier(obs.cell_area, mad_k)
    obs["flag_counts_outlier"] = _mad_outlier(obs.total_counts, mad_k) & ~obs.fail_low_counts
    obs["flag_any_control"] = obs.control_counts > 0
    obs["pass_qc"] = ~(obs.fail_low_counts | obs.fail_no_nucleus)
    flag_cols = [c for c in obs.columns if c.startswith(("fail_", "flag_", "pass_"))]
    fs = obs[flag_cols].agg(["sum", "mean"]).T.rename(columns={"sum": "n_cells", "mean": "fraction"})
    fs["n_cells"] = fs["n_cells"].astype(int)
    out.table(fs.rename_axis("filter").reset_index(), "qc_cell_filter_summary")
    M["filters"] = fs.to_dict(orient="index")
    M["n_cells"] = int(len(obs)); M["n_pass"] = int(obs.pass_qc.sum())
    log.info("cells passing hard filters: %d / %d (%.2f%%)", M["n_pass"], M["n_cells"], 100 * obs.pass_qc.mean())

    s = sample_rows(obs, 150_000)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.scatter(s.x_centroid, s.y_centroid, s=0.1, c="lightgrey", rasterized=True)
    f = s[~s.pass_qc]
    ax.scatter(f.x_centroid, f.y_centroid, s=0.5, c="red", rasterized=True, label="removed")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.legend(markerscale=10)
    ax.set_title(f"Cells removed by hard filters ({(~obs.pass_qc).mean():.2%})")
    out.figure(fig, "qc_removed_cells")

    keep = ["x_centroid", "y_centroid", "transcript_counts", "total_counts", "n_genes_by_counts", "cell_area",
            "nucleus_area", "nucleus_frac", "control_counts", "segmentation_method"] + flag_cols
    obs[[c for c in keep if c in obs]].rename_axis("cell_id").reset_index().to_parquet(
        out.outdir / "cells_qc.parquet", index=False)
    log.info("wrote cells_qc.parquet")

    # ---- gene filter
    filt = rna[obs.pass_qc.to_numpy()].copy()
    n_before = filt.n_vars
    sc.pp.filter_genes(filt, min_cells=min_cells_per_gene)
    M.update({"genes_before_filter": int(n_before), "genes_after_filter": int(filt.n_vars),
              "genes_removed": sorted(set(rna.var_names) - set(filt.var_names))})
    log.info("genes: %d -> %d (min %d cells)", n_before, filt.n_vars, min_cells_per_gene)

    # ---- segmentation spillover
    with timed(log, "segmentation spillover"):
        excl = [tuple(p) for p in cfg.get("spillover_exclusive", []) if all(g in filt.var_names for g in p)]
        pos = [tuple(p) for p in cfg.get("spillover_positive", []) if all(g in filt.var_names for g in p)]
        if excl:
            spill_genes = sorted({g for pr in excl + pos for g in pr})
            nuc = (b.transcripts().filter((pl.col("qv") >= 20) & pl.col("feature_name").is_in(spill_genes)
                                          & (pl.col("overlaps_nucleus") == 1) & (pl.col("cell_id") != "UNASSIGNED"))
                     .group_by("cell_id", "feature_name").agg(n=pl.len()).collect().to_pandas()
                     .pivot(index="cell_id", columns="feature_name", values="n"))
            W = filt[:, spill_genes].X
            whole = pd.DataFrame(W.toarray() if sp.issparse(W) else W, index=filt.obs_names, columns=spill_genes)
            nuc = nuc.reindex(index=whole.index, columns=spill_genes).fillna(0)
            spill = pd.DataFrame([{"A_high": a, "B": c, "type": "positive control" if (a, c) in pos else "should be exclusive",
                                   "n_A_high": int((whole[a] >= 3).sum()),
                                   "leak_whole_cell": _leak(whole, a, c), "leak_nucleus_only": _leak(nuc, a, c)}
                                  for a, c in excl + pos])
            spill["nucleus_over_whole"] = spill.leak_nucleus_only / spill.leak_whole_cell
            out.table(spill, "qc_segmentation_spillover")
            M["spillover"] = spill.round(3).to_dict(orient="records")
        else:
            log.warning("no spillover gene pairs present in this panel; skipped")

    # ---- normalization (Salas et al. 2025: library size to 100, log1p)
    if write_h5ad:
        with timed(log, "normalizing and writing filtered.h5ad"):
            filt.layers["counts"] = filt.X.copy()
            sc.pp.normalize_total(filt, target_sum=target_sum)
            sc.pp.log1p(filt)
            filt.obs = filt.obs[[c for c in keep if c in filt.obs]].copy()
            filt.uns["xenium_qc"] = {"min_counts": min_counts, "min_cells_per_gene": min_cells_per_gene,
                                     "target_sum": target_sum}
            filt.write_h5ad(out.outdir / "filtered.h5ad", compression="gzip")
    out.write_metrics()
