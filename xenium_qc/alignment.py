"""H&E-to-Xenium alignment QC: keypoint residuals, leave-one-out error, comparison with a reference transform,
and a ranking of panel genes that make good landmarks for placing keypoints.

Transforms follow Xenium Explorer's convention: a 3x3 matrix mapping H&E pixels -> DAPI (morphology) pixels.
Keypoint CSVs exported by Xenium Explorer have fixedX/fixedY in DAPI pixels and alignmentX/alignmentY in H&E pixels.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

from .common import Bundle, StepOutput, setup_logging, timed


def fit_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Least-squares similarity transform a -> b (Horn/Umeyama; reflection allowed), as a 3x3 matrix."""
    ma, mb = a.mean(0), b.mean(0)
    A, B = a - ma, b - mb
    U, S, Vt = np.linalg.svd(B.T @ A / len(a))
    R = U @ Vt
    s = S.sum() / (A ** 2).sum(1).mean()
    M = np.eye(3)
    M[:2, :2] = s * R
    M[:2, 2] = mb - s * R @ ma
    return M


def apply(M: np.ndarray, p: np.ndarray) -> np.ndarray:
    return p @ M[:2, :2].T + M[:2, 2]


def describe(M: np.ndarray) -> dict:
    return {"scale": float(np.sqrt(abs(np.linalg.det(M[:2, :2])))),
            "rotation_deg": float(np.degrees(np.arctan2(M[1, 0], M[0, 0]))),
            "reflected": bool(np.linalg.det(M[:2, :2]) < 0), "shift_x_px": float(M[0, 2]), "shift_y_px": float(M[1, 2])}


def _dapi_extent(b: Bundle) -> tuple[float, float]:
    """DAPI image width/height in pixels (from the image, or the transcript extent as a fallback)."""
    import tifffile

    ch = b.morphology_channels()
    if ch:
        with tifffile.TiffFile(next(iter(ch.values()))) as tf:
            shp = tf.series[0].levels[0].shape
            return float(shp[-1]), float(shp[-2])
    ext = b.transcripts().select(pl.col("x_location").max(), pl.col("y_location").max()).collect().row(0)
    return ext[0] / b.pixel_size, ext[1] / b.pixel_size


def rank_landmark_genes(b: Bundle, bin_um: int = 50, cover_um: int = 1000, min_per_tile: int = 200) -> pd.DataFrame:
    """Score genes by abundance x spatial contrast x tissue coverage: good keypoint-finding guides."""
    tx = (b.transcripts().filter((pl.col("qv") >= 20) & pl.col("is_gene"))
            .select(pl.col("feature_name").alias("gene"), (pl.col("x_location") // bin_um).cast(pl.Int32).alias("bx"),
                    (pl.col("y_location") // bin_um).cast(pl.Int32).alias("by"))
            .group_by("gene", "bx", "by").agg(n=pl.len()).collect().to_pandas())
    n_tissue = tx[["bx", "by"]].drop_duplicates().shape[0]
    cell = cover_um // bin_um
    tiles_all = tx.assign(tx_=tx.bx // cell, ty_=tx.by // cell)[["tx_", "ty_"]].drop_duplicates().shape[0]
    rows = []
    for gene, g in tx.groupby("gene"):
        c = np.sort(g.n.to_numpy())[::-1]
        k = max(1, int(0.05 * n_tissue))
        cov = (g.assign(tx_=g.bx // cell, ty_=g.by // cell).groupby(["tx_", "ty_"]).n.sum() >= min_per_tile).sum()
        rows.append((gene, int(c.sum()), c[:k].sum() / c.sum(), cov / tiles_all))
    df = pd.DataFrame(rows, columns=["gene", "count", "frac_in_top5pct_area", "tissue_coverage"])
    r = lambda s: s.rank() / len(df)  # noqa: E731
    df["score"] = r(df.frac_in_top5pct_area) * r(np.log10(df["count"])) * r(df.tissue_coverage)
    return df.sort_values("score", ascending=False)


def run(bundle: str, outdir: str, keypoints: str | None = None, he_alignment: str | None = None,
        reference_alignment: str | None = None):
    """keypoints: Xenium Explorer keypoints.csv. he_alignment: your exported matrix (optional if keypoints given).
    reference_alignment: a second matrix to compare against (e.g. the 10x supplemental alignment)."""
    log = setup_logging("alignment", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "alignment", log)
    M = out.metrics
    px = b.pixel_size

    with timed(log, "ranking landmark genes"):
        rank = rank_landmark_genes(b)
        out.table(rank, "keypoint_gene_ranking")
        M["top_landmark_genes"] = rank.gene.head(15).tolist()

    my_M = np.loadtxt(he_alignment, delimiter=",", encoding="utf-8-sig") if he_alignment else None
    ref_M = np.loadtxt(reference_alignment, delimiter=",", encoding="utf-8-sig") if reference_alignment else None
    transforms = {}

    if keypoints:
        kp = pd.read_csv(keypoints)
        he = kp[["alignmentX", "alignmentY"]].to_numpy(float)
        fx = kp[["fixedX", "fixedY"]].to_numpy(float)
        fit = fit_similarity(he, fx)
        transforms["keypoint fit"] = fit
        if my_M is None:
            my_M = fit
        err = np.linalg.norm(apply(fit, he) - fx, axis=1) * px
        loo = np.array([np.linalg.norm(apply(fit_similarity(np.delete(he, i, 0), np.delete(fx, i, 0)), he[i:i + 1])[0] - fx[i]) * px
                        for i in range(len(he))])
        tab = pd.DataFrame({"keypoint": np.arange(1, len(he) + 1), "dapi_x_um": fx[:, 0] * px, "dapi_y_um": fx[:, 1] * px,
                            "residual_um": err, "loo_error_um": loo})
        if ref_M is not None:
            tab["reference_error_um"] = np.linalg.norm(apply(ref_M, he) - fx, axis=1) * px
        out.table(tab.round(3), "alignment_keypoint_errors")
        rmse = lambda e: float(np.sqrt(np.mean(np.square(e))))  # noqa: E731
        M.update({"n_keypoints": int(len(he)), "residual_rmse_um": rmse(err), "loo_median_um": float(np.median(loo)),
                  "loo_max_um": float(loo.max()), "worst_keypoint": int(tab.keypoint[loo.argmax()]),
                  "reference_rmse_um": rmse(tab.reference_error_um) if ref_M is not None else None})
        # nearest-neighbour spacing: are keypoints spread or clustered?
        d = np.linalg.norm((fx[:, None] - fx[None]) * px, axis=-1); np.fill_diagonal(d, np.inf)
        M["keypoint_nn_spacing_median_um"] = float(np.median(d.min(1)))
        log.info("keypoints %d | residual RMSE %.2f µm | LOO median %.2f, max %.2f µm (keypoint %d)",
                 len(he), M["residual_rmse_um"], M["loo_median_um"], M["loo_max_um"], M["worst_keypoint"])

        dens = (b.transcripts().filter(pl.col("is_gene") & (pl.int_range(pl.len()) % 200 == 0))
                  .select("x_location", "y_location").collect().to_pandas())
        fig, axs = plt.subplots(1, 2, figsize=(15, 4.5), gridspec_kw={"width_ratios": [2.3, 1]})
        axs[0].scatter(dens.x_location, dens.y_location, s=0.05, c="lightgrey", rasterized=True)
        sca = axs[0].scatter(tab.dapi_x_um, tab.dapi_y_um, c=tab.loo_error_um, s=90, cmap="plasma", edgecolors="k",
                             vmin=0, vmax=max(5, tab.loo_error_um.max()))
        for _, r in tab.iterrows():
            axs[0].annotate(int(r.keypoint), (r.dapi_x_um, r.dapi_y_um), xytext=(4, 4), textcoords="offset points", fontsize=7)
        axs[0].set_aspect("equal"); axs[0].invert_yaxis(); fig.colorbar(sca, ax=axs[0], shrink=0.7, label="LOO error (µm)")
        axs[0].set_title("Keypoints over the tissue, coloured by leave-one-out error")
        w = 0.4
        axs[1].bar(tab.keypoint - w / 2, tab.residual_um, w, label="fit residual")
        axs[1].bar(tab.keypoint + w / 2, tab.loo_error_um, w, label="leave-one-out")
        if ref_M is not None:
            axs[1].plot(tab.keypoint, tab.reference_error_um, "k.", label="reference matrix")
        axs[1].axhline(5, color="r", ls="--", lw=0.8); axs[1].set_xlabel("keypoint"); axs[1].set_ylabel("error (µm)")
        axs[1].legend(fontsize=8); axs[1].set_title("Per-keypoint error (red: ~half a nucleus, 5 µm)")
        out.figure(fig, "alignment_keypoint_errors")

    if my_M is not None:
        transforms.setdefault("my alignment", my_M)
    if ref_M is not None:
        transforms["reference alignment"] = ref_M
    if transforms:
        tt = pd.DataFrame({k: describe(v) for k, v in transforms.items()}).T.rename_axis("transform").reset_index()
        out.table(tt, "alignment_transforms")
        M["transforms"] = {k: describe(v) for k, v in transforms.items()}

    if my_M is not None and ref_M is not None:
        with timed(log, "comparing transforms over the capture area"):
            W, H = _dapi_extent(b)
            gx, gy = np.meshgrid(np.linspace(0, W, 80), np.linspace(0, H, 40))
            grid = np.c_[gx.ravel(), gy.ravel()]
            gap = np.linalg.norm(apply(my_M, apply(np.linalg.inv(ref_M), grid)) - grid, axis=1) * px
            corners = np.array([[0, 0], [W, 0], [0, H], [W, H]], float)
            gap_c = np.linalg.norm(apply(my_M, apply(np.linalg.inv(ref_M), corners)) - corners, axis=1) * px
            M.update({"vs_reference_median_um": float(np.median(gap)), "vs_reference_max_um": float(gap.max()),
                      "vs_reference_corners_um": gap_c.round(2).tolist()})
            log.info("my vs reference transform: median %.2f µm, max %.2f µm", np.median(gap), gap.max())
            fig, ax = plt.subplots(figsize=(10, 4.5))
            im = ax.imshow(gap.reshape(gy.shape), extent=(0, W * px, H * px, 0), cmap="viridis")
            fig.colorbar(im, ax=ax, shrink=0.8, label="displacement (µm)")
            ax.set_title("Difference between my and the reference transform across the capture area")
            out.figure(fig, "alignment_transform_difference")
    if not keypoints and not he_alignment:
        log.info("no keypoints or alignment matrix supplied; only the landmark-gene ranking was produced")
    out.write_metrics()
