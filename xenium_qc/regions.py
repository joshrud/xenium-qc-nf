"""Annotation regions and broad cell types on QC-passing cells: the hand-off point to differential expression.

- Pathology GeoJSON (drawn on the H&E, in H&E pixels) -> Xenium µm via the H&E alignment matrix and pixel size.
- Cells assigned to regions by centroid; cells inside two regions are labelled ambiguous.
- Leiden clusters labelled with the highest-scoring marker set; composition and transcript depth per region.
- Xenium Explorer cell-group CSVs (cell_id,group).
"""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import shapely

from .common import Bundle, StepOutput, load_markers, sample_rows, setup_logging, timed


def load_annotations(geojson: str, matrix: np.ndarray | None, pixel_size: float, units: str = "he_px") -> pd.DataFrame:
    """Return one row per polygon with a unique region id, in Xenium µm."""
    feats = json.loads(open(geojson).read())["features"]

    def conv(ring):
        r = np.asarray(ring, float)[:, :2]
        if units == "um":
            return r
        if matrix is None:
            raise ValueError("annotation in H&E pixels needs an H&E alignment matrix")
        return (r @ matrix[:2, :2].T + matrix[:2, 2]) * pixel_size

    rows = []
    for f in feats:
        g, props = f["geometry"], f.get("properties", {})
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"] if g["type"] == "MultiPolygon" else []
        geom = shapely.unary_union([shapely.Polygon(conv(p[0]), [conv(h) for h in p[1:]]) for p in polys])
        cls = (props.get("classification") or {}).get("name")
        rows.append({"name": props.get("name") or cls or "region", "classification": cls or "Unclassified", "geometry": geom})
    ann = pd.DataFrame(rows)
    ann["area_mm2"] = [g.area / 1e6 for g in ann.geometry]
    ann["cx"] = [g.centroid.x for g in ann.geometry]
    ann = ann.sort_values(["name", "cx"]).reset_index(drop=True)
    ann["region"] = ann.name
    dup = ann.name.duplicated(keep=False)
    ann.loc[dup, "region"] = ann[dup].name + " " + (ann[dup].groupby("name").cumcount() + 1).astype(str)
    return ann


def write_geojson_um(ann: pd.DataFrame, path):
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"region": r.region, "name": r.name, "classification": r.classification,
                                           "area_mm2": round(r.area_mm2, 4)},
         "geometry": shapely.geometry.mapping(r.geometry)} for r in ann.itertuples()],
          "note": "Coordinates in Xenium microns"}
    path.write_text(json.dumps(fc))


def run(bundle: str, outdir: str, filtered_h5ad: str, annotation: str | None = None, he_alignment: str | None = None,
        annotation_units: str = "he_px", markers: str | None = None, leiden_resolution: float = 1.0,
        n_pcs: int = 30, n_neighbors: int = 15):
    log = setup_logging("regions", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "regions", log)
    M = out.metrics
    cfg = load_markers(markers)
    ad = sc.read_h5ad(filtered_h5ad)
    log.info("%d QC-passing cells x %d genes", ad.n_obs, ad.n_vars)

    # ---- cell types
    with timed(log, f"PCA ({n_pcs}), kNN ({n_neighbors}), Leiden (res {leiden_resolution})"):
        sc.pp.pca(ad, n_comps=n_pcs)
        sc.pp.neighbors(ad, n_neighbors=n_neighbors)
        sc.tl.leiden(ad, resolution=leiden_resolution, flavor="igraph", n_iterations=2, key_added="leiden")
    markers_ok = {k: [g for g in v if g in ad.var_names] for k, v in cfg.get("cell_types", {}).items()}
    markers_ok = {k: v for k, v in markers_ok.items() if len(v) >= 2}
    for t, g in markers_ok.items():
        sc.tl.score_genes(ad, g, score_name=f"score_{t}")
    cl = ad.obs.groupby("leiden", observed=True)[[f"score_{t}" for t in markers_ok]].mean()
    ad.obs["cell_type"] = ad.obs.leiden.map(cl.idxmax(1).str.replace("score_", "")).astype("category")
    lab = pd.DataFrame({"type": cl.idxmax(1).str.replace("score_", ""), "n_cells": ad.obs.leiden.value_counts().sort_index()})
    out.table(lab.rename_axis("leiden").reset_index(), "cluster_labels")
    M["cell_type_counts"] = ad.obs.cell_type.value_counts().to_dict()
    log.info("cell types: %s", M["cell_type_counts"])
    sc.pl.dotplot(ad, markers_ok, groupby="leiden", standard_scale="var", show=False)
    fig = plt.gcf()
    out.figure(fig, "cluster_markers")

    # ---- regions
    ad.obs["region"] = "unannotated"
    if annotation:
        matrix = np.loadtxt(he_alignment, delimiter=",", encoding="utf-8-sig") if he_alignment else None
        ann = load_annotations(annotation, matrix, b.pixel_size, annotation_units)
        write_geojson_um(ann, out.outdir / "annotations_xenium_um.geojson")
        out.table(ann.drop(columns="geometry").round(3), "annotation_regions")
        x, y = ad.obs.x_centroid.to_numpy(), ad.obs.y_centroid.to_numpy()
        inside = pd.DataFrame({r: shapely.contains_xy(g, x, y) for r, g in zip(ann.region, ann.geometry)},
                              index=ad.obs_names)
        n = inside.sum(1)
        ad.obs["region"] = np.where(n == 1, inside.idxmax(1), np.where(n > 1, "ambiguous", "unannotated"))
        M["region_counts"] = ad.obs.region.value_counts().to_dict()
        log.info("cells per region: %s", M["region_counts"])

        fig, ax = plt.subplots(figsize=(11, 6))
        s = sample_rows(ad.obs, 80_000)
        regs = [r for r in ann.region if (s.region == r).any()]
        cmap = plt.get_cmap("tab20")
        ax.scatter(s.x_centroid, s.y_centroid, s=0.2, c="lightgrey", rasterized=True)
        for i, r in enumerate(regs):
            d = s[s.region == r]
            ax.scatter(d.x_centroid, d.y_centroid, s=0.2, color=cmap(i % 20), label=r, rasterized=True)
        for r in ann.itertuples():
            for g in getattr(r.geometry, "geoms", [r.geometry]):
                ax.plot(*g.exterior.xy, lw=0.8, c="k")
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.legend(markerscale=20, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
        ax.set_title("Cells assigned to annotation regions (sample of 80k)")
        out.figure(fig, "regions_selection")

    # ---- composition and depth per region (depth differences bias naive DE between regions)
    grp = ad.obs[ad.obs.region.ne("ambiguous")]
    comp = pd.crosstab(grp.cell_type, grp.region, normalize="columns") * 100
    out.table(comp.round(3).reset_index(), "cell_type_composition")
    if comp.shape[1] > 1:
        ax = comp.T.plot.barh(stacked=True, figsize=(10, 0.45 * comp.shape[1] + 1.5), colormap="tab20", width=0.8)
        ax.legend(bbox_to_anchor=(1, 1), fontsize=8); ax.set_xlabel("% of cells in region")
        out.figure(ax.get_figure(), "cell_type_composition")
    depth = grp.pivot_table(index="cell_type", columns="region", values="total_counts", aggfunc="median", observed=True)
    out.table(depth.round(1).reset_index(), "depth_by_region_cell_type")
    region_depth = grp.groupby("region").total_counts.median()
    M["median_counts_by_region"] = region_depth.round(1).to_dict()
    if len(region_depth) > 1:
        M["depth_ratio_max_min"] = float(region_depth.max() / max(region_depth.min(), 1))
        fig, ax = plt.subplots(figsize=(8, 0.35 * len(region_depth) + 1.2))
        region_depth.sort_values().plot.barh(ax=ax)
        ax.set_xlabel("median transcripts per cell"); ax.set_title("Transcript depth by region (match depth before DE)")
        out.figure(fig, "depth_by_region")
        log.info("median transcripts/cell by region: %s", M["median_counts_by_region"])

    # ---- outputs for downstream DE and Xenium Explorer
    labels = ad.obs[["region", "cell_type", "leiden", "x_centroid", "y_centroid", "total_counts"]].rename_axis("cell_id").reset_index()
    labels.to_csv(out.outdir / "tables" / "cell_labels.csv", index=False)
    cg = out.outdir / "cell_groups"
    cg.mkdir(exist_ok=True)
    labels.assign(group=labels.cell_type)[["cell_id", "group"]].to_csv(cg / "cell_groups_cell_types.csv", index=False)
    if annotation:
        labels.assign(group=labels.region)[["cell_id", "group"]].to_csv(cg / "cell_groups_regions.csv", index=False)
        labels.assign(group=labels.region + " | " + labels.cell_type.astype(str))[["cell_id", "group"]].to_csv(
            cg / "cell_groups_region_by_cell_type.csv", index=False)
    log.info("wrote cell labels and Xenium Explorer cell groups")
    out.write_metrics()
