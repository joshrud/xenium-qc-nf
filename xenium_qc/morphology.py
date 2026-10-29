"""Marker transcripts against morphology: on-stain enrichment and marker crops (DAPI/protein/optional H&E)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import tifffile
from scipy.ndimage import map_coordinates

from .common import Bundle, StepOutput, load_markers, setup_logging, timed


def _read_region(page, y0, y1, x0, x1):
    """Decode only the tiles of a tiled TIFF page overlapping [y0:y1, x0:x1]."""
    y0, x0 = max(y0, 0), max(x0, 0)
    y1, x1 = min(y1, page.shape[0]), min(x1, page.shape[1])
    if not page.is_tiled:
        return page.asarray()[y0:y1, x0:x1]
    th, tw = page.tilelength, page.tilewidth
    ntx = -(-page.shape[1] // tw)
    out = np.zeros((y1 - y0, x1 - x0) + page.shape[2:], page.dtype)
    fh = page.parent.filehandle
    for ty in range(y0 // th, (y1 - 1) // th + 1):
        for tx in range(x0 // tw, (x1 - 1) // tw + 1):
            i = ty * ntx + tx
            fh.seek(page.dataoffsets[i])
            tile = page.decode(fh.read(page.databytecounts[i]), i)[0].reshape((th, tw) + page.shape[2:])
            ys, xs = ty * th, tx * tw
            sy0, sy1, sx0, sx1 = max(y0, ys), min(y1, ys + th), max(x0, xs), min(x1, xs + tw)
            out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = tile[sy0 - ys:sy1 - ys, sx0 - xs:sx1 - xs]
    return out


def _norm(a, lo=1, hi=99.5):
    a0, a1 = np.percentile(a, [lo, hi])
    return np.clip((a - a0) / (a1 - a0 + 1e-9), 0, 1)


class Morphology:
    def __init__(self, b: Bundle):
        self.b, self.channels = b, b.morphology_channels()
        self.names = list(self.channels)
        self._lower = {n.lower(): n for n in self.names}

    def resolve(self, name: str | None) -> str | None:
        """Case-insensitive channel lookup ('CD20' matches a file-derived 'cd20')."""
        return self._lower.get(str(name).lower()) if name else None

    def _page(self, tf, name, level):
        series = tf.series[0]
        lev = series.levels[min(level, len(series.levels) - 1)]
        pages = lev.pages
        return pages[self.names.index(name)] if len(pages) > 1 else pages[0]

    def level_image(self, name, level=2):
        with tifffile.TiffFile(self.channels[name]) as tf:
            return self._page(tf, name, level).asarray()

    def crop(self, name, x0, y0, size):
        px = self.b.pixel_size
        with tifffile.TiffFile(self.channels[name]) as tf:
            return _read_region(self._page(tf, name, 0), int(y0 / px), int((y0 + size) / px),
                                int(x0 / px), int((x0 + size) / px)).astype(float)


def _he_crop(he_path, align_inv, px, x0, y0, size, out_px=500):
    s = size / out_px
    cc, rr = np.meshgrid(np.arange(out_px), np.arange(out_px))
    dapi = np.stack([(x0 + cc * s) / px, (y0 + rr * s) / px, np.ones(cc.shape)])
    hx, hy, _ = np.tensordot(align_inv, dapi, axes=1)
    with tifffile.TiffFile(he_path) as tf:
        page = tf.pages[0]
        bx0, by0 = max(int(hx.min()) - 2, 0), max(int(hy.min()) - 2, 0)
        bx1, by1 = min(int(hx.max()) + 3, page.shape[1]), min(int(hy.max()) + 3, page.shape[0])
        region = _read_region(page, by0, by1, bx0, bx1)
    if region.ndim == 2:
        region = np.repeat(region[..., None], 3, -1)
    return np.clip(np.stack([map_coordinates(region[..., k].astype(float), [hy - by0, hx - bx0], order=1)
                             for k in range(3)], -1) / 255.0, 0, 1)


def run(bundle: str, outdir: str, markers: str | None = None, he_image: str | None = None,
        he_alignment: str | None = None, crop_um: float = 120, pick_um: float = 50, level: int = 2):
    log = setup_logging("morphology", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "morphology", log)
    M = out.metrics
    cfg = load_markers(markers)
    mo = Morphology(b)
    if not mo.channels:
        log.warning("no morphology_focus/ch*.ome.tif images in bundle; morphology checks skipped")
        M["skipped"] = True
        out.write_metrics()
        return
    log.info("morphology channels: %s", mo.names)
    dapi = next((n for n in mo.names if n.upper().startswith("DAPI")), mo.names[0])
    genes_q20 = b.transcripts().filter((pl.col("qv") >= 20) & pl.col("is_gene"))
    present = set(genes_q20.select(pl.col("feature_name").unique()).collect().to_series())

    # ---- 1. do marker transcripts land on bright pixels of their own stain?
    pairs = [(g, mo.resolve(c)) for g, c in cfg.get("image_pairs", []) if g in present and mo.resolve(c)]
    if pairs:
        with timed(log, f"transcript-on-stain enrichment for {len(pairs)} pairs"):
            mt = (genes_q20.filter(pl.col("feature_name").is_in([g for g, _ in pairs]))
                    .select("feature_name", "x_location", "y_location").collect().to_pandas())
            mt = mt.groupby("feature_name", group_keys=False).apply(lambda d: d.sample(min(len(d), 50_000), random_state=0))
            bg = (genes_q20.select("x_location", "y_location").filter(pl.int_range(pl.len()) % 300 == 0)
                    .collect().to_pandas())
            ds = 2 ** level

            def lookup(img, x, y):
                r = np.clip((y / b.pixel_size / ds).astype(int), 0, img.shape[0] - 1)
                c = np.clip((x / b.pixel_size / ds).astype(int), 0, img.shape[1] - 1)
                return img[r, c].astype(float)

            rows = []
            for gene, ch in pairs:
                img = mo.level_image(ch, level)
                bgv = lookup(img, bg.x_location.to_numpy(), bg.y_location.to_numpy())
                m = mt[mt.feature_name == gene]
                on = lookup(img, m.x_location.to_numpy(), m.y_location.to_numpy())
                bright = np.quantile(bgv, 0.95)
                rows.append({"gene": gene, "channel": ch, "n_transcripts": len(m),
                             "enrichment": on.mean() / max(bgv.mean(), 1e-9), "frac_on_bright": (on > bright).mean()})
                log.info("%s on %s: frac on top-5%% pixels %.3f", gene, ch, rows[-1]["frac_on_bright"])
                del img
            fid = pd.DataFrame(rows).sort_values("enrichment", ascending=False)
            out.table(fid, "qc_transcript_on_protein_image")
            M["weak_on_stain"] = fid[fid.frac_on_bright < 0.15][["gene", "channel"]].to_dict(orient="records")
            fig, ax = plt.subplots(figsize=(6, 0.3 * len(fid) + 1))
            ax.barh(fid.gene + " / " + fid.channel, fid.frac_on_bright)
            ax.axvline(0.05, color="k", ls="--", lw=1); ax.invert_yaxis()
            ax.set_xlabel("fraction of marker transcripts on top-5% stain pixels (random = 0.05)")
            ax.set_title("Marker transcripts on their matching stain")
            out.figure(fig, "qc_transcript_on_protein_image")
    else:
        log.info("no gene/channel image pairs available; on-stain enrichment skipped")

    # ---- 2. marker crops
    crops = []
    for c in cfg.get("morphology_crops", []):
        genes = [g for g in c["genes"] if g in present]
        if genes:
            crops.append((c["title"], genes, mo.resolve(c.get("channel")) or dapi))
    if not crops:
        out.write_metrics()
        return
    use_he = bool(he_image and he_alignment)
    align_inv = np.linalg.inv(np.loadtxt(he_alignment, delimiter=",", encoding="utf-8-sig")) if use_he else None
    bnd = pd.read_parquet(b.file("cell_boundaries.parquet")) if b.has("cell_boundaries.parquet") else None
    ncol = 3 if use_he else 2
    with timed(log, f"{len(crops)} marker crops (H&E column: {use_he})"):
        fig, axs = plt.subplots(len(crops), ncol, figsize=(4.4 * ncol, 4.3 * len(crops)), squeeze=False)
        for row, (title, genes, ch) in zip(axs, crops):
            hot = (genes_q20.filter(pl.col("feature_name").is_in(genes))
                     .group_by((pl.col("x_location") // pick_um).alias("bx"), (pl.col("y_location") // pick_um).alias("by"))
                     .agg(n=pl.len()).sort("n", descending=True).head(1).collect().row(0, named=True))
            x0 = hot["bx"] * pick_um + pick_um / 2 - crop_um / 2
            y0 = hot["by"] * pick_um + pick_um / 2 - crop_um / 2
            ext = (x0, x0 + crop_um, y0 + crop_um, y0)
            d = _norm(mo.crop(dapi, x0, y0, crop_um))
            p = _norm(mo.crop(ch, x0, y0, crop_um)) if ch != dapi else np.zeros_like(d)
            rgb = np.stack([p, 0.2 * p, np.maximum(d, p)], -1)
            tx = (genes_q20.filter(pl.col("feature_name").is_in(genes) & pl.col("x_location").is_between(x0, x0 + crop_um)
                                   & pl.col("y_location").is_between(y0, y0 + crop_um))
                    .select("x_location", "y_location").collect().to_pandas())
            axes = list(row)
            if use_he:
                axes[0].imshow(_he_crop(he_image, align_inv, b.pixel_size, x0, y0, crop_um), extent=ext)
                axes[0].set_title(f"{title}: H&E")
                axes = axes[1:]
            for ax in axes:
                ax.imshow(rgb, extent=ext)
                if bnd is not None:
                    sb = bnd[bnd.vertex_x.between(x0 - 20, x0 + crop_um + 20) & bnd.vertex_y.between(y0 - 20, y0 + crop_um + 20)]
                    for _, poly in sb.groupby("cell_id"):
                        ax.plot(poly.vertex_x, poly.vertex_y, color="white", lw=0.4, alpha=0.7)
            axes[0].set_title(f"{title}: {ch} (magenta) + DAPI" if ch != dapi else f"{title}: DAPI")
            axes[1].scatter(tx.x_location, tx.y_location, s=3, c="yellow", lw=0)
            axes[1].set_title(f"+ {'/'.join(genes)} transcripts (n={len(tx)})")
            for ax in row:
                ax.set_xlim(x0, x0 + crop_um); ax.set_ylim(y0 + crop_um, y0); ax.set_xticks([]); ax.set_yticks([])
            row[0].set_ylabel(f"x={x0:.0f}, y={y0:.0f} µm")
        fig.tight_layout()
        out.figure(fig, "qc_marker_morphology_crops", dpi=120)
    out.write_metrics()
