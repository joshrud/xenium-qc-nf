"""Build a small test XOA bundle (zip) by cropping a window out of a full Xenium output bundle.

The crop keeps the XOA file layout so every pipeline step runs on it:
transcripts.parquet, cells.parquet, cell_boundaries.parquet, nucleus_boundaries.parquet, cell_feature_matrix.h5,
metrics_summary.csv, experiment.xenium, gene_panel.json, and a few morphology_focus channels.
Morphology images are taken from pyramid level 2 (4x downsampled) and the bundle's pixel_size is updated to match,
so the test bundle stays small. Coordinates are shifted so the window starts at (0, 0); the H&E alignment matrix,
keypoints and annotation file (if given) are adjusted accordingly.

Usage:
  python tests/make_test_bundle.py --bundle /path/to/outs --out tests/data \
      [--he-alignment M.csv --keypoints keypoints.csv --annotation ann.geojson] [--x0 4000 --y0 1000 --w 2500 --h 2000]
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import polars as pl
import tifffile

CHANNELS = ["DAPI", "CD20", "CD3E", "CD31", "CD68", "PanCK"]
LEVEL = 2


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--bundle", required=True); a.add_argument("--out", required=True)
    a.add_argument("--he-alignment"); a.add_argument("--keypoints"); a.add_argument("--annotation")
    a.add_argument("--x0", type=float, default=4000); a.add_argument("--y0", type=float, default=1000)
    a.add_argument("--w", type=float, default=2500); a.add_argument("--h", type=float, default=2000)
    args = a.parse_args()
    src, out = Path(args.bundle), Path(args.out)
    x0, y0, x1, y1 = args.x0, args.y0, args.x0 + args.w, args.y0 + args.h
    name = "test_bundle_outs"
    dst = out / name
    if dst.exists():
        shutil.rmtree(dst)
    (dst / "morphology_focus").mkdir(parents=True)

    exp = json.loads((src / "experiment.xenium").read_text())
    px0 = exp["pixel_size"]
    px = px0 * 2 ** LEVEL
    exp["pixel_size"] = px
    exp["run_name"] = exp.get("run_name", "") + " (test crop)"
    (dst / "experiment.xenium").write_text(json.dumps(exp, indent=2))
    for f in ["metrics_summary.csv", "gene_panel.json"]:
        shutil.copy(src / f, dst / f)

    inwin = lambda xc, yc: pl.col(xc).is_between(x0, x1) & pl.col(yc).is_between(y0, y1)  # noqa: E731
    (pl.scan_parquet(src / "transcripts.parquet").filter(inwin("x_location", "y_location"))
       .with_columns(pl.col("x_location") - x0, pl.col("y_location") - y0)
       .sink_parquet(dst / "transcripts.parquet"))
    cells = pl.read_parquet(src / "cells.parquet").filter(inwin("x_centroid", "y_centroid"))
    ids = cells["cell_id"]
    cells.with_columns(pl.col("x_centroid") - x0, pl.col("y_centroid") - y0).write_parquet(dst / "cells.parquet")
    for f in ["cell_boundaries.parquet", "nucleus_boundaries.parquet"]:
        if (src / f).exists():
            (pl.read_parquet(src / f).filter(pl.col("cell_id").is_in(ids.to_list()))
               .with_columns(pl.col("vertex_x") - x0, pl.col("vertex_y") - y0).write_parquet(dst / f))

    # 10x HDF5 matrix: keep the selected barcodes (CSC: columns = cells)
    with h5py.File(src / "cell_feature_matrix.h5") as fi, h5py.File(dst / "cell_feature_matrix.h5", "w") as fo:
        for k, v in fi.attrs.items():
            fo.attrs[k] = v
        m = fi["matrix"]
        bc = m["barcodes"][:].astype(str)
        wanted = set(ids.to_list())  # set lookup: np.isin on string arrays is quadratic
        keep = np.flatnonzero([b in wanted for b in bc])
        indptr, indices, data = m["indptr"][:], m["indices"][:], m["data"][:]
        new_ind, new_dat, new_ptr = [], [], [0]
        for j in keep:
            s, e = indptr[j], indptr[j + 1]
            new_ind.append(indices[s:e]); new_dat.append(data[s:e]); new_ptr.append(new_ptr[-1] + e - s)
        g = fo.create_group("matrix")
        g.create_dataset("barcodes", data=m["barcodes"][:][keep])
        g.create_dataset("data", data=np.concatenate(new_dat), compression="gzip")
        g.create_dataset("indices", data=np.concatenate(new_ind), compression="gzip")
        g.create_dataset("indptr", data=np.array(new_ptr, dtype=indptr.dtype))
        g.create_dataset("shape", data=np.array([m["shape"][0], len(keep)], dtype=m["shape"].dtype))
        fi.copy(m["features"], g, "features")

    # morphology: crop pyramid level 2 of selected channels into single-channel pyramidal OME-TIFFs
    files = sorted((src / "morphology_focus").glob("ch*.ome.tif"))
    with tifffile.TiffFile(files[0]) as tf:
        names = re.findall(r'Channel ID="[^"]+" Name="([^"]+)"', tf.ome_metadata)
    r0, r1, c0, c1 = int(y0 / px), int(y1 / px), int(x0 / px), int(x1 / px)
    for k, ch in enumerate([c for c in CHANNELS if c in names]):
        i = names.index(ch)
        with tifffile.TiffFile(files[i]) as tf:
            img = tf.series[0].levels[LEVEL].pages[i].asarray()[r0:r1, c0:c1]
        meta = {"axes": "YX", "Channel": {"Name": [ch]}, "PhysicalSizeX": px, "PhysicalSizeY": px}
        with tifffile.TiffWriter(dst / "morphology_focus" / f"ch{k:04d}_{ch.lower()}.ome.tif", bigtiff=True, ome=True) as tw:
            tw.write(img, subifds=2, tile=(256, 256), compression="zlib", metadata=meta)
            for lv in (1, 2):
                tw.write(img[:: 2 ** lv, :: 2 ** lv], subfiletype=1, tile=(256, 256), compression="zlib")

    # alignment matrix (H&E px -> DAPI px): shift by the window and rescale to the new pixel size
    extra = out / "test_inputs"
    extra.mkdir(exist_ok=True)
    T = np.array([[px0 / px, 0, -x0 / px], [0, px0 / px, -y0 / px], [0, 0, 1]])
    if args.he_alignment:
        M = np.loadtxt(args.he_alignment, delimiter=",", encoding="utf-8-sig")
        np.savetxt(extra / "he_imagealignment.csv", T @ M, delimiter=",")
    if args.keypoints:
        kp = pd.read_csv(args.keypoints)
        kp["fixedX"] = (kp.fixedX * px0 - x0) / px; kp["fixedY"] = (kp.fixedY * px0 - y0) / px
        kp.to_csv(extra / "keypoints.csv", index=False)
    if args.annotation:
        shutil.copy(args.annotation, extra / "annotation.geojson")

    zpath = out / f"{name}.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_STORED) as z:
        for f in dst.rglob("*"):
            if f.is_file():
                z.write(f, f.relative_to(out))
    shutil.rmtree(dst)
    print(f"wrote {zpath} ({zpath.stat().st_size / 1e6:.1f} MB), {len(ids)} cells; extra inputs in {extra}")


if __name__ == "__main__":
    main()
