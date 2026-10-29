"""Unpack an XOA output bundle (.zip) into a folder, skipping very large files QC does not need."""
from __future__ import annotations

import zipfile
from pathlib import Path

from .common import find_bundle_root, setup_logging, timed

# Large files not read by any QC step (morphology z-stack, Xenium Explorer zarr archives, tarballs).
SKIP_DEFAULT = ("morphology.ome.tif", ".zarr.zip", ".tar.gz", "_mex/", "cell_feature_matrix/")


def run(input_path: str, outdir: str, extract_all: bool = False):
    out = Path(outdir)
    log = setup_logging("unpack", out)
    src = Path(input_path)
    target = out / "bundle"
    target.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        root = find_bundle_root(src)
        log.info("input is a folder; linking files from %s", root)
        for f in root.iterdir():
            dst = target / f.name
            if not dst.exists():
                dst.symlink_to(f.resolve())
        return
    if not zipfile.is_zipfile(src):
        raise ValueError(f"{src} is neither a folder nor a zip archive")
    with zipfile.ZipFile(src) as z, timed(log, f"extracting {src.name}"):
        members = z.infolist()
        keep = [m for m in members if extract_all or not any(s in m.filename for s in SKIP_DEFAULT)]
        skipped = len(members) - len(keep)
        size = sum(m.file_size for m in keep) / 1e9
        log.info("extracting %d of %d entries (%.2f GB); skipped %d large/unused entries",
                 len(keep), len(members), size, skipped)
        for m in keep:
            z.extract(m, target)
    root = find_bundle_root(target)
    if root != target:  # flatten a single nested folder so downstream steps see bundle/experiment.xenium
        log.info("bundle root nested at %s; flattening", root.relative_to(target))
        for f in root.iterdir():
            f.rename(target / f.name)
    log.info("bundle ready: %s", target)
