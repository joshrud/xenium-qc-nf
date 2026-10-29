"""Shared helpers: logging, bundle access, step outputs, plotting utilities."""
from __future__ import annotations

import contextlib
import json
import logging
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402
import yaml  # noqa: E402

from . import __version__  # noqa: E402

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
DEFAULT_MARKERS = Path(__file__).parent / "data" / "markers_default.yaml"


# --------------------------------------------------------------------------- logging
def setup_logging(step: str, outdir: Path, level: str = "INFO") -> logging.Logger:
    """Log to stderr and to <outdir>/logs/<step>.log. Python warnings are routed into the log."""
    logdir = Path(outdir) / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    log = logging.getLogger(f"xqc.{step}")
    log.setLevel(level)
    log.handlers.clear()
    fmt = logging.Formatter(LOG_FORMAT, "%Y-%m-%d %H:%M:%S")
    for h in (logging.StreamHandler(sys.stderr), logging.FileHandler(logdir / f"{step}.log", mode="w")):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.propagate = False
    logging.captureWarnings(True)
    wlog = logging.getLogger("py.warnings")
    wlog.handlers = log.handlers
    warnings.filterwarnings("ignore", category=FutureWarning)
    # gene and protein features share names (e.g. CD8A); each step separates them by feature type
    warnings.filterwarnings("ignore", message="Variable names are not unique")
    log.info("xenium_qc %s | step %s", __version__, step)
    return log


@contextlib.contextmanager
def timed(log: logging.Logger, label: str):
    t0 = time.perf_counter()
    log.info("%s ...", label)
    yield
    log.info("%s done (%.1f s)", label, time.perf_counter() - t0)


# --------------------------------------------------------------------------- bundle
def find_bundle_root(path: Path) -> Path:
    """Return the folder containing experiment.xenium (the bundle itself or a nested folder)."""
    path = Path(path)
    if (path / "experiment.xenium").exists():
        return path
    hits = sorted(path.glob("**/experiment.xenium"), key=lambda p: len(p.parts))
    if not hits:
        raise FileNotFoundError(f"No experiment.xenium found under {path}")
    return hits[0].parent


class Bundle:
    """Read access to an unpacked XOA output bundle."""

    def __init__(self, path: Path):
        self.path = find_bundle_root(Path(path))
        self.experiment = json.loads((self.path / "experiment.xenium").read_text())
        mpath = self.path / "metrics_summary.csv"
        self.metrics = pd.read_csv(mpath).iloc[0] if mpath.exists() else pd.Series(dtype=float)
        self.pixel_size = float(self.experiment.get("pixel_size", 0.2125))

    def file(self, name: str) -> Path:
        return self.path / name

    def has(self, name: str) -> bool:
        return (self.path / name).exists()

    def transcripts(self) -> pl.LazyFrame:
        return pl.scan_parquet(self.path / "transcripts.parquet")

    def cells(self) -> pd.DataFrame:
        return pd.read_parquet(self.path / "cells.parquet").set_index("cell_id")

    def morphology_channels(self) -> dict[str, Path]:
        """Channel name -> OME-TIFF file in morphology_focus/ (XOA >= 2.0 multi-file layout)."""
        import re

        import tifffile

        files = sorted((self.path / "morphology_focus").glob("ch*.ome.tif"))
        if not files:
            return {}
        with tifffile.TiffFile(files[0]) as tf:
            names = re.findall(r'Channel ID="[^"]+" Name="([^"]+)"', tf.ome_metadata or "")
        if len(names) == len(files):
            return dict(zip(names, files))
        # fall back to file names: ch0001_cd20.ome.tif -> CD20
        return {f.name.split("_", 1)[1].replace(".ome.tif", ""): f for f in files}

    def summary(self) -> dict:
        e = self.experiment
        return {k: e.get(k) for k in ["run_name", "region_name", "panel_name", "panel_num_targets_predesigned",
                                      "panel_num_targets_custom", "num_protein_targets", "analysis_sw_version",
                                      "instrument_sw_version", "preservation_method", "num_cells", "pixel_size",
                                      "segmentation_stain"]}


# --------------------------------------------------------------------------- outputs
class StepOutput:
    """Collects figures, tables and metrics for one step under a common output folder."""

    def __init__(self, outdir: Path, step: str, log: logging.Logger):
        self.outdir, self.step, self.log = Path(outdir), step, log
        for sub in ("figures", "tables", "metrics"):
            (self.outdir / sub).mkdir(parents=True, exist_ok=True)
        self.metrics: dict = {"step": step, "version": __version__}

    def figure(self, fig, name: str, dpi: int = 150):
        p = self.outdir / "figures" / f"{name}.png"
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        self.log.info("figure  %s", p.name)
        return p

    def table(self, df: pd.DataFrame, name: str, index: bool = False):
        p = self.outdir / "tables" / f"{name}.csv"
        df.to_csv(p, index=index)
        self.log.info("table   %s (%d rows)", p.name, len(df))
        return p

    def write_metrics(self):
        p = self.outdir / "metrics" / f"{self.step}.json"
        p.write_text(json.dumps(_jsonable(self.metrics), indent=2))
        self.log.info("metrics %s", p.name)


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, float) and np.isnan(o):
        return None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def load_markers(path: str | Path | None) -> dict:
    return yaml.safe_load(Path(path or DEFAULT_MARKERS).read_text())


# --------------------------------------------------------------------------- plotting
def to_grid(df: pd.DataFrame, col: str) -> np.ndarray:
    g = np.full((int(df.by.max()) + 1, int(df.bx.max()) + 1), np.nan)
    g[df.by.astype(int), df.bx.astype(int)] = df[col]
    return g


def show_grid(ax, grid, binsize, title, cmap="viridis", vmin=None, vmax=None):
    ny, nx = grid.shape
    im = ax.imshow(grid, cmap=cmap, vmin=vmin, vmax=vmax, extent=(0, nx * binsize, ny * binsize, 0),
                   interpolation="nearest")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, shrink=0.7)
    return im


def sample_rows(df: pd.DataFrame, n: int, seed: int = 0) -> pd.DataFrame:
    return df.sample(min(n, len(df)), random_state=seed)
