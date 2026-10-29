"""Command-line entry point: `xqc <step> ...`. Each step writes figures/, tables/, metrics/ and logs/ under --outdir."""
from __future__ import annotations

import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(prog="xqc", description="Xenium (XOA) QC and alignment checks")
    sub = p.add_subparsers(dest="step", required=True)

    s = sub.add_parser("unpack", help="unzip an XOA bundle (or link a folder) into OUTDIR/bundle")
    s.add_argument("--input", required=True); s.add_argument("--outdir", required=True)
    s.add_argument("--all", action="store_true", help="extract every file (default skips morphology.ome.tif and zarr archives)")

    def common(sp, markers=False):
        sp.add_argument("--bundle", required=True, help="unpacked XOA output folder")
        sp.add_argument("--outdir", required=True)
        if markers:
            sp.add_argument("--markers", help="marker YAML (default: built-in)")

    common(sub.add_parser("run-metrics", help="grade metrics_summary.csv against 10x thresholds"))
    s = sub.add_parser("transcripts", help="transcript-level QC"); common(s)
    s.add_argument("--bin-um", type=int, default=20); s.add_argument("--seam-window", type=int, default=31)

    s = sub.add_parser("cells", help="cell-level QC, filters, spillover, normalization"); common(s, markers=True)
    s.add_argument("--min-counts", type=int, default=10); s.add_argument("--min-cells-per-gene", type=int, default=5)
    s.add_argument("--mad-k", type=float, default=5.0); s.add_argument("--target-sum", type=float, default=100.0)

    s = sub.add_parser("protein", help="protein QC (Xenium Protein runs)"); common(s, markers=True)
    s.add_argument("--cells-qc", required=True, help="cells_qc.parquet from the cells step")

    s = sub.add_parser("morphology", help="marker transcripts against morphology images"); common(s, markers=True)
    s.add_argument("--he-image"); s.add_argument("--he-alignment")

    s = sub.add_parser("alignment", help="H&E alignment accuracy and landmark genes"); common(s)
    s.add_argument("--keypoints"); s.add_argument("--he-alignment"); s.add_argument("--reference-alignment")

    s = sub.add_parser("regions", help="annotation regions, cell types, composition, depth"); common(s, markers=True)
    s.add_argument("--filtered-h5ad", required=True); s.add_argument("--annotation"); s.add_argument("--he-alignment")
    s.add_argument("--annotation-units", choices=["he_px", "um"], default="he_px")
    s.add_argument("--leiden-resolution", type=float, default=1.0)

    s = sub.add_parser("report", help="decisions table and HTML report")
    s.add_argument("--step-dirs", nargs="+", required=True); s.add_argument("--outdir", required=True)
    s.add_argument("--sample", required=True); s.add_argument("--multiqc-report"); s.add_argument("--params-json")

    a = p.parse_args(argv)
    if a.step == "unpack":
        from .unpack import run
        run(a.input, a.outdir, a.all)
    elif a.step == "run-metrics":
        from .run_metrics import run
        run(a.bundle, a.outdir)
    elif a.step == "transcripts":
        from .transcripts import run
        run(a.bundle, a.outdir, bin_um=a.bin_um, seam_window=a.seam_window)
    elif a.step == "cells":
        from .cells import run
        run(a.bundle, a.outdir, a.markers, a.min_counts, a.min_cells_per_gene, a.mad_k, a.target_sum)
    elif a.step == "protein":
        from .protein import run
        run(a.bundle, a.outdir, a.cells_qc, a.markers)
    elif a.step == "morphology":
        from .morphology import run
        run(a.bundle, a.outdir, a.markers, a.he_image, a.he_alignment)
    elif a.step == "alignment":
        from .alignment import run
        run(a.bundle, a.outdir, a.keypoints, a.he_alignment, a.reference_alignment)
    elif a.step == "regions":
        from .regions import run
        run(a.bundle, a.outdir, a.filtered_h5ad, a.annotation, a.he_alignment, a.annotation_units, a.markers,
            a.leiden_resolution)
    elif a.step == "report":
        from .report import run
        run(a.step_dirs, a.outdir, a.sample, a.multiqc_report, a.params_json)


if __name__ == "__main__":
    sys.exit(main())
