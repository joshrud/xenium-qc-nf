"""Run-level metrics from metrics_summary.csv graded against 10x warning/error thresholds."""
from __future__ import annotations

import pandas as pd

from .common import Bundle, StepOutput, setup_logging

# metric: (direction, warning, error, benchmark note). Thresholds follow 10x's analysis-summary alerts.
THRESHOLDS = {
    "fraction_transcripts_decoded_q20": ("min", 0.60, 0.50, "mean 0.81 (0.72-0.91), Salas et al. 2025"),
    "fraction_predesigned_transcripts_decoded_q20": ("min", 0.60, 0.50, ""),
    "fraction_custom_transcripts_decoded_q20": ("min", 0.60, 0.50, ""),
    "nuclear_transcripts_per_100um2": ("min", 10, 1, "panel/tissue dependent"),
    "negative_control_probe_counts_per_control_per_cell": ("max", 0.025, 0.05, ""),
    "genomic_control_probe_counts_per_control_per_cell": ("max", 0.01, 0.02, "Xenium Prime only"),
    "fraction_empty_cells": ("max", None, 0.10, ""),
    "segmented_cell_stain_frac": ("min", 0.70, 0.50, ""),
    "fraction_transcripts_assigned": ("min", 0.50, None, "mean 0.768, Salas et al. 2025"),
}
INFO = ["num_cells_detected", "median_transcripts_per_cell", "median_genes_per_cell",
        "decoded_transcripts_per_100um2", "cells_per_100um2", "adjusted_negative_control_probe_rate",
        "adjusted_negative_control_codeword_rate", "estimated_number_of_false_positive_transcripts_per_cell",
        "thickness_transcripts_high_quality", "region_area", "total_high_quality_decoded_transcripts"]


def grade(value, direction, warn, err):
    if value is None or pd.isna(value):
        return "n/a"
    worse = (lambda v, t: v < t) if direction == "min" else (lambda v, t: v > t)
    if err is not None and worse(value, err):
        return "ERROR"
    if warn is not None and worse(value, warn):
        return "WARNING"
    return "pass"


def run(bundle: str, outdir: str):
    log = setup_logging("run_metrics", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "run_metrics", log)
    m = b.metrics
    run_qc = pd.DataFrame([{"metric": k, "value": m.get(k), "warning": w, "error": e,
                            "status": grade(m.get(k), d, w, e), "benchmark": note}
                           for k, (d, w, e, note) in THRESHOLDS.items()])
    run_info = pd.DataFrame([{"metric": k, "value": m.get(k)} for k in INFO if k in m.index])
    for _, r in run_qc.iterrows():
        log.info("%-55s %-12s %s", r.metric, r.value, r.status)
    out.table(run_qc, "qc_run_metrics")
    out.table(run_info, "qc_run_info")
    out.metrics.update({
        "sample": b.summary(),
        "n_pass": int((run_qc.status == "pass").sum()), "n_warning": int((run_qc.status == "WARNING").sum()),
        "n_error": int((run_qc.status == "ERROR").sum()),
        "fraction_transcripts_decoded_q20": m.get("fraction_transcripts_decoded_q20"),
        "fraction_transcripts_assigned": m.get("fraction_transcripts_assigned"),
        "num_cells_detected": m.get("num_cells_detected"),
        "thickness_transcripts_high_quality": m.get("thickness_transcripts_high_quality"),
        "region_area": m.get("region_area"),
    })
    out.write_metrics()
