"""Retain/drop/flag decisions and a single self-contained HTML QC report."""
from __future__ import annotations

import base64
import datetime as dt
import html
import io
import json
import platform
from pathlib import Path

import pandas as pd
from PIL import Image

from . import __version__
from .common import StepOutput, setup_logging

SECTIONS = [
    ("run_metrics", "Run metrics", "XOA run-level metrics from <code>metrics_summary.csv</code>, graded against the warning/error "
     "thresholds used in the 10x analysis summary.", ["qc_run_metrics", "qc_run_info"], []),
    ("transcripts", "Transcript quality", "Decoding quality (QV), genes against negative-control background, the z-distribution of "
     "transcripts, per-FOV consistency, transcript density, stitching seams, negative-control hotspots, section tilt and "
     "local thickness, and transcript assignment to cells.",
     ["qc_codeword_categories", "qc_fov"],
     ["qc_qv_distribution", "qc_gene_vs_negative_controls", "qc_z_distribution", "qc_fov_maps", "qc_transcript_density",
      "qc_fov_seam_profile", "qc_negative_control_hotspots", "qc_z_maps", "qc_assignment"]),
    ("cells", "Cell segmentation and filtering", "Segmentation method, per-cell count/gene/area distributions and their spatial "
     "pattern, hard filters (cells removed) and soft flags, and segmentation spillover: how much of a lineage marker leaks "
     "into cells of a mutually exclusive lineage (whole cell vs nucleus-only counts; positive control sets the ceiling).",
     ["qc_segmentation_methods", "qc_cell_filter_summary", "qc_segmentation_spillover"],
     ["qc_cell_distributions", "qc_counts_vs_area", "qc_spatial_maps", "qc_removed_cells"]),
    ("protein", "Protein (Xenium Protein runs)", "Per-protein detection, per-cell RNA-protein concordance for matched markers, and "
     "co-occurrence of mutually exclusive lineage proteins.",
     ["qc_protein_summary", "qc_rna_protein_concordance", "qc_protein_exclusivity"], ["qc_rna_protein_concordance"]),
    ("morphology", "Markers against morphology", "Whether marker transcripts fall on bright pixels of their matching stain "
     "(random expectation 0.05), and crops of the densest 50 µm bin per marker with stain, cell boundaries, transcripts and, "
     "when an aligned H&E is supplied, the H&E.", ["qc_transcript_on_protein_image"],
     ["qc_transcript_on_protein_image", "qc_marker_morphology_crops"]),
    ("alignment", "H&E alignment", "Per-keypoint fit residual and leave-one-out (LOO) error of a least-squares similarity "
     "transform (the model Xenium Explorer fits), comparison with a reference transform across the capture area, and panel "
     "genes ranked as landmark guides (abundant, spatially sharp, spread across the tissue).",
     ["alignment_transforms", "alignment_keypoint_errors", "keypoint_gene_ranking"],
     ["alignment_keypoint_errors", "alignment_transform_difference"]),
    ("regions", "Regions and cell types", "Annotation polygons converted to Xenium µm and cells assigned by centroid; Leiden "
     "clusters labelled by the highest-scoring marker set; composition and median transcript depth per region. Depth "
     "differences between regions must be matched (e.g. by downsampling) before differential expression.",
     ["annotation_regions", "cluster_labels", "cell_type_composition", "depth_by_region_cell_type"],
     ["regions_selection", "cluster_markers", "cell_type_composition", "depth_by_region"]),
]


def build_decisions(m: dict) -> pd.DataFrame:
    """Rule-based retain/drop/flag table from step metrics. Each row states the rule applied."""
    rm, tr, ce, pr, mo, al, rg = (m.get(k, {}) for k in ["run_metrics", "transcripts", "cells", "protein", "morphology", "alignment", "regions"])
    rows = []

    def add(level, check, result, decision, rule):
        rows.append({"level": level, "check": check, "result": result, "decision": decision, "rule": rule})

    if rm:
        add("Run", "10x metrics", f"{rm['n_pass']} pass, {rm['n_warning']} warn, {rm['n_error']} error",
            "Retain" if rm["n_error"] == 0 else "Flag", "Flag if any metric is past its 10x error threshold")
    if tr:
        add("Transcripts", "QV < 20", f"{1 - tr['gene_frac_q20']:.1%} of gene spots", "Drop", "Always excluded (XOA excludes them from the matrix)")
        add("Transcripts", "control codewords/probes", f"{tr['n_control_q20']:,} Q20 spots", "Drop", "Background estimate only, never features")
        low = tr.get("genes_below_background", [])
        add("Genes", "below negative-control background", ", ".join(low) or "none", "Flag" if low else "Retain",
            "Gene Q20 count <= 99th percentile of negative controls")
        add("FOVs", "low Q20 / high controls", f"{tr['fov_low_q20']} / {tr['fov_high_control']} of {tr['n_fov']}",
            "Flag" if tr["fov_low_q20"] + tr["fov_high_control"] else "Retain", "3 MAD from the FOV median")
        add("FOVs", "low transcript density", f"{tr['fov_low_density']} FOVs", "Review" if tr["fov_low_density"] else "Retain",
            "3 MAD below median; often partly empty edge FOVs")
        add("Space", "stitching seams", f"max deviation {tr['seam_max_deviation']:.1%} within 8 µm of a border",
            "Retain" if tr["seam_max_deviation"] < 0.1 else "Flag", "Flag if density deviates >= 10% at FOV borders")
        add("Space", "negative-control hotspots", f"{tr['control_hotspots']} of {tr['control_bins']} bins",
            "Retain" if tr["control_hotspots"] == 0 else "Flag", "Poisson test, Bonferroni, >2x expected, >1,000 spots")
        add("Space", "thick / folded tissue", f"{tr['thick_bin_fraction']:.2%} of tissue bins",
            "Flag" if tr["thick_bin_fraction"] > 0.005 else "Retain", "z IQR > median + 5 MAD; flag if > 0.5% of tissue")
    if ce:
        f = ce.get("filters", {})
        for key, label, dec in [("fail_low_counts", "< min transcripts", "Drop"), ("fail_no_nucleus", "no nucleus", "Drop"),
                                ("flag_area_outlier", "area outlier (MAD)", "Review"), ("flag_counts_outlier", "count outlier (MAD)", "Review")]:
            if key in f:
                add("Cells", label, f"{int(f[key]['n_cells']):,} ({f[key]['fraction']:.2%})", dec,
                    "Hard filter" if dec == "Drop" else "Soft flag; inspect if spatially clustered")
        if "nuclear_expansion_fraction" in ce:
            add("Cells", "nuclear-expansion segmentation", f"{ce['nuclear_expansion_fraction']:.1%}",
                "Flag" if ce["nuclear_expansion_fraction"] > 0.05 else "Retain", "Flag if > 5% of cells")
        sp = [s for s in ce.get("spillover", []) if s["type"] == "should be exclusive" and (s["leak_whole_cell"] or 0) >= 0.15]
        add("Cells", "boundary spillover", ", ".join(f"{s['A_high']}->{s['B']} {s['leak_whole_cell']:.2f}" for s in sp) or "none >= 0.15",
            "Flag" if sp else "Retain", "Leak ratio >= 0.15 for a mutually exclusive pair")
    if pr and not pr.get("skipped"):
        w = pr.get("weak_concordance", [])
        add("Protein", "weak RNA-protein concordance", ", ".join(w) or "none", "Flag" if w else "Retain", "Spearman ρ < 0.15")
    if mo and not mo.get("skipped"):
        w = mo.get("weak_on_stain", [])
        add("Morphology", "marker transcripts off their stain", ", ".join(f"{x['gene']}/{x['channel']}" for x in w) or "none",
            "Review" if w else "Retain", "< 15% of spots on top-5% stain pixels (broad stains cannot pass this test)")
    if al.get("n_keypoints"):
        add("Alignment", "keypoint LOO error", f"median {al['loo_median_um']:.2f} µm, max {al['loo_max_um']:.2f} µm (#{al['worst_keypoint']})",
            "Retain" if al["loo_max_um"] < 5 else "Review", "Review keypoints with LOO error >= 5 µm")
    if al.get("vs_reference_max_um") is not None:
        add("Alignment", "vs reference transform", f"max {al['vs_reference_max_um']:.2f} µm", "Retain" if al["vs_reference_max_um"] < 5 else "Review",
            "Review if the transforms differ by >= 5 µm anywhere")
    if rg.get("depth_ratio_max_min"):
        add("Regions", "transcript depth imbalance", f"max/min median {rg['depth_ratio_max_min']:.2f}x",
            "Flag" if rg["depth_ratio_max_min"] > 1.2 else "Retain", "Flag if > 1.2x; match depth before comparing regions")
    return pd.DataFrame(rows)


def _img_b64(path: Path, max_w: int = 1500) -> str:
    im = Image.open(path).convert("RGB")
    if im.width > max_w:
        im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=85, optimize=True)
    return base64.b64encode(buf.getvalue()).decode()


def _table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    more = f"<p class='note'>Showing {max_rows} of {len(df)} rows; full table in tables/.</p>" if len(df) > max_rows else ""
    return df.head(max_rows).to_html(index=False, border=0, classes="tbl", float_format=lambda v: f"{v:,.4g}", na_rep="") + more


CSS = """
:root{--fg:#1d1d1f;--muted:#5f6368;--bg:#fff;--card:#f6f8fa;--line:#d0d7de;--accent:#1f4e79;--pass:#1a7f37;--warn:#9a6700;--err:#cf222e}
@media (prefers-color-scheme: dark){:root{--fg:#e6edf3;--muted:#9da7b3;--bg:#0d1117;--card:#161b22;--line:#30363d;--accent:#79c0ff}}
body{font:15px/1.5 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg);margin:0}
main{max-width:1180px;margin:auto;padding:24px 16px 64px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:20px;margin:36px 0 8px;color:var(--accent);border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:15px;margin:18px 0 6px}.muted,.note{color:var(--muted)}.note{font-size:12px}
nav{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);padding:8px 0;margin-bottom:8px;font-size:13px;overflow-x:auto;white-space:nowrap}
nav a{margin-right:14px;color:var(--accent);text-decoration:none}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(170px,1fr));gap:10px;margin:14px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 12px}.card b{font-size:20px;display:block}
.tbl{border-collapse:collapse;font-size:12.5px;margin:6px 0 14px;display:block;overflow-x:auto;max-width:100%}
.tbl th,.tbl td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;white-space:nowrap}.tbl th{background:var(--card)}
.Retain,.pass{color:var(--pass);font-weight:600}.Flag,.WARNING,.Review{color:var(--warn);font-weight:600}.Drop,.ERROR{color:var(--err);font-weight:600}
figure{margin:12px 0 20px}figure img{max-width:100%;border:1px solid var(--line);border-radius:6px;background:#fff}
figcaption{font-size:12px;color:var(--muted)}code{font-size:12.5px}
"""


def run(step_dirs: list[str], outdir: str, sample: str, multiqc_report: str | None = None, params_json: str | None = None):
    log = setup_logging("report", outdir)
    out = StepOutput(outdir, "report", log)
    dirs = [Path(d) for d in step_dirs]
    metrics, figs, tabs = {}, {}, {}
    for d in dirs:
        for p in (d / "metrics").glob("*.json") if (d / "metrics").exists() else []:
            metrics[p.stem] = json.loads(p.read_text())
        for p in (d / "figures").glob("*.png") if (d / "figures").exists() else []:
            figs[p.stem] = p
        for p in (d / "tables").glob("*.csv") if (d / "tables").exists() else []:
            tabs[p.stem] = p
    log.info("collected %d metrics files, %d figures, %d tables from %d step folders", len(metrics), len(figs), len(tabs), len(dirs))

    dec = build_decisions(metrics)
    out.table(dec, "qc_decisions")
    rm, tr, ce, al = (metrics.get(k, {}) for k in ["run_metrics", "transcripts", "cells", "alignment"])
    info = rm.get("sample", {})

    def card(label, value):
        return f"<div class='card'><span class='muted'>{html.escape(label)}</span><b>{html.escape(str(value))}</b></div>"

    cards = [card("Cells detected", f"{ce.get('n_cells', info.get('num_cells', '–')):,}" if ce else info.get("num_cells", "–"))]
    if ce:
        cards.append(card("Cells passing QC", f"{ce['n_pass']:,} ({ce['n_pass'] / ce['n_cells']:.1%})"))
        cards.append(card("Median transcripts / cell", f"{ce['median_counts']:.0f}"))
    if tr:
        cards += [card("Gene transcripts ≥ Q20", f"{tr['gene_frac_q20']:.1%}"),
                  card("Genes above background", f"{tr['n_genes_above_background']} / {tr['n_genes']}"),
                  card("Tissue area", f"{tr['tissue_area_mm2']:.1f} mm²")]
    if rm:
        cards.append(card("10x metrics", f"{rm['n_pass']} pass · {rm['n_warning']} warn · {rm['n_error']} error"))
    if al.get("n_keypoints"):
        cards.append(card("Alignment LOO error", f"{al['loo_median_um']:.2f} µm (max {al['loo_max_um']:.2f})"))

    dec_html = dec.to_html(index=False, border=0, classes="tbl", escape=True)
    for d in ["Retain", "Drop", "Flag", "Review"]:
        dec_html = dec_html.replace(f"<td>{d}</td>", f"<td class='{d}'>{d}</td>")
    nav, body = ["<a href='#summary'>Summary</a>"], []
    for key, title, desc, tnames, fnames in SECTIONS:
        if key not in metrics or metrics[key].get("skipped"):
            continue
        nav.append(f"<a href='#{key}'>{title}</a>")
        parts = [f"<h2 id='{key}'>{title}</h2><p>{desc}</p>"]
        for f in fnames:
            if f in figs:
                parts.append(f"<figure><img alt='{f}' src='data:image/jpeg;base64,{_img_b64(figs[f])}'>"
                             f"<figcaption>figures/{f}.png</figcaption></figure>")
        for t in tnames:
            if t in tabs:
                df = pd.read_csv(tabs[t])
                if t == "qc_run_metrics":
                    h = _table_html(df)
                    for s in ["pass", "WARNING", "ERROR"]:
                        h = h.replace(f"<td>{s}</td>", f"<td class='{s}'>{s}</td>")
                    parts.append(f"<h3>{t}</h3>{h}")
                else:
                    parts.append(f"<h3>{t}</h3>{_table_html(df)}")
        body.append("\n".join(parts))

    meta = {"Sample": sample, "Run": info.get("run_name"), "Region": info.get("region_name"), "Panel": info.get("panel_name"),
            "XOA version": info.get("analysis_sw_version"), "Preservation": info.get("preservation_method")}
    meta_html = " · ".join(f"<b>{k}:</b> {html.escape(str(v))}" for k, v in meta.items() if v)
    links = []
    if multiqc_report:
        links.append(f"<li><a href='{html.escape(multiqc_report)}'>MultiQC report (core Xenium module + Xenium Extra plugin)</a></li>")
    params = json.loads(Path(params_json).read_text()) if params_json and Path(params_json).exists() else {}
    soft = {"xenium_qc": __version__, "python": platform.python_version()}
    for mod in ["scanpy", "polars", "pandas", "numpy", "tifffile", "shapely"]:
        try:
            soft[mod] = __import__(mod).__version__
        except Exception:  # noqa: BLE001
            pass
    nav.append("<a href='#provenance'>Provenance</a>")
    doc = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>Xenium QC · {html.escape(sample)}</title><style>{CSS}</style></head><body><main>
<h1>Xenium QC report</h1><p class='muted'>{meta_html}</p>
<nav>{''.join(nav)}</nav>
<h2 id='summary'>Summary</h2><div class='cards'>{''.join(cards)}</div>
<h3>Retain, drop or flag</h3><p class='note'>Rule-based from the metrics below; "Review" means inspect before deciding.</p>{dec_html}
{''.join(body)}
<h2 id='provenance'>Provenance</h2>
<p class='note'>Generated {dt.datetime.now():%Y-%m-%d %H:%M} by xenium-qc-nf / xenium_qc {__version__}. Per-step logs are in <code>logs/</code>;
Nextflow trace, timeline and report in <code>pipeline_info/</code>.</p><ul>{''.join(links)}</ul>
<h3>Software</h3>{pd.DataFrame(list(soft.items()), columns=['software', 'version']).to_html(index=False, border=0, classes='tbl')}
{('<h3>Parameters</h3>' + pd.DataFrame(list(params.items()), columns=['parameter', 'value']).to_html(index=False, border=0, classes='tbl')) if params else ''}
</main></body></html>"""
    p = Path(outdir) / f"{sample}_xenium_qc_report.html"
    p.write_text(doc)
    log.info("report written: %s (%.1f MB)", p.name, p.stat().st_size / 1e6)
    out.metrics["decisions"] = dec.to_dict(orient="records")
    out.write_metrics()
