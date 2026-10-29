# xenium-qc-nf

QC and H&E-alignment checks for **10x Xenium Onboard Analysis (XOA)** output bundles, as a Nextflow pipeline.
Point it at the `outs` **.zip** (or the unpacked folder) and it produces every QC figure and table, per-step logs,
and a single self-contained **HTML QC report**. It stops at the hand-off to differential expression: QC-filtered,
normalized cells (`filtered.h5ad`) with region and cell-type labels.

It is a companion to [nf-core/spatialaxe](https://github.com/nf-core/spatialaxe) (re-segmentation and processing of
Xenium data). It uses the same samplesheet columns (`sample,bundle,image`) and runs the same MultiQC-based QC that
spatialaxe's `--mode qc` runs ([MultiQC](https://multiqc.info) + the
[Xenium Extra plugin](https://github.com/MultiQC/xenium-extra)), then adds spatial, segmentation, morphology,
alignment and region checks that spatialaxe does not cover.

## Quick start

```bash
git clone https://github.com/joshrud/xenium-qc-nf && cd xenium-qc-nf

# Python environment (or use -profile conda / -profile docker below)
python3.11 -m venv .venv && .venv/bin/pip install -e ".[multiqc]"

# One bundle, straight from the XOA zip
nextflow run main.nf -profile local --venv .venv \
    --bundle /data/Xenium_..._outs.zip --sample my_sample --outdir results

# With the post-run H&E, your Xenium Explorer alignment export and a pathology annotation
nextflow run main.nf -profile local --venv .venv \
    --bundle /data/Xenium_..._outs.zip --sample my_sample \
    --he_image /data/he_image.ome.tif \
    --he_alignment /data/alignment_files/matrix.csv --keypoints /data/alignment_files/keypoints.csv \
    --reference_alignment /data/Xenium_..._he_imagealignment.csv \
    --annotation /data/annotation.geojson

# Several samples
nextflow run main.nf -profile local --venv .venv --input samplesheet.csv --outdir results
```

Open `results/<sample>/report/<sample>_xenium_qc_report.html`.

**Test run** (a small cropped bundle, zip input, all optional inputs):

```bash
python tests/make_test_bundle.py --bundle /data/Xenium_..._outs --out tests/data \
    --he-alignment /data/he_imagealignment.csv --keypoints /data/alignment_files/keypoints.csv \
    --annotation /data/annotation.geojson
nextflow run main.nf -profile test,local --venv .venv
```

### Samplesheet

Only `sample` and `bundle` are required; leave other cells empty when not available.

| column | meaning |
|---|---|
| `sample` | sample id (output folder name) |
| `bundle` | XOA output `.zip` or unpacked `outs/` folder |
| `image` | post-run H&E OME-TIFF (adds H&E panels to the marker crops) |
| `alignment` | your 3×3 H&E→DAPI pixel matrix (Xenium Explorer *Download Image Alignment File*, `matrix.csv`) |
| `keypoints` | Xenium Explorer `keypoints.csv` (fixedX/fixedY = DAPI px, alignmentX/alignmentY = H&E px) |
| `annotation` | pathology GeoJSON drawn on the H&E (or in µm with `--annotation_units um`) |
| `reference_alignment` | a second matrix to compare against, e.g. 10x's supplemental `*_he_imagealignment.csv` |

The H&E crops and the annotation→Xenium transform use `alignment` if given, otherwise `reference_alignment`.

### Profiles

| profile | environment |
|---|---|
| `local` | uses `xqc` on `PATH`; add `--venv <path>` to use a virtualenv |
| `conda` | builds `environment.yml` |
| `docker` | `docker build -t xenium-qc-nf:0.1.0 .` first (Dockerfile in the repo root) |
| `test` | cropped test bundle (`tests/make_test_bundle.py`) |

## What it checks

| Step (`results/<sample>/<step>/`) | Checks | Main outputs |
|---|---|---|
| `run_metrics` | `metrics_summary.csv` graded against the 10x analysis-summary warning/error thresholds | `qc_run_metrics`, `qc_run_info` |
| `transcripts` | QV distribution; per-gene Q20 counts vs negative-control background; z-distribution (section thickness); per-FOV Q20/density/control outliers (3 MAD); transcript density map; stitching-seam profile (4 µm strips, 31-strip running median); negative-control hotspots (Poisson, Bonferroni); z median/IQR maps (tilt, folds); transcript assignment map and distance to nucleus | `qc_qv_distribution`, `qc_gene_vs_negative_controls`, `qc_z_distribution`, `qc_fov_maps`, `qc_transcript_density`, `qc_fov_seam_profile`, `qc_negative_control_hotspots`, `qc_z_maps`, `qc_assignment`; tables `qc_codeword_categories`, `qc_gene_detection`, `qc_fov`, `qc_negative_control_bins` |
| `cells` | segmentation methods; count/gene/area distributions and maps; hard filters (< 10 counts, no nucleus) and soft flags (5-MAD area/count outliers, any control counts); gene filter (≥ 5 cells); segmentation spillover leak ratios (whole cell vs nucleus-only, with a positive control); library-size normalization to 100 + log1p | `qc_cell_distributions`, `qc_counts_vs_area`, `qc_spatial_maps`, `qc_removed_cells`; tables `qc_segmentation_methods`, `qc_cell_filter_summary`, `qc_segmentation_spillover`; `cells_qc.parquet`, `filtered.h5ad` |
| `protein` | (Xenium Protein runs) per-protein detection, RNA–protein concordance, lineage-exclusive co-detection | `qc_protein_summary`, `qc_rna_protein_concordance`, `qc_protein_exclusivity` |
| `morphology` | do marker transcripts land on bright pixels of their matching stain; crops of the densest bin per marker with stain, DAPI, cell outlines, transcripts and (optional) aligned H&E | `qc_transcript_on_protein_image`, `qc_marker_morphology_crops` |
| `alignment` | least-squares similarity fit to the keypoints (the model Xenium Explorer fits), per-keypoint residual and leave-one-out error mapped over the tissue, transform parameters, difference from a reference transform across the capture area; panel genes ranked as landmark guides | `alignment_keypoint_errors`, `alignment_transform_difference`; tables `alignment_keypoint_errors`, `alignment_transforms`, `keypoint_gene_ranking` |
| `regions` | annotation polygons → Xenium µm; cells assigned by centroid; Leiden clusters labelled by marker score; composition and transcript depth per region (depth must be matched before comparing regions) | `regions_selection`, `cluster_markers`, `cell_type_composition`, `depth_by_region`; `annotations_xenium_um.geojson`, `tables/cell_labels.csv`, `cell_groups/*.csv` (Xenium Explorer cell groups) |
| `multiqc` | MultiQC core Xenium module + Xenium Extra plugin | `multiqc_report.html` |
| `report` | rule-based retain/drop/flag/review table and the HTML report | `<sample>_xenium_qc_report.html`, `tables/qc_decisions.csv` |

Marker genes, spillover pairs, protein–gene pairs and crop definitions are in
[`xenium_qc/data/markers_default.yaml`](xenium_qc/data/markers_default.yaml). Genes or channels missing from a panel are
skipped. Copy and edit it for another tissue, then pass `--markers my_markers.yaml`.

Large files that no step reads (`morphology.ome.tif`, Xenium Explorer `*.zarr.zip`, tarballs) are not extracted from
the zip unless `--unzip_all` is set.

## Logging

- Each step writes `results/<sample>/<step>/logs/<step>.log` (timestamps, inputs, parameters, timings, key numbers,
  Python warnings) and the same lines to the Nextflow task log.
- `results/pipeline_info/` holds the Nextflow execution report, timeline, trace and DAG.
- The HTML report's Provenance section lists software versions and all parameters.

## Parameters

| parameter | default | |
|---|---|---|
| `--input` / `--bundle` | – | samplesheet, or a single bundle (`--sample`, `--he_image`, `--he_alignment`, `--keypoints`, `--annotation`, `--reference_alignment`) |
| `--outdir` | `results` | |
| `--markers` | built-in | marker YAML |
| `--min_counts` | 10 | cells below are removed |
| `--min_cells_per_gene` | 5 | |
| `--mad_k` | 5 | soft-flag threshold for area/count outliers |
| `--target_sum` | 100 | normalization target (Salas et al. 2025) |
| `--bin_um` / `--seam_window` | 20 / 31 | density-map bin size; seam running-median window (strips of 4 µm) |
| `--leiden_resolution` | 1.0 | |
| `--annotation_units` | `he_px` | `um` if the GeoJSON is already in Xenium microns |
| `--skip_multiqc` / `--skip_he_crops` | false | |
| `--unzip_all` | false | |
| `--max_cpus` / `--max_memory` / `--max_time` | 8 / 48.GB / 24.h | caps per task |

## Validation

Run on the 10x *FFPE Human Renal Cell Carcinoma* Xenium Protein dataset (465,534 cells, 77 M decoded spots,
XOA 4.0) on a 10-core / 64 GB laptop: **~2.5 min end to end** (Leiden clustering in `regions` is the slowest step).
The outputs reproduce the original notebook analysis exactly, e.g. 81.8% of gene spots ≥ Q20, 447,626 cells passing
filters, 6.1% maximum seam deviation, 0 negative-control hotspots, and a keypoint fit RMSE of 1.82 µm
(leave-one-out max 6.45 µm) within 0.93 µm of the 10x supplemental alignment everywhere in the capture area.
The `test` profile (a 2.5 × 2 mm crop, zip input) runs in under a minute.

`-profile local` is the tested route; the `conda` environment and the Dockerfile have not been built or tested yet.

## Running steps without Nextflow

Every step is a subcommand of the `xqc` CLI (`xqc --help`), e.g.

```bash
xqc transcripts --bundle outs/ --outdir qc/transcripts
xqc alignment --bundle outs/ --outdir qc/alignment --keypoints keypoints.csv --he-alignment matrix.csv
xqc report --step-dirs qc/* --outdir qc/report --sample my_sample
```

## Development

```bash
.venv/bin/pip install -e ".[multiqc,test]"
.venv/bin/pytest -q tests/test_units.py
```

## References

- Salas SM et al. Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows. *Nat Methods* (2025).
- 10x Genomics, Xenium Onboard Analysis documentation: output files, metrics and analysis-summary alerts.
- nf-core/spatialaxe; MultiQC Xenium Extra plugin.

## License

MIT
