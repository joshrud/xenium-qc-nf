/*
 * Per-sample QC workflow. Input tuples:
 *   [ meta, bundle(.zip|dir), he_image, alignment, keypoints, annotation, reference_alignment ]  ([] = not provided)
 * The alignment used for H&E crops and annotation regions is `alignment` if given, else `reference_alignment`.
 */
include { UNPACK         } from '../modules/local/unpack/main'
include { RUN_METRICS    } from '../modules/local/run_metrics/main'
include { TRANSCRIPT_QC  } from '../modules/local/transcript_qc/main'
include { CELL_QC        } from '../modules/local/cell_qc/main'
include { PROTEIN_QC     } from '../modules/local/protein_qc/main'
include { MORPHOLOGY_QC  } from '../modules/local/morphology_qc/main'
include { ALIGNMENT_QC   } from '../modules/local/alignment_qc/main'
include { REGIONS        } from '../modules/local/regions/main'
include { MULTIQC_XENIUM } from '../modules/local/multiqc_xenium/main'
include { QC_REPORT      } from '../modules/local/qc_report/main'

workflow XENIUM_QC {
    take:
    ch_samples

    main:
    markers = params.markers ? file(params.markers, checkIfExists: true) : []

    // ---- unpack .zip bundles; folders are used in place
    ch_split = ch_samples.branch { zip: it[0].is_zip; dir: true }
    UNPACK(ch_split.zip.map { [it[0], it[1]] })
    ch_bundle = UNPACK.out.bundle.mix(ch_split.dir.map { [it[0], it[1]] })

    // side inputs keyed by sample id
    ch_extra = ch_samples.map { meta, b, img, aln, kp, ann, ref -> [meta.id, img, aln ?: ref, kp, ann, ref] }
    ch_in = ch_bundle.map { meta, b -> [meta.id, meta, b] }.join(ch_extra)
                     .map { id, meta, b, img, aln, kp, ann, ref -> [meta, b, img, aln, kp, ann, ref] }

    RUN_METRICS(ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta, b] })
    TRANSCRIPT_QC(ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta, b] })
    CELL_QC(ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta, b] }, markers)

    ch_cells = CELL_QC.out.cells_qc.map { meta, p -> [meta.id, p] }
    ch_h5ad  = CELL_QC.out.h5ad.map { meta, p -> [meta.id, p] }
    ch_keyed = ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta.id, meta, b, img, aln, kp, ann, ref] }

    PROTEIN_QC(ch_keyed.join(ch_cells).map { id, meta, b, img, aln, kp, ann, ref, cq -> [meta, b, cq] }, markers)
    MORPHOLOGY_QC(
        ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta, b, params.skip_he_crops ? [] : img, params.skip_he_crops ? [] : aln] },
        markers)
    ALIGNMENT_QC(ch_samples.map { meta, b, img, aln, kp, ann, ref -> [meta.id, aln, kp, ref] }
                     .join(ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta.id, meta, b] })
                     .map { id, aln, kp, ref, meta, b -> [meta, b, kp, aln, ref] })
    REGIONS(ch_keyed.join(ch_h5ad).map { id, meta, b, img, aln, kp, ann, ref, h5 -> [meta, b, h5, ann, aln] }, markers)

    ch_mqc = Channel.empty()
    if (!params.skip_multiqc) {
        MULTIQC_XENIUM(ch_in.map { meta, b, img, aln, kp, ann, ref -> [meta, b] })
        ch_mqc = MULTIQC_XENIUM.out.report.map { meta, p -> [meta.id, p] }
    }

    // ---- collect every step folder per sample, then build the report
    ch_steps = RUN_METRICS.out.dir
        .mix(TRANSCRIPT_QC.out.dir, CELL_QC.out.dir, PROTEIN_QC.out.dir, MORPHOLOGY_QC.out.dir,
             ALIGNMENT_QC.out.dir, REGIONS.out.dir)
        .map { meta, d -> [meta.id, meta, d] }
        .groupTuple(size: 7)
        .map { id, metas, dirs -> [id, metas[0], dirs] }
    ch_report_in = params.skip_multiqc
        ? ch_steps.map { id, meta, dirs -> [meta, dirs, []] }
        : ch_steps.join(ch_mqc).map { id, meta, dirs, mqc -> [meta, dirs, mqc] }
    QC_REPORT(ch_report_in)

    emit:
    report = QC_REPORT.out.html
}
