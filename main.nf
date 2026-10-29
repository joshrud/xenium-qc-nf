#!/usr/bin/env nextflow
/*
 * xenium-qc-nf: QC and H&E alignment checks for 10x Xenium Onboard Analysis (XOA) output bundles.
 * Companion to nf-core/spatialaxe (same samplesheet columns: sample,bundle,image).
 *
 *   nextflow run main.nf --bundle /path/to/outs.zip --outdir results          (single sample)
 *   nextflow run main.nf --input samplesheet.csv --outdir results             (one or more samples)
 */
nextflow.enable.dsl = 2

include { XENIUM_QC } from './workflows/xenium_qc'

def optFile(v) {
    // Optional path inputs are passed as [] when absent (nf-core convention)
    return (v && v.toString().trim()) ? file(v.toString().trim(), checkIfExists: true) : []
}

def toRow(Map r) {
    def meta = [id: r.sample ?: 'sample']
    def bundle = file(r.bundle.toString().trim(), checkIfExists: true)
    meta.is_zip = !bundle.isDirectory()
    return [meta, bundle, optFile(r.image), optFile(r.alignment), optFile(r.keypoints),
            optFile(r.annotation), optFile(r.reference_alignment)]
}

workflow {
    log.info """
    ==============================================
     xenium-qc-nf  v${workflow.manifest.version}
    ==============================================
     input        : ${params.input ?: params.bundle}
     outdir       : ${params.outdir}
     markers      : ${params.markers ?: 'built-in default'}
     multiqc      : ${!params.skip_multiqc}
     profile      : ${workflow.profile}
    ----------------------------------------------
    """.stripIndent()

    if (params.input) {
        ch_samples = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { toRow(it) }
    } else if (params.bundle) {
        ch_samples = Channel.of(toRow([sample: params.sample, bundle: params.bundle, image: params.he_image,
                                       alignment: params.he_alignment, keypoints: params.keypoints,
                                       annotation: params.annotation, reference_alignment: params.reference_alignment]))
    } else {
        error "Provide --input samplesheet.csv or --bundle <XOA .zip or folder>"
    }

    XENIUM_QC(ch_samples)
}

workflow.onComplete {
    log.info(workflow.success
        ? "Done. Reports: ${params.outdir}/<sample>/report/<sample>_xenium_qc_report.html | pipeline_info: ${params.outdir}/pipeline_info"
        : "Pipeline failed: ${workflow.errorMessage}. See .nextflow.log and ${params.outdir}/<sample>/<step>/logs/")
}
