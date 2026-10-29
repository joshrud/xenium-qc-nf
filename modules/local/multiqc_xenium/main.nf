/*
 * MultiQC core Xenium module + Xenium Extra plugin (the QC used by nf-core/spatialaxe --mode qc).
 * Failures are logged but do not stop the pipeline; the main HTML report is built either way.
 */
process MULTIQC_XENIUM {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle')

    output:
    tuple val(meta), path("multiqc/multiqc_report.html"), emit: report
    tuple val(meta), path("multiqc"), emit: dir

    script:
    """
    mkdir -p multiqc
    python -m multiqc bundle/ --outdir multiqc --filename multiqc_report.html --force --title "${meta.id}" \\
        > multiqc/multiqc.log 2>&1 \\
      || { echo "<html><body><p>MultiQC failed; see multiqc.log</p></body></html>" > multiqc/multiqc_report.html; }
    """
}
