process UNPACK {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(zip)

    output:
    tuple val(meta), path("bundle"), emit: bundle
    tuple val(meta), path("unpack"), emit: dir

    script:
    def all = params.unzip_all ? '--all' : ''
    """
    xqc unpack --input ${zip} --outdir unpack ${all}
    mv unpack/bundle bundle
    """
}
