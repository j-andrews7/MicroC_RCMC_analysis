process PAIRIX_INDEX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path(pairs), path("*.px2"), emit: index
    tuple val("${task.process}"), val('pairix'), eval("pairix --help 2>&1 | sed -n 's/^Version: //p'"), emit: versions_pairix, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    pairix ${args} ${pairs}
    """

    stub:
    """
    touch ${pairs}.px2
    """
}
