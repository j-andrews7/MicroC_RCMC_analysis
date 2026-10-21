process CAT_LOOPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    // meta must contain `out_name`; files must already be in the desired (chromosome) order
    tuple val(meta), path(files, stageAs: 'parts/*')

    output:
    tuple val(meta), path("${meta.out_name}"), emit: loops

    script:
    def file_list = files instanceof List ? files : [files]
    """
    awk 'FNR == 1 && NR != 1 { next } { print }' ${file_list.join(' ')} > ${meta.out_name}
    """

    stub:
    def file_list = files instanceof List ? files : [files]
    """
    awk 'FNR == 1 && NR != 1 { next } { print }' ${file_list.join(' ')} > ${meta.out_name}
    """
}
