process COOLER_CLOAD_PAIRS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(pairs), path(index)
    path chromsizes
    val assembly
    val binsize

    output:
    tuple val(meta), path("${prefix}.cool"), emit: cool
    tuple val(meta), path("*.log")         , emit: log
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // Same command as the original Snakemake `rep_cools` / `condition_cools` rules.
    """
    if ! (bgzip -cd -@ ${task.cpus} ${pairs} \\
        | cooler cload pairs \\
            -c1 2 -p1 3 -c2 4 -p2 5 \\
            --assembly ${assembly} \\
            ${args} \\
            ${chromsizes}:${binsize} \\
            - \\
            ${prefix}.cool) 2> ${prefix}.log; then
        cat ${prefix}.log >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cool
    touch ${prefix}.log
    """
}
