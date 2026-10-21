process COOLER_ZOOMIFY_BALANCE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(cool)
    path blacklist      // RCMC only; pass [] for genome-wide balancing
    val resolutions

    output:
    tuple val(meta), path("${prefix}.mcool"), emit: mcool
    tuple val(meta), path("*.log")          , emit: log
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

    script:
    def args         = task.ext.args ?: ''
    prefix           = task.ext.prefix ?: "${meta.id}"
    // Original hard-coded `--nproc 60` for RCMC balancing; the result does not depend on nproc.
    def balance_args = blacklist ? "--balance-args '--nproc ${task.cpus} --blacklist ${blacklist}' " : ''
    // Same command as the original Snakemake `rep_mcools` / `condition_mcools` rules.
    """
    if ! (cooler zoomify \\
            --nproc ${task.cpus} \\
            --balance ${balance_args}\\
            --out ${prefix}.mcool \\
            --resolutions ${resolutions} \\
            ${args} \\
            ${cool}) 2> ${prefix}.log; then
        cat ${prefix}.log >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mcool
    touch ${prefix}.log
    """
}
