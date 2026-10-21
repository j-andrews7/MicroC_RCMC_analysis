process PAIRTOOLS_SORT_LANE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    // Chunks MUST be provided in their original read order: `pairtools sort` is stable, so the order
    // of tied records (and therefore which duplicate `pairtools dedup` keeps) follows the input order.
    tuple val(meta), path(chunks, stageAs: 'chunks/*')

    output:
    tuple val(meta), path("${prefix}.pairs.gz"), emit: pairs
    tuple val(meta), path("*.log")             , emit: log
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*version //'"), emit: versions_pairtools, topic: versions

    script:
    def args       = task.ext.args ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    def chunk_list = chunks instanceof List ? chunks : [chunks]
    def sort_mem   = task.memory ? "${Math.max(1, (task.memory.toGiga() * 0.8) as int)}G" : '2G'
    """
    mkdir -p sort_tmp
    if ! (pairtools merge --concatenate --nproc-in ${task.cpus} ${chunk_list.join(' ')} \\
        | pairtools sort \\
            --tmpdir sort_tmp \\
            --nproc ${task.cpus} \\
            --memory ${sort_mem} \\
            ${args} \\
            -o ${prefix}.pairs.gz) 2> ${prefix}.sort.log; then
        cat ${prefix}.sort.log >&2
        exit 1
    fi
    rm -rf sort_tmp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    touch ${prefix}.sort.log
    """
}
