process PAIRTOOLS_MERGE_REPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(pairs, stageAs: 'reps/*')

    output:
    tuple val(meta), path("${prefix}_nodups.pairs.gz"), emit: pairs
    tuple val(meta), path("*.log")                    , emit: log
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*version //'"), emit: versions_pairtools, topic: versions

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def pair_list = pairs instanceof List ? pairs : [pairs]
    // Same command as the original Snakemake `merge_reps_by_condition` rule.
    """
    mkdir -p sort_tmp
    if ! (pairtools merge \\
            --tmpdir sort_tmp \\
            --nproc ${task.cpus} \\
            --output ${prefix}_nodups.pairs.gz \\
            ${args} \\
            ${pair_list.join(' ')}) 2> ${prefix}.log; then
        cat ${prefix}.log >&2
        exit 1
    fi
    rm -rf sort_tmp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_nodups.pairs.gz
    touch ${prefix}.log
    """
}
