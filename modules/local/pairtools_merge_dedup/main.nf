process PAIRTOOLS_MERGE_DEDUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(pairs, stageAs: 'lanes/*')

    output:
    tuple val(meta), path("${prefix}_nodups.pairs.gz")  , emit: nodups
    tuple val(meta), path("${prefix}_unmapped.pairs.gz"), emit: unmapped
    tuple val(meta), path("${prefix}_dups.pairs.gz")    , emit: dups
    tuple val(meta), path("${prefix}.dedup.stats")      , emit: stats
    tuple val(meta), path("*.log")                      , emit: log
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*version //'"), emit: versions_pairtools, topic: versions

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def pair_list = pairs instanceof List ? pairs : [pairs]
    // Same command as the original Snakemake `merge_and_dedup_lanes` rule.
    """
    mkdir -p sort_tmp
    if ! (pairtools merge --tmpdir sort_tmp --nproc ${task.cpus} ${pair_list.join(' ')} \\
        | pairtools dedup \\
            --max-mismatch 1 \\
            --mark-dups \\
            --output ${prefix}_nodups.pairs.gz \\
            --output-unmapped ${prefix}_unmapped.pairs.gz \\
            --output-dups ${prefix}_dups.pairs.gz \\
            --output-stats ${prefix}.dedup.stats \\
            ${args} \\
        | cat) 2> ${prefix}.log; then
        cat ${prefix}.log >&2
        exit 1
    fi
    rm -rf sort_tmp
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_nodups.pairs.gz
    echo "" | gzip > ${prefix}_unmapped.pairs.gz
    echo "" | gzip > ${prefix}_dups.pairs.gz
    touch ${prefix}.dedup.stats
    touch ${prefix}.log
    """
}
