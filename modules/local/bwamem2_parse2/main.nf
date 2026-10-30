process BWAMEM2_PARSE2 {
    tag "${meta.n_chunks && meta.n_chunks > 1 ? "${meta.id}.part_${meta.chunk}" : meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path chromsizes
    val assembly
    val mapq_min

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    tuple val(meta), path("*.log")     , emit: log
    tuple val("${task.process}"), val('bwamem2'), eval('bwa-mem2 version 2>/dev/null | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem2, topic: versions
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*version //'"), emit: versions_pairtools, topic: versions

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same command as the original Snakemake `make_pairs` rule, minus the trailing `pairtools sort`
    // (sorting happens per lane in PAIRTOOLS_SORT_LANE after all chunks are re-joined in read order).
    // The index prefix is located with bash globs rather than `find`: the Singularity build of the
    // container is a minimal image without findutils. The index is staged either as loose files
    // (index next to the FASTA) or as a directory (--bwamem2_index, or BWAMEM2_INDEX output).
    """
    shopt -s nullglob
    INDEX_FILES=( ./*.0123 ./*/*.0123 )
    shopt -u nullglob
    if [ \${#INDEX_FILES[@]} -ne 1 ]; then
        echo "ERROR: expected exactly one bwa-mem2 index (*.0123) in the staged index, found: \${INDEX_FILES[*]:-none}" >&2
        exit 1
    fi
    INDEX=\${INDEX_FILES[0]%.0123}

    if ! (bwa-mem2 mem -t ${task.cpus} -SP ${args} \$INDEX ${reads[0]} ${reads[1]} \\
        | pairtools parse2 \\
            --add-columns mapq \\
            -c ${chromsizes} \\
            --assembly ${assembly} \\
            --min-mapq ${mapq_min} \\
            --drop-sam \\
            --drop-readid \\
            --expand \\
            --nproc-in ${task.cpus} \\
            ${args2} \\
            -o ${prefix}.pairs.gz) 2> ${prefix}.log; then
        cat ${prefix}.log >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    touch ${prefix}.log
    """
}
