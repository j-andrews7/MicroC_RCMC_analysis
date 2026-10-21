process MUSTACHE {
    tag "${meta.id}:${meta.resolution}:${meta.chrom}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b7f6d46e94715a3714c926a741977658ddf5182357df3db824e629aebf34c9b7/data'
        : 'community.wave.seqera.io/library/python_numpy_scipy_pandas_pruned:f3778d159bb6fa30'}"

    input:
    // meta must contain `resolution` (bp) and `chrom`
    tuple val(meta), path(mcool)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: loops
    tuple val(meta), path("*.log")        , emit: log
    tuple val("${task.process}"), val('mustache'), eval("python -c 'from mustache._version import __version__; print(__version__)'"), emit: versions_mustache, topic: versions

    script:
    def args          = task.ext.args ?: ''
    def allow_failure = task.ext.allow_failure == null ? true : task.ext.allow_failure
    prefix            = task.ext.prefix ?: "${meta.id}_${meta.resolution}bp.${meta.chrom}"
    // Exit codes >= 128 (signals, e.g. OOM kills) are always propagated so Nextflow can retry with more memory.
    """
    set +e
    mustache \\
        -f ${mcool} \\
        -r ${meta.resolution} \\
        -ch ${meta.chrom} \\
        -p ${task.cpus} \\
        -o ${prefix}.tsv \\
        ${args} > ${prefix}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -ne 0 ]; then
        cat ${prefix}.log >&2
        if [ \$rc -ge 128 ] || [ "${allow_failure}" != "true" ]; then
            exit \$rc
        fi
        echo "WARNING: mustache exited with status \$rc for ${meta.chrom} at ${meta.resolution} bp; writing an empty loop file." | tee -a ${prefix}.log >&2
        rm -f ${prefix}.tsv
    fi
    if [ ! -s ${prefix}.tsv ]; then
        printf 'BIN1_CHR\\tBIN1_START\\tBIN1_END\\tBIN2_CHROMOSOME\\tBIN2_START\\tBIN2_END\\tFDR\\tDETECTION_SCALE\\n' > ${prefix}.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.resolution}bp.${meta.chrom}"
    """
    printf 'BIN1_CHR\\tBIN1_START\\tBIN1_END\\tBIN2_CHROMOSOME\\tBIN2_START\\tBIN2_END\\tFDR\\tDETECTION_SCALE\\n' > ${prefix}.tsv
    touch ${prefix}.log
    """
}
