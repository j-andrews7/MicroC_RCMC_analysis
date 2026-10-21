process DIFF_MUSTACHE {
    tag "${meta.id}:${meta.resolution}:${meta.chrom}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b7f6d46e94715a3714c926a741977658ddf5182357df3db824e629aebf34c9b7/data'
        : 'community.wave.seqera.io/library/python_numpy_scipy_pandas_pruned:f3778d159bb6fa30'}"

    input:
    // meta must contain `resolution` (bp) and `chrom`; mcool1/mcool2 = condition1/condition2
    tuple val(meta), path(mcool1, stageAs: 'map1/*'), path(mcool2, stageAs: 'map2/*')

    output:
    tuple val(meta), path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2"), emit: loops
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('mustache'), eval("python -c 'from mustache._version import __version__; print(__version__)'"), emit: versions_mustache, topic: versions

    script:
    def args          = task.ext.args ?: ''
    def allow_failure = task.ext.allow_failure == null ? true : task.ext.allow_failure
    prefix            = task.ext.prefix ?: "${meta.id}_${meta.resolution}bp.${meta.chrom}"
    // diff_mustache.py is not exposed as a console script by the mustache-hic package and does
    // `from mustache import ...`, so it must be run as a script from inside the installed package dir.
    """
    DIFF_MUSTACHE=\$(python -c 'import mustache, os; print(os.path.join(os.path.dirname(mustache.__file__), "diff_mustache.py"))')

    set +e
    python "\$DIFF_MUSTACHE" \\
        -f1 ${mcool1} \\
        -f2 ${mcool2} \\
        -r ${meta.resolution} \\
        -ch ${meta.chrom} \\
        -p ${task.cpus} \\
        -o ${prefix} \\
        ${args} > ${prefix}.log 2>&1
    rc=\$?
    set -e

    if [ \$rc -ne 0 ]; then
        cat ${prefix}.log >&2
        if [ \$rc -ge 128 ] || [ "${allow_failure}" != "true" ]; then
            exit \$rc
        fi
        echo "WARNING: diff_mustache exited with status \$rc for ${meta.chrom} at ${meta.resolution} bp; writing empty loop files." | tee -a ${prefix}.log >&2
        rm -f ${prefix}.loop1 ${prefix}.loop2 ${prefix}.diffloop1 ${prefix}.diffloop2
    fi
    for ext in loop1 loop2 diffloop1 diffloop2; do
        if [ ! -s ${prefix}.\$ext ]; then
            printf 'BIN1_CHR\\tBIN1_START\\tBIN1_END\\tBIN2_CHROMOSOME\\tBIN2_START\\tBIN2_END\\tFDR\\tDETECTION_SCALE\\n' > ${prefix}.\$ext
        fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.resolution}bp.${meta.chrom}"
    """
    for ext in loop1 loop2 diffloop1 diffloop2; do
        printf 'BIN1_CHR\\tBIN1_START\\tBIN1_END\\tBIN2_CHROMOSOME\\tBIN2_START\\tBIN2_END\\tFDR\\tDETECTION_SCALE\\n' > ${prefix}.\$ext
    done
    touch ${prefix}.log
    """
}
