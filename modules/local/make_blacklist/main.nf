process MAKE_BLACKLIST {
    tag "${regions}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data'
        : 'community.wave.seqera.io/library/python_bwa-mem2_pairtools_cooler_pruned:b6affcbd92754aec'}"

    input:
    path chromsizes
    path regions, stageAs: 'regions/*'

    output:
    path "chroms.bed"   , emit: chroms
    path "blacklist.bed", emit: blacklist
    tuple val("${task.process}"), val('bedops'), eval("bedops --version 2>&1 | sed -n 's/^ *version: *//p' | cut -d' ' -f1"), emit: versions_bedops, topic: versions
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed 's/bedtools v//'"), emit: versions_bedtools, topic: versions

    script:
    // chroms.bed: equivalent of the original pandas-based `make_chrombed` (chrom, 0, length).
    // blacklist.bed: same command as the original Snakemake `make_blacklist` rule.
    """
    awk 'BEGIN { OFS = "\\t" } { sub(/\\r\$/, "") } NF >= 2 { print \$1, 0, \$2 }' ${chromsizes} > chroms.bed

    bedops --partition <(bedtools sort -i chroms.bed) <(bedtools sort -i ${regions}) \\
        | bedops --not-element-of 1 - <(bedtools sort -i ${regions}) > blacklist.bed
    """

    stub:
    """
    touch chroms.bed
    touch blacklist.bed
    """
}
