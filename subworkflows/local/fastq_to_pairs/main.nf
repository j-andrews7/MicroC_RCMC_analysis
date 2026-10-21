//
// FASTQ -> per-replicate deduplicated pairs -> per-condition merged pairs
//

include { SEQKIT_SPLIT2         } from '../../../modules/nf-core/seqkit/split2/main'
include { BWAMEM2_PARSE2        } from '../../../modules/local/bwamem2_parse2/main'
include { PAIRTOOLS_SORT_LANE   } from '../../../modules/local/pairtools_sort_lane/main'
include { PAIRTOOLS_MERGE_DEDUP } from '../../../modules/local/pairtools_merge_dedup/main'
include { PAIRTOOLS_MERGE_REPS  } from '../../../modules/local/pairtools_merge_reps/main'

workflow FASTQ_TO_PAIRS {

    take:
    ch_reads      // channel: [ meta, [ fastq_1, fastq_2 ] ] (one item per lane)
    ch_index      // value channel: [ meta, bwa-mem2 index ]
    ch_chromsizes // value channel: chromosome sizes file
    assembly      // string
    mapq_min      // integer
    split_fastq   // integer: read pairs per chunk (0 = no chunking)

    main:
    //
    // Split each lane into chunks of `split_fastq` read pairs and align them in parallel
    //
    def ch_chunks
    if (split_fastq > 0) {
        SEQKIT_SPLIT2(ch_reads)
        ch_chunks = SEQKIT_SPLIT2.out.reads.flatMap { meta, files -> pairChunks(meta, files) }
    } else {
        ch_chunks = ch_reads.map { meta, reads -> [meta + [chunk: 1, n_chunks: 1], reads] }
    }

    BWAMEM2_PARSE2(ch_chunks, ch_index, ch_chromsizes, assembly, mapq_min)

    //
    // Re-join chunks in read order and sort per lane (equivalent to the original `make_pairs` rule)
    //
    def ch_lane_chunks = BWAMEM2_PARSE2.out.pairs
        .map { meta, pairs ->
            def lane_meta = meta.findAll { key, _value -> !(key in ['chunk', 'n_chunks']) }
            [groupKey(lane_meta, meta.n_chunks), [meta.chunk, pairs]]
        }
        .groupTuple()
        .map { key, chunks -> [key.getGroupTarget(), sortByIndex(chunks)] }

    PAIRTOOLS_SORT_LANE(ch_lane_chunks)

    //
    // Merge lanes per replicate and deduplicate (original `merge_and_dedup_lanes` rule)
    //
    def ch_rep_lanes = PAIRTOOLS_SORT_LANE.out.pairs
        .map { meta, pairs ->
            def rep_meta = [id: "${meta.condition}_${meta.rep}".toString(), condition: meta.condition, rep: meta.rep, n_reps: meta.n_reps]
            [groupKey(rep_meta, meta.n_lanes), [meta.lane_order, pairs]]
        }
        .groupTuple()
        .map { key, lanes -> [key.getGroupTarget(), sortByIndex(lanes)] }

    PAIRTOOLS_MERGE_DEDUP(ch_rep_lanes)

    //
    // Merge replicates per condition (original `merge_reps_by_condition` rule)
    //
    def ch_condition_reps = PAIRTOOLS_MERGE_DEDUP.out.nodups
        .map { meta, pairs ->
            def condition_meta = [id: "${meta.condition}_merged".toString(), condition: meta.condition]
            [groupKey(condition_meta, meta.n_reps), [meta.rep, pairs]]
        }
        .groupTuple()
        .map { key, reps -> [key.getGroupTarget(), sortByIndex(reps)] }

    PAIRTOOLS_MERGE_REPS(ch_condition_reps)

    emit:
    rep_pairs       = PAIRTOOLS_MERGE_DEDUP.out.nodups // channel: [ meta, {condition}_{rep}_nodups.pairs.gz ]
    dedup_stats     = PAIRTOOLS_MERGE_DEDUP.out.stats  // channel: [ meta, {condition}_{rep}.dedup.stats ]
    condition_pairs = PAIRTOOLS_MERGE_REPS.out.pairs   // channel: [ meta, {condition}_merged_nodups.pairs.gz ]
}

//
// Sort [ index, file ] pairs by index and return the files
//
def sortByIndex(items) {
    return items.sort(false) { a, b -> a[0] <=> b[0] }.collect { item -> item[1] }
}

//
// Pair up the R1/R2 chunk files written by `seqkit split2` (<name>.part_NNN.<ext>) in part order
//
def pairChunks(meta, files) {
    def file_list = files instanceof List ? files : [files]
    def r1_base = meta.fastq1_name.replaceAll(/\.f(ast)?q(\.gz)?$/, '')
    def r2_base = meta.fastq2_name.replaceAll(/\.f(ast)?q(\.gz)?$/, '')
    def parts = [:]
    file_list.each { chunk_file ->
        def name = chunk_file.name
        def matcher = name =~ /\.part_(\d+)\.[^\/]+$/
        def part = matcher.find() ? (matcher.group(1) as Integer) : 1
        def mate = -1
        if (name.startsWith("${r1_base}.part_") || name == meta.fastq1_name) {
            mate = 0
        } else if (name.startsWith("${r2_base}.part_") || name == meta.fastq2_name) {
            mate = 1
        } else {
            error("Cannot assign split FASTQ chunk '${name}' to read 1 or read 2 of ${meta.id}")
        }
        if (!parts.containsKey(part)) {
            parts[part] = [null, null]
        }
        parts[part][mate] = chunk_file
    }
    def part_ids = parts.keySet().sort()
    return part_ids.collect { part ->
        def pair = parts[part]
        if (pair[0] == null || pair[1] == null) {
            error("Missing read 1 or read 2 chunk for part ${part} of ${meta.id}")
        }
        [meta + [chunk: part, n_chunks: part_ids.size()], pair]
    }
}
