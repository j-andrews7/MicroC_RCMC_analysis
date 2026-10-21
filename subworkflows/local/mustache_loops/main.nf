//
// Loop calling (mustache) and differential loop calling (diff_mustache), parallelised per chromosome
//

include { MUSTACHE                       } from '../../../modules/local/mustache/main'
include { DIFF_MUSTACHE                  } from '../../../modules/local/diff_mustache/main'
include { CAT_LOOPS as CAT_MUSTACHE      } from '../../../modules/local/cat_loops/main'
include { CAT_LOOPS as CAT_DIFF_MUSTACHE } from '../../../modules/local/cat_loops/main'

workflow MUSTACHE_LOOPS {

    take:
    ch_mcools    // channel: [ meta (condition), mcool ]
    ch_contrasts // channel: [ condition1, condition2 ]
    resolutions  // list of integers
    chromosomes  // list of chromosome names (defines output order)
    call_loops   // boolean
    diff_loops   // boolean

    main:
    def res_chroms = []
    resolutions.each { res ->
        chromosomes.eachWithIndex { chrom, idx -> res_chroms << [res, chrom, idx] }
    }
    def n_chroms = chromosomes.size()

    def ch_loops = channel.empty()
    if (call_loops) {
        def ch_mustache_in = ch_mcools
            .combine(channel.fromList(res_chroms))
            .map { meta, mcool, res, chrom, idx ->
                [meta + [resolution: res, chrom: chrom, chrom_idx: idx], mcool]
            }

        MUSTACHE(ch_mustache_in)

        def ch_cat_loops = MUSTACHE.out.loops
            .map { meta, loops ->
                def key = [
                    id        : "${meta.id}_${meta.resolution}bp".toString(),
                    condition : meta.condition,
                    resolution: meta.resolution,
                    out_name  : "${meta.id}_${meta.resolution}bp.mustache.tsv".toString(),
                ]
                [groupKey(key, n_chroms), [meta.chrom_idx, loops]]
            }
            .groupTuple()
            .map { key, items -> [key.getGroupTarget(), sortChromItems(items)] }

        CAT_MUSTACHE(ch_cat_loops)
        ch_loops = CAT_MUSTACHE.out.loops
    }

    def ch_diff_loops = channel.empty()
    if (diff_loops) {
        def ch_by_condition = ch_mcools.map { meta, mcool -> [meta.condition, mcool] }

        def ch_diff_in = ch_contrasts
            .combine(ch_by_condition, by: 0)
            .map { c1, c2, mcool1 -> [c2, c1, mcool1] }
            .combine(ch_by_condition, by: 0)
            .map { c2, c1, mcool1, mcool2 ->
                [[id: "${c1}_vs_${c2}".toString(), condition1: c1, condition2: c2], mcool1, mcool2]
            }
            .combine(channel.fromList(res_chroms))
            .map { meta, mcool1, mcool2, res, chrom, idx ->
                [meta + [resolution: res, chrom: chrom, chrom_idx: idx], mcool1, mcool2]
            }

        DIFF_MUSTACHE(ch_diff_in)

        def ch_cat_diff = DIFF_MUSTACHE.out.loops
            .flatMap { meta, loop1, loop2, diffloop1, diffloop2 ->
                [['loop1', loop1], ['loop2', loop2], ['diffloop1', diffloop1], ['diffloop2', diffloop2]].collect { ext, loop_file ->
                    def key = [
                        id        : "${meta.id}_${meta.resolution}bp.${ext}".toString(),
                        contrast  : meta.id,
                        resolution: meta.resolution,
                        out_name  : "${meta.id}_${meta.resolution}bp.${ext}".toString(),
                    ]
                    [groupKey(key, n_chroms), [meta.chrom_idx, loop_file]]
                }
            }
            .groupTuple()
            .map { key, items -> [key.getGroupTarget(), sortChromItems(items)] }

        CAT_DIFF_MUSTACHE(ch_cat_diff)
        ch_diff_loops = CAT_DIFF_MUSTACHE.out.loops
    }

    emit:
    loops      = ch_loops      // channel: [ meta, {condition}_merged_{res}bp.mustache.tsv ]
    diff_loops = ch_diff_loops // channel: [ meta, {c1}_vs_{c2}_{res}bp.{loop1,loop2,diffloop1,diffloop2} ]
}

def sortChromItems(items) {
    return items.sort(false) { a, b -> a[0] <=> b[0] }.collect { item -> item[1] }
}
