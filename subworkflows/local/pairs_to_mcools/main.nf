//
// Index pairs and build balanced multi-resolution contact maps for replicates and conditions
//

include { PAIRIX_INDEX as PAIRIX_REPS                         } from '../../../modules/local/pairix_index/main'
include { PAIRIX_INDEX as PAIRIX_CONDITIONS                   } from '../../../modules/local/pairix_index/main'
include { COOLER_CLOAD_PAIRS as COOLER_CLOAD_REPS             } from '../../../modules/local/cooler_cload_pairs/main'
include { COOLER_CLOAD_PAIRS as COOLER_CLOAD_CONDITIONS       } from '../../../modules/local/cooler_cload_pairs/main'
include { COOLER_ZOOMIFY_BALANCE as COOLER_ZOOMIFY_REPS       } from '../../../modules/local/cooler_zoomify_balance/main'
include { COOLER_ZOOMIFY_BALANCE as COOLER_ZOOMIFY_CONDITIONS } from '../../../modules/local/cooler_zoomify_balance/main'

workflow PAIRS_TO_MCOOLS {

    take:
    ch_rep_pairs       // channel: [ meta, nodups.pairs.gz ]
    ch_condition_pairs // channel: [ meta, merged_nodups.pairs.gz ]
    ch_chromsizes      // value channel: chromosome sizes file
    ch_blacklist       // value channel: blacklist.bed or []
    assembly           // string
    min_resolution     // integer: base resolution for `cooler cload`
    resolutions        // string: comma-separated resolutions for `cooler zoomify`

    main:
    // original `index_reps` / `index_conditions`
    PAIRIX_REPS(ch_rep_pairs)
    PAIRIX_CONDITIONS(ch_condition_pairs)

    // original `rep_cools` / `condition_cools`
    COOLER_CLOAD_REPS(PAIRIX_REPS.out.index, ch_chromsizes, assembly, min_resolution)
    COOLER_CLOAD_CONDITIONS(PAIRIX_CONDITIONS.out.index, ch_chromsizes, assembly, min_resolution)

    // original `rep_mcools` / `condition_mcools`
    COOLER_ZOOMIFY_REPS(COOLER_CLOAD_REPS.out.cool, ch_blacklist, resolutions)
    COOLER_ZOOMIFY_CONDITIONS(COOLER_CLOAD_CONDITIONS.out.cool, ch_blacklist, resolutions)

    emit:
    rep_mcools       = COOLER_ZOOMIFY_REPS.out.mcool       // channel: [ meta, {condition}_{rep}.mcool ]
    condition_mcools = COOLER_ZOOMIFY_CONDITIONS.out.mcool // channel: [ meta, {condition}_merged.mcool ]
}
