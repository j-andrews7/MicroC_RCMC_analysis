//
// Prepare the bwa-mem2 index and (RCMC only) the balancing blacklist
//

include { BWAMEM2_INDEX  } from '../../../modules/nf-core/bwamem2/index/main'
include { MAKE_BLACKLIST } from '../../../modules/local/make_blacklist/main'

workflow PREPARE_GENOME {

    take:
    fasta         // string: path to genome FASTA (may be null if bwamem2_index is given)
    bwamem2_index // string: path to directory with a pre-built bwa-mem2 index (may be null)
    chromsizes    // string: path to chromosome sizes file
    regions       // string: path to RCMC capture regions BED (null for Micro-C)

    main:
    def bwamem2_exts = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
    def ch_index
    if (bwamem2_index) {
        ch_index = channel.value([[id: 'bwamem2_index'], file(bwamem2_index, checkIfExists: true)])
    } else {
        def fasta_file = file(fasta, checkIfExists: true)
        // Reuse an index sitting next to the FASTA (the layout created by the original Snakemake pipeline)
        def existing = bwamem2_exts.collect { ext -> file("${fasta}${ext}") }
        if (existing.every { index_file -> index_file.exists() }) {
            log.info "Using existing bwa-mem2 index next to the FASTA: ${fasta}.{${bwamem2_exts.collect { ext -> ext.substring(1) }.join(',')}}"
            ch_index = channel.value([[id: fasta_file.name], existing])
        } else {
            BWAMEM2_INDEX(channel.value([[id: fasta_file.name], fasta_file]))
            ch_index = BWAMEM2_INDEX.out.index
        }
    }

    def ch_blacklist = channel.value([])
    if (regions) {
        MAKE_BLACKLIST(
            file(chromsizes, checkIfExists: true),
            file(regions, checkIfExists: true)
        )
        ch_blacklist = MAKE_BLACKLIST.out.blacklist
    }

    emit:
    index     = ch_index     // value channel: [ meta, index ]
    blacklist = ch_blacklist // value channel: blacklist.bed or [] (Micro-C)
}
