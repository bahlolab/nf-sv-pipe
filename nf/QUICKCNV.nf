
params.caller = 'QUICKCNV'
params.bin_size = 10000

include { MOSDEPTH } from './QUICKCNV/mosdepth'
include { SNORM    } from './QUICKCNV/snorm'
include { BNORM    } from './QUICKCNV/bnorm'
include { CALL     } from './QUICKCNV/call'

workflow QUICKCNV {
    take:
    ref
    fam_bam_ch

    main:
        
    sam_bam_ch = fam_bam_ch.map { it.drop(1) }

    MOSDEPTH(
        sam_bam_ch
    )

    SNORM(
        MOSDEPTH.out
    )

    BNORM(
        SNORM.out.bins.toSortedList { it.name }.map{it[0]},
        SNORM.out.coverage.toSortedList { it.name },
    )

    CALL(
        BNORM.out
            .flatten()
            .map { [it.name.replaceAll('.bnorm.rds', ''), it] }
    )

    emit:
    null
}
