
params.caller = 'QUICKCNV'
params.bin_size = 500
params.n_phases = 4
params.n_shards = 50

include { MOSDEPTH   } from './QUICKCNV/mosdepth'
include { NUC         } from './QUICKCNV/nuc'
include { SNORM      } from './QUICKCNV/snorm'
include { BNORM      } from './QUICKCNV/bnorm'
include { CALL       } from './QUICKCNV/call'
include { BPT_DEPTH  } from './QUICKCNV/bpt_depth'
include { REFINE     } from './QUICKCNV/refine'
include { MERGE      } from './QUICKCNV/merge'
include { FIXVCF     } from './QUICKCNV/fixvcf'


workflow QUICKCNV {
    take:
    ref
    fam_bam_ch

    main:
        
    sam_bam_ch = fam_bam_ch.map { it.drop(1) }
    
    NUC(
        ref
    )

    MOSDEPTH(
        sam_bam_ch
    )
    
    SNORM(
        MOSDEPTH.out,
        NUC.out
    )

    bins = SNORM.out.bins
        .flatten()
        .map { [(it.name =~ /\.shard_([0-9]+)\./)[0][1], it] }
        .groupTuple(by:0)
        .map { [it[0], it[1].sort { it.name }[0] ] }
    
    snorm = SNORM.out.snorm 
        .flatten()
        .map { [(it.name =~ /\.shard_([0-9]+)\./)[0][1], it] }
        .groupTuple(by:0)
        .map { [it[0], it[1].sort { it.name } ] }

    BNORM(
        bins.combine(snorm, by: 0)
    )

    CALL(
        BNORM.out
            .flatten()
            .map { [(it.name =~ /(.+)\.shard_[0-9]+\.bnorm\.rds/)[0][1], it] }
            .groupTuple(by:0)
            .map { [it[0], it[1].sort { it.name } ] }
    )

    // BNORM.out.flatten().map { it.toString() }.collectFile(name: './bnorm.txt', newLine: true)

    BPT_DEPTH(
        sam_bam_ch.combine(CALL.out.reg, by:0)
    )

    REFINE(
        CALL.out.calls.combine(BPT_DEPTH.out, by:0)
    )

    MERGE(
        REFINE.out.collect()
    )

    FIXVCF(
        MERGE.out,
        ref
    )

    emit:
    null
}
