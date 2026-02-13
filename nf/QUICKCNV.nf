
params.caller = 'QUICKCNV'
params.bin_size = 500
params.n_phases = 4
params.n_shards = 50

include { MOSDEPTH      } from './QUICKCNV/mosdepth'
include { SNORM         } from './QUICKCNV/snorm'
include { BNORM         } from './QUICKCNV/bnorm'
include { CALL          } from './QUICKCNV/call'
include { VCF_HEADER    } from './QUICKCNV/vcf_header'
include { jasmine_merge as JASMINE } from './common/jasmine_merge'
include { publish_vcf as PUBLISH   } from './common/publish_vcf'

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

    VCF_HEADER(
        CALL.out,
        ref
    )

    VCF_HEADER.out
            .toSortedList()
            .map {it.transpose() }
            
    JASMINE(
        VCF_HEADER.out
            .toSortedList()
            .map {it.transpose() }
    )

    PUBLISH(
        JASMINE.out
    )

    emit:
    null
}
