
include { qdnaseq_bins } from './qdnaseq_bins'
include { mosdepth } from './mosdepth'

workflow QDNASEQ {
    take:
        ref
        fam_bam_ch

    main:
        fam_bam_ch |
            combine(qdnaseq_bins()) |
            mosdepth

//    emit:
}
