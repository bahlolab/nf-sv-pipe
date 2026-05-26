
process MAKE_CALL_REGIONS {
    label 'bcftools'
    label 'C2M2T2'

    input:
    tuple path(ref_fa), path(ref_idx)
    val(chrs)
    path(excl_tsv)

    output:
    tuple path("call_regions.bed.gz"), path("call_regions.bed.gz.tbi")

    script:
    """
    make_call_regions.awk -v chrs="${chrs.join(',')}" ${ref_fa}.fai $excl_tsv \\
        | sort -k1,1 -k2,2n \\
        | bgzip > call_regions.bed.gz
    tabix -p bed call_regions.bed.gz
    """
}
