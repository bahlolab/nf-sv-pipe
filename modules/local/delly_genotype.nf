
process DELLY_GENOTYPE {
    label 'delly'
    label 'C2M8T8'
    tag { sam }

    input:
    tuple val(sam), path(bam), path(bai), path(sites_bcf), path(sites_csi)
    tuple path(ref_fa), path(ref_idx)
    path(excl)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DELLY.geno.bcf"
    """
    delly call -h ${task.cpus} -g $ref_fa -x $excl -v $sites_bcf -o $out_bcf $bam
    """
}
