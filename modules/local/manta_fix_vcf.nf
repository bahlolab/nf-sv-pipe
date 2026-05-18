
process MANTA_FIX_VCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { id }
    publishDir "${params.progdir}/MANTA/fix_vcf", mode: 'symlink'

    input:
    tuple val(id), path(vcf), path(tbi)

    output:
    tuple val(id), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = vcf.name.replaceAll('.vcf.gz', '.fixed.bcf')
    """
    bcftools view -h $vcf | \\
        sed 's:##INFO=<ID=SVLEN,Number=.,:##INFO=<ID=SVLEN,Number=1,:' > header.txt
    bcftools reheader $vcf -h header.txt | bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
