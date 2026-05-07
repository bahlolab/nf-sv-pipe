
process SVTYPE_SPLIT {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/svtype_split", mode: 'symlink'

    input:
        tuple val(sam), path(vcf), path(tbi)

    output:
        tuple val(sam),
              path(jasmine_vcf), path("${jasmine_vcf}.tbi"),
              path(truvari_vcf), path("${truvari_vcf}.tbi")

    script:
        jasmine_vcf = "${sam}.jasmine.vcf.gz"
        truvari_vcf = "${sam}.truvari.vcf.gz"
        """
        bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \\
            $vcf -Oz -o $jasmine_vcf
        tabix $jasmine_vcf

        bcftools view -i 'INFO/SVTYPE="INS" || INFO/SVTYPE="BND" || INFO/SVTYPE="TRA"' \\
            $vcf -Oz -o $truvari_vcf
        tabix $truvari_vcf
        """
}
