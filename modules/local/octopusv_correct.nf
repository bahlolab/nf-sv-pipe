
process OCTOPUSV_CORRECT {
    label 'octopusv'
    label 'C2M4T4'
    tag { "${caller}/${sam}" }
    publishDir "${params.progdir}/octopusv_correct", mode: 'symlink'

    input:
        tuple val(caller), val(sam), path(vcf)

    output:
        tuple val(caller), val(sam), path(out_vcf)

    script:
        out_vcf = "${caller}_${sam}.corrected.vcf.gz"
        """
        zcat $vcf > input.vcf
        octopusv correct input.vcf corrected.svcf
        octopusv svcf2vcf -i corrected.svcf -o corrected.vcf
        gzip -c corrected.vcf > $out_vcf
        rm input.vcf corrected.svcf corrected.vcf
        """
}
