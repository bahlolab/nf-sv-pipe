
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
    out_vcf = "${sam}.${caller}.corrected.vcf.gz"
    """
    zcat $vcf > input.vcf
    octopusv correct input.vcf corrected.svcf
    octopusv svcf2vcf -i corrected.svcf -o corrected.vcf
    awk 'BEGIN { OFS = "\\t" }
        /^#/ { print; next }
        {
            svtype = \$8; sub(/.*SVTYPE=/, "", svtype); sub(/;.*/, "", svtype)
            \$3 = "${caller}_" svtype "_" ++n[svtype]
            sub(/SVMETHOD=OctopuSV/, "SVMETHOD=${caller}")
            print
        }' corrected.vcf | gzip -c > $out_vcf
    rm input.vcf corrected.svcf corrected.vcf
    """
}
