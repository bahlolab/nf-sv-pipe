
process SVTYPE_SPLIT {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/svtype_split", mode: 'symlink'

    input:
    tuple val(sam), path(vcf)

    output:
    tuple val(sam), path(jasmine_vcf), path(truvari_vcf)

    script:
    jasmine_vcf = "${sam}.jasmine.vcf.gz"
    truvari_vcf = "${sam}.truvari.vcf.gz"
    """
    echo "$sam" > sample_name.txt

    bcftools view --threads $task.cpus -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \\
        $vcf -Oz | bcftools reheader -s sample_name.txt -o $jasmine_vcf

    bcftools view --threads $task.cpus -i 'INFO/SVTYPE="INS" || INFO/SVTYPE="BND" || INFO/SVTYPE="TRA"' \\
        $vcf -Oz | bcftools reheader -s sample_name.txt -o $truvari_vcf

    rm sample_name.txt
    """
}
