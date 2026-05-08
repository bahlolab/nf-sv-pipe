
process MANTA_CALL {
    label 'manta'
    label 'C8M16T48'
    tag { fam }
    publishDir "${params.progdir}/MANTA/call", mode: 'symlink'

    input:
    tuple val(fam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${fam}.MANTA.vcf.gz"
    """
    configManta.py \\
        --referenceFasta $ref_fa \\
        --runDir `pwd -P` \\
        --bam ${bam.join(' --bam ')}

    ./runWorkflow.py -j $task.cpus -g ${task.memory.toGiga().intValue()}

    mv ./results/variants/diploidSV.vcf.gz $out_vcf
    mv ./results/variants/diploidSV.vcf.gz.tbi ${out_vcf}.tbi
    """
}
