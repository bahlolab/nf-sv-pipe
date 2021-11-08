
process manta_call {
    cpus 8
    memory '8 GB'
    time '24 h'
    tag { fam }
    publishDir "progress/manta_call", mode: 'symlink'

    input:
    tuple val(fam), path(bam), path(bai), path(ref_fa), path(ref_fai)

    output:
    tuple val(fam), path("${fam}.*.vcf.gz")

    script:
    """
    configManta.py \\
        --referenceFasta $ref_fa \\
        --runDir `pwd -P` \\
        --bam ${bam.join(' --bam ')}
    
    ./runWorkflow.py -j 8 -g 16

    for VCF in ./results/variants/*.vcf.gz; do mv \$VCF "$fam.`basename \$VCF`"; done
    """
}

