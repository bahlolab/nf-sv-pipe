
process manta {
    cpus 8
    memory '16 GB'
    time '48 h'
    tag { fam }
    publishDir "progress/manta", mode: 'symlink'

    input:
    tuple val(fam), path(bam), path(bai), path(ref_fa), path(ref_fai)

    output:
    tuple val(fam), path("*.diploidSV.vcf.gz"), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSV.vcf.gz")

    script:
    """
    /stornext/HPCScratch/home/munro.j/software/manta/bin/configManta.py \\
        --referenceFasta $ref_fa \\
        --runDir `pwd -P` \\
        --bam ${bam.join(' --bam ')}
    
    #./runWorkflow.py -j 8 -g 4

    #for VCF in ./results/variants/*.vcf.gz; do mv \$VCF "$fam.`basename \$VCF`"; done
    """
}