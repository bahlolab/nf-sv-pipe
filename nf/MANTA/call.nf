
process call {
    cpus 8
    memory '8 GB'
    time '24 h'
    tag { fam }
    publishDir "${params.progdir}/MANTA/call", mode: 'symlink'
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'

    input:
        tuple val(fam), path(bam), path(bai), path(ref_fa), path(ref_fai)

    output:
        tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${fam}.manta.vcf.gz"
        """
        configManta.py \\
            --referenceFasta $ref_fa \\
            --runDir `pwd -P` \\
            --bam ${bam.join(' --bam ')}
        
        ./runWorkflow.py -j $task.cpus -g 8
    
        mv ./results/variants/diploidSV.vcf.gz $out_vcf
        mv ./results/variants/diploidSV.vcf.gz.tbi ${out_vcf}.tbi
        """
}

