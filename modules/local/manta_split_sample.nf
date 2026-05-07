
process MANTA_SPLIT_SAMPLE {
    label 'bcftools'
    label 'C2M2T2'
    tag { sample_id }
    publishDir "${params.progdir}/MANTA/split_sample", mode: 'symlink'

    input:
    tuple val(fam), path(vcf), path(tbi), val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}.MANTA.vcf.gz"), path("${sample_id}.MANTA.vcf.gz.tbi")

    script:
    """
    bcftools view --threads $task.cpus -s ${sample_id} --min-ac 1 -Oz -o ${sample_id}.MANTA.vcf.gz ${vcf}
    bcftools index -t ${sample_id}.MANTA.vcf.gz
    """
}
