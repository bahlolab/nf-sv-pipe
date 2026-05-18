
process MANTA_SPLIT_SAMPLE {
    label 'bcftools'
    label 'C2M2T2'
    tag { sample_id }
    publishDir "${params.progdir}/MANTA/split_sample", mode: 'symlink'

    input:
    tuple val(fam), path(vcf), path(idx), val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}.MANTA.bcf"), path("${sample_id}.MANTA.bcf.csi")

    script:
    """
    bcftools view --threads $task.cpus -s ${sample_id} --min-ac 1 -Ob -o ${sample_id}.MANTA.bcf ${vcf}
    bcftools index --threads ${task.cpus} ${sample_id}.MANTA.bcf
    """
}
