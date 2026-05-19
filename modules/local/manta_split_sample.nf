
process MANTA_SPLIT_SAMPLE {
    label 'bcftools'
    label 'C2M2T2'
    tag { sample_id }

    input:
    tuple val(fam), path(vcf), path(idx), val(sample_id)

    output:
    tuple val(sample_id), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sample_id}.MANTA.bcf"
    """
    bcftools view --threads $task.cpus -s ${sample_id} --min-ac 1 -Ob -o ${out_bcf} ${vcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
