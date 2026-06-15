
process MANTA_FIX_VCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.outdir}/MANTA", mode: 'copy'

    input:
    tuple val(sam), path(in_file, stageAs: 'staged/*'), path(in_idx, stageAs: 'staged/*')

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.MANTA.bcf"
    """
    bcftools view -h ${in_file} | \\
        sed 's:##INFO=<ID=SVLEN,Number=.,:##INFO=<ID=SVLEN,Number=1,:' > header.txt
    bcftools reheader ${in_file} -h header.txt \\
        | bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
