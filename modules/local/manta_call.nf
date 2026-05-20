
process MANTA_CALL {
    label 'manta'
    label 'C8M16T48'
    tag { fam }
    storeDir params.cachedir ? "${params.cachedir}/MANTA_CALL" : null

    input:
    tuple val(fam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)
    tuple path(call_regions), path(call_regions_tbi)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${fam}.MANTA.vcf.gz"
    """
    configManta.py \\
        --referenceFasta $ref_fa \\
        --runDir `pwd -P` \\
        --callRegions $call_regions \\
        --bam ${bam.join(' --bam ')}

    ./runWorkflow.py -j $task.cpus -g ${task.memory.toGiga().intValue()}

    MANTA_ROOT=\$(dirname \$(dirname \$(readlink -f \$(which configManta.py))))
    export PATH=\$PATH:\$MANTA_ROOT/libexec
    convertInversion.py \$MANTA_ROOT/libexec/samtools $ref_fa ./results/variants/diploidSV.vcf.gz \\
        | bgzip --threads $task.cpus > $out_vcf
    tabix $out_vcf
    """
}
