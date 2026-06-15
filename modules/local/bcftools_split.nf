
process BCFTOOLS_SPLIT {
    label 'bcftools'
    label 'C2M2T2'
    tag { id }
    publishDir { "${params.outdir}/${suffix}" }, mode: 'copy', saveAs: { fn -> publish ? fn.replaceFirst('^split/', '') : null }

    input:
    tuple val(id), path(vcf), path(idx)
    val  suffix
    val  publish

    output:
    tuple val(id), path("split/*.${suffix}.bcf"), path("split/*.${suffix}.bcf.csi")

    script:
    """
    mkdir -p split
    bcftools +split ${vcf} \\
        -Ob -W \\
        -o split
    for bcf in split/*.bcf; do
        sam=\$(basename "\$bcf" .bcf)
        mv "\$bcf"          "split/\${sam}.${suffix}.bcf"
        mv "\${bcf}.csi"    "split/\${sam}.${suffix}.bcf.csi"
    done
    """
}
