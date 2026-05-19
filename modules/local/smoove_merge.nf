
process SMOOVE_MERGE {
    label 'smoove'
    label 'C2M8T2'
    tag { fam }

    input:
    tuple val(fam), path(vcfs), path(indices)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(fam), path(out_vcf)

    script:
    pref    = "${fam}.SMOOVE"
    out_vcf = "${pref}.sites.vcf.gz"
    """
    mkdir tmp && export TMPDIR=tmp
    smoove merge $vcfs \\
        --name $pref \\
        --fasta $ref_fa
    """
}
