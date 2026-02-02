/*
    Note: INV/BND not currently possible as not handled well by jasmine
*/
process jasmine_merge {
    cpus 2
    memory '32 GB'
    time '2 h'
    publishDir "${params.progdir}/jasmine_merge", mode: 'symlink'

    input:
    tuple path(vcfs), path(indices)

    output:
    tuple path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "$params.id" + '.' + "$params.caller" + '.merged.vcf.gz'
    """
    jasmine_merge.py $vcfs \\
         --jasmine-dir \$(dirname \$(which jasmine)) \\
         --output $out_vcf \\
         --sv-type INS \\
         --sv-type DEL \\
         --sv-type DUP \\
         --args 'max_dist_linear=0.20' \\
         --args 'max_dist=500' \\
         --args 'min_overlap=0.80' \\
         --args 'threads=2'
    """
}

