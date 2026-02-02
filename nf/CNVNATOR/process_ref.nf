
process process_ref {
    label 'C1M1T1'
    publishDir "${params.progdir}/CNVNATOR/process_ref", mode: 'symlink'

    input:
        tuple path(ref), path(fai), val(chrs)

    output:
        tuple path('ref'), path('ref.fai')

    script:
        """
        mkdir ref
        ${chrs.collect { "samtools faidx $ref $it > ref/${it}.fa" }.join ('\n')}
        grep -P '${chrs.collect { it + '\\t' }.join('|')}' $fai > ref.fai
        """
}
