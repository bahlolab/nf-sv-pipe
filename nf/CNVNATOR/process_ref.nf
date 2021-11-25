
process process_ref {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir 'progress/CNVNATOR/process_ref', mode: 'symlink'

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
