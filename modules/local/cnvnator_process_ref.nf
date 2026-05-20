
process CNVNATOR_PROCESS_REF {
    label 'samtools'
    label 'C2M4T4'

    input:
    tuple path(ref), path(fai)
    val(chrs)

    output:
    tuple path('ref'), path('ref.fai')

    script:
    if (chrs) {
        """
        mkdir ref
        ${chrs.collect { "samtools faidx $ref $it > ref/${it}.fa" }.join('\n')}
        awk 'BEGIN{split("${chrs.join(' ')}",a); for(k in a) s[a[k]]=1} \$1 in s' $fai > ref.fai
        """
    } else {
        """
        mkdir ref
        while IFS=\$'\\t' read -r contig _; do
            samtools faidx $ref "\$contig" > "ref/\$contig.fa"
        done < $fai
        cp $fai ref.fai
        """
    }
}
