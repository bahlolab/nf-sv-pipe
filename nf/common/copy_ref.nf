
process copy_ref {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "${params.progdir}/copy_ref", mode: 'symlink'
    container null

    input:
        tuple path(ref), path(fai)

    output:
        tuple path(ref), path(fai)

    script:
        """
        cp `readlink $ref` tmp && mv tmp $ref
        cp `readlink $fai` tmp && mv tmp $fai
        """
}
