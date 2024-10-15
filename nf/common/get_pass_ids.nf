
process get_pass_ids {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "${params.progdir}/${params.caller}/set_id", mode: 'symlink'

    input:
    tuple val(fam), path(vcf), path(index)

    output:
    path(out)

    script:
    out = fam + '.pass_id.gz'
    """
    bcftools view -f "PASS,." $vcf -Ou |
         bcftools query -f '%ID\\n' | 
         gzip > $out
    """
}

