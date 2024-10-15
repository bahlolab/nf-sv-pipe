
process split_sv_types {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "${params.progdir}/MANTA/split_sv_types", mode: 'symlink'\
    tag { fam }

    input:
        tuple val(fam), path(vcf), path(index)

    output:
        tuple val(fam), path("$pref*.vcf.gz"), path("$pref*.vcf.gz.tbi")

    script:
        pref = vcf.name.replaceAll(/(\.vcf\.gz)|(\.bcf)$/, '')
        """
        bcftools view --no-version $vcf |
             split_sv_types.py --input - \\
                --index \\
                --pref $pref
        """
}
