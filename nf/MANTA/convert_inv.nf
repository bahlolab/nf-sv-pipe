
process convert_inv {
    cpus 1
    memory '1 GB'
    time '1 h'
    tag { fam }
    publishDir "progress/MANTA/convert_inv", mode: 'symlink'
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'

    input:
        tuple val(fam), path(vcf), path(tbi), path(ref_fa), path(ref_fai)

    output:
        tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${fam}.manta.ci.vcf.gz"
        """
        DIR=\$(dirname \$(dirname \$(readlink \$(which configManta.py) -f)))/libexec
        PATH=\$DIR:\$PATH
        convertInversion.py \$DIR/samtools $ref_fa $vcf |
            bgzip > $out_vcf
        tabix $out_vcf
        """
}

