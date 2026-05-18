#!/usr/bin/awk -f
# Adapted from https://raw.githubusercontent.com/abyzovlab/CNVnator/refs/heads/master/cnvnator2VCF.pl
# Usage: awk [-v prefix=PREFIX] [-v sample_name=SAMPLE] -f cnvnator2VCF.awk file.calls

BEGIN {
    count = 0
}

FNR == 1 {
    if (sample_name != "") {
        pop_id = sample_name
    } else {
        split(FILENAME, _a, ".")
        pop_id = _a[1]
    }

    print "Reading calls ..." > "/dev/stderr"

    "date '+%Y%m%d'" | getline _date
    close("date '+%Y%m%d'")

    print "##fileformat=VCFv4.1"
    print "##fileDate=" _date
    print "##source=CNVnator"
    print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    print "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"
    print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    print "##INFO=<ID=natorRD,Number=1,Type=Float,Description=\"Normalized RD\">"
    print "##INFO=<ID=natorP1,Number=1,Type=Float,Description=\"e-val by t-test\">"
    print "##INFO=<ID=natorP2,Number=1,Type=Float,Description=\"e-val by Gaussian tail\">"
    print "##INFO=<ID=natorP3,Number=1,Type=Float,Description=\"e-val by t-test (middle)\">"
    print "##INFO=<ID=natorP4,Number=1,Type=Float,Description=\"e-val by Gaussian tail (middle)\">"
    print "##INFO=<ID=natorQ0,Number=1,Type=Float,Description=\"Fraction of reads with 0 mapping quality\">"
    print "##INFO=<ID=natorPE,Number=1,Type=Integer,Description=\"Number of paired-ends support the event\">"
    print "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"Sample genotyped to have the variant\">"
    print "##ALT=<ID=DEL,Description=\"Deletion\">"
    print "##ALT=<ID=DUP,Description=\"Duplication\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">"
    print "##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-ends that support the event\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" pop_id
}

{
    type=$1; coor=$2; len=$3; rd=$4; p1=$5; p2=$6; p3=$7; p4=$8; q0=$9; pe=$10

    isDel = (type == "deletion")
    isDup = (type == "duplication")

    if (!isDel && !isDup) {
        print "Skipping unrecognized event type '" type "'." > "/dev/stderr"
        next
    }

    split(coor, _c, ":")
    chrom = _c[1]
    split(_c[2], _se, "-")
    start = _se[1]
    end   = _se[2]

    count++
    id = (prefix != "" ? prefix "_" : "") "CNVnator_" (isDel ? "del" : "dup") "_" count

    INFO = "END=" end
    if (isDel) INFO = INFO ";SVTYPE=DEL;SVLEN=-" int(len)
    else       INFO = INFO ";SVTYPE=DUP;SVLEN="  int(len)
    INFO = INFO ";IMPRECISE"
    if (rd != "") INFO = INFO ";natorRD=" rd
    if (p1 != "") INFO = INFO ";natorP1=" p1
    if (p2 != "") INFO = INFO ";natorP2=" p2
    if (p3 != "") INFO = INFO ";natorP3=" p3
    if (p4 != "") INFO = INFO ";natorP4=" p4
    if (q0 != "") INFO = INFO ";natorQ0=" q0
    if (pe != "") INFO = INFO ";natorPE=" pe

    if (rd != "") {
        fmt = "GT:CN"
        if (pe != "") fmt = fmt ":PE"

        if      (isDel && rd+0 <  0.25)                  sample = "1/1:0"
        else if (isDel && rd+0 >= 0.25)                  sample = "0/1:1"
        else if (isDup && rd+0 <= 1.75)                  sample = "0/1:2"
        else if (isDup && rd+0 >  1.75 && rd+0 <= 2.25) sample = "1/1:2"
        else if (isDup && rd+0 >  2.25)                  sample = "./1:" sprintf("%.0f", rd+0)
        else                                             { fmt = "GT"; sample = "./." }

        if (pe != "") sample = sample ":" pe
    } else {
        fmt = "GT"
        sample = "./."
    }

    printf "%s\t%s\t%s\tN\t%s\t.\tPASS\t%s\t%s\t%s\n",
        chrom, start, id, (isDel ? "<DEL>" : "<DUP>"), INFO, fmt, sample
}
