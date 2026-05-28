#!/usr/bin/awk -f
# Adapted from https://raw.githubusercontent.com/abyzovlab/CNVnator/refs/heads/master/cnvnator2VCF.pl
# Usage: awk [-v prefix=PREFIX] [-v sample_name=SAMPLE] [-v exclude_tsv=FILE] [-v min_overlap=FRAC] -f cnvnator2VCF.awk file.calls
#   exclude_tsv: delly-style exclude TSV; 1-field lines = whole-chr exclusion, 4-field = interval exclusion
#   min_overlap: fraction of call length that must overlap a single exclude interval to drop the call (default 0.5)

BEGIN {
    count = 0
    if (min_overlap == "") min_overlap = 0.5
    if (exclude_tsv != "") {
        while ((getline line < exclude_tsv) > 0) {
            n = split(line, f, "\t")
            if (n == 1 && f[1] != "") {
                excl_chr[f[1]] = 1
            } else if (n >= 3) {
                chr = f[1]
                excl_n[chr]++
                excl_s[chr, excl_n[chr]] = f[2] + 0
                excl_e[chr, excl_n[chr]] = f[3] + 0
            }
        }
        close(exclude_tsv)
    }
}

function is_excluded(chr, start, end,    call_len, i, ov_s, ov_e) {
    if (chr in excl_chr) return 1
    if (!(chr in excl_n)) return 0
    call_len = end - start
    if (call_len <= 0) return 0
    for (i = 1; i <= excl_n[chr]; i++) {
        ov_s = (start > excl_s[chr, i]) ? start : excl_s[chr, i]
        ov_e = (end   < excl_e[chr, i]) ? end   : excl_e[chr, i]
        if (ov_e > ov_s && (ov_e - ov_s) / call_len > min_overlap) return 1
    }
    return 0
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

    svtype = isDel ? "DEL" : "DUP"
    n_total[chrom, svtype]++
    if (exclude_tsv != "" && is_excluded(chrom, start + 0, end + 0)) {
        n_dropped[chrom, svtype]++
        next
    }

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

END {
    if (exclude_tsv == "") next
    tot_del = 0; tot_dup = 0; excl_del = 0; excl_dup = 0
    for (k in n_total) {
        split(k, kk, SUBSEP)
        chrs_seen[kk[1]] = 1
        if (kk[2] == "DEL") tot_del += n_total[k]
        else                 tot_dup += n_total[k]
    }
    for (k in n_dropped) {
        split(k, kk, SUBSEP)
        if (kk[2] == "DEL") excl_del += n_dropped[k]
        else                 excl_dup += n_dropped[k]
    }
    print "Excluded " excl_del "/" tot_del " DELs and " excl_dup "/" tot_dup " DUPs:" > "/dev/stderr"
    cmd = "sort -V > /dev/stderr"
    for (chr in chrs_seen) {
        del_tot  = n_total[chr, "DEL"] + 0
        dup_tot  = n_total[chr, "DUP"] + 0
        del_excl = n_dropped[chr, "DEL"] + 0
        dup_excl = n_dropped[chr, "DUP"] + 0
        if (del_tot > 0 && dup_tot > 0)
            line = "  " chr ": " del_excl "/" del_tot " DEL, " dup_excl "/" dup_tot " DUP"
        else if (del_tot > 0)
            line = "  " chr ": " del_excl "/" del_tot " DEL"
        else
            line = "  " chr ": " dup_excl "/" dup_tot " DUP"
        print line | cmd
    }
    close(cmd)
}
