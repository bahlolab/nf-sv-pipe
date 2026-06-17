#!/usr/bin/awk -f
# Add or fix SVLEN in SV VCF/BCF records passed through stdin.
# - Adds SVLEN when absent, computed from END - POS (negative for DEL).
# - Corrects sign when wrong: DEL must be negative, DUP/INV positive.
# - Replaces any existing ##INFO SVLEN header (description may be INS-only)
#   with a general one, or inserts one before #CHROM if absent.
# Only acts on SVTYPE in {DEL, DUP, INV}; all other records pass through unchanged.

BEGIN { OFS = "\t"; svlen_header_seen = 0 }

/^##INFO=<ID=SVLEN,/ {
    svlen_header_seen = 1
    print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    next
}

/^#CHROM/ {
    if (!svlen_header_seen)
        print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    print; next
}

/^#/ { print; next }

{
    svtype = ""; svlen = ""; end = 0; svlen_idx = 0
    n = split($8, info, ";")
    for (i = 1; i <= n; i++) {
        if (info[i] ~ /^SVTYPE=/) svtype = substr(info[i], 8)
        if (info[i] ~ /^SVLEN=/)  { svlen = substr(info[i], 7) + 0; svlen_idx = i }
        if (info[i] ~ /^END=/)    end = substr(info[i], 5) + 0
    }

    if (svtype !~ /^(DEL|DUP|INV)$/ || end == 0) { print; next }

    sv_size = end - $2

    if (svlen_idx == 0) {
        $8 = $8 ";SVLEN=" ((svtype == "DEL") ? -sv_size : sv_size)
    } else if (svtype == "DEL" && svlen > 0) {
        info[svlen_idx] = "SVLEN=" (-svlen)
        $8 = info[1]; for (i = 2; i <= n; i++) $8 = $8 ";" info[i]
    } else if (svtype != "DEL" && svlen < 0) {
        info[svlen_idx] = "SVLEN=" (-svlen)
        $8 = info[1]; for (i = 2; i <= n; i++) $8 = $8 ";" info[i]
    }

    print
}
