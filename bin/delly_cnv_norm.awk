#!/usr/bin/awk -f
# Normalise delly cnv output: recode SVTYPE=CNV records to DEL/DUP based on
# the sample CN field, fix the symbolic allele, and update GT.
# CN=2 (diploid) records are dropped; all other records pass through unchanged.
BEGIN { OFS = "\t" }

/^#CHROM/ {
    if (NF != 10) {
        print "delly_cnv_norm.awk: expected exactly 1 sample, got " (NF - 9) > "/dev/stderr"
        exit 1
    }
    print; next
}

/^#/ { print; next }

{
    n = split($9, fmt, ":")
    cn_idx = 0
    for (i = 1; i <= n; i++) if (fmt[i] == "CN") { cn_idx = i; break }

    if ($8 !~ /SVTYPE=CNV/ || cn_idx == 0) { print; next }

    split($10, samp, ":")
    cn = int(samp[cn_idx] + 0.5)

    if (cn == 2) next

    if (cn <= 1) {
        sub(/SVTYPE=CNV/, "SVTYPE=DEL", $8)
        sub(/<CNV>/, "<DEL>", $5)
        gt = (cn == 0) ? "0/0" : "0/1"
    } else {
        sub(/SVTYPE=CNV/, "SVTYPE=DUP", $8)
        sub(/<CNV>/, "<DUP>", $5)
        gt = (cn >= 4) ? "1/1" : "0/1"
    }

    samp[1] = gt
    $10 = samp[1]; for (i = 2; i <= n; i++) $10 = $10 ":" samp[i]

    print
}
