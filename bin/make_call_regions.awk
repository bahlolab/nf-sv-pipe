#!/usr/bin/awk -f
# Build a BED of callable regions from a .fai and a delly-style exclude TSV.
# Usage: make_call_regions.awk -v chrs=chr1,chr2 ref.fai excl.tsv > out.bed
#   chrs: comma-separated chromosome names; empty string = include all
#   Output is unsorted BED3; pipe through: sort -k1,1 -k2,2n | bgzip; tabix -p bed

BEGIN {
    nofilter = (chrs == "")
    if (!nofilter) {
        n = split(chrs, a, ",")
        for (i = 1; i <= n; i++) want[a[i]] = 1
    }
}

FNR == NR {
    if (nofilter || ($1 in want)) len[$1] = $2
    next
}

NF == 4 && ($1 in len) {
    chr = $1; start = $2 + 0; end = $3 + 0
    if (start > prev_end[chr]) print chr "\t" prev_end[chr] "\t" start
    if (end > prev_end[chr]) prev_end[chr] = end
}

END {
    for (chr in len)
        if (prev_end[chr] < len[chr]) print chr "\t" prev_end[chr] "\t" len[chr]
}
