#!/usr/bin/env bash
# Generates test fixtures for nf-sv-pipe test profile.
# Requires: samtools, bgzip, wget (install via: conda activate sv_pipe)
# BAMs are streamed directly from 1000G via HTTPS (no full download).
# Idempotent — skips files that already exist.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FIXTURES="${SCRIPT_DIR}/fixtures"
mkdir -p "$FIXTURES"
cd "$FIXTURES"

# ── Reference ──────────────────────────────────────────────────────────────────
# UCSC hg19 chr11 + chr20, sequences renamed to GRCh37 naming (no "chr" prefix).
# Two chromosomes give smoove/lumpy enough inter-chr discordant pairs and
# split reads to avoid crashing on empty lumpy output.
REF="hg19.chr11_chr20.fa.gz"
if [[ ! -f $REF ]]; then
    echo "[1/6] Downloading chr11 + chr20 reference..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz
    { bgzip -cd chr11.fa.gz | sed 's/>chr11\b/>11/';
      bgzip -cd chr20.fa.gz | sed 's/>chr20\b/>20/'; } | bgzip -c > "$REF"
    rm chr11.fa.gz chr20.fa.gz
else
    echo "[1/6] Reference already exists, skipping."
fi
[[ -f ${REF}.fai ]] || samtools faidx "$REF"

# ── BAMs ───────────────────────────────────────────────────────────────────────
# Stream 100k reads per chromosome directly from 1000G via HTTPS.
# The chrom-specific source BAMs are already single-chromosome so no region
# filter is needed. Strip the full b37 header to the relevant SN only, drop
# reads whose mate maps elsewhere (Manta crashes on mate_tid=-1 with stale
# coordinates), round-trip through SAM to re-encode reference IDs correctly,
# then merge chr11+chr20 mini-BAMs into a single sorted BAM per sample.
#
# Samples: YRI trio (NA19238/NA19239/NA19240) + NA12878 singleton (CEU).
# NA12878 exercises the SMOOVE singleton path (skip merge/genotype).
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
READS=100000
STEP=1
for ENTRY in "NA19238:YRI:20130415" "NA19239:YRI:20130415" "NA19240:YRI:20130415" "NA12878:CEU:20121211"; do
    SAMPLE="${ENTRY%%:*}"; REST="${ENTRY#*:}"; POP="${REST%%:*}"; DATE="${REST##*:}"
    STEP=$((STEP + 1))
    OUT="${SAMPLE}.hg19.bam"
    if [[ ! -f $OUT ]]; then
        echo "[${STEP}/6] Streaming ${SAMPLE} (${READS} reads × chr11 + chr20)..."
        for CHR in 11 20; do
            case $CHR in
                11) REGION="11:60000-50783853"  ;;
                20) REGION="20:60000-26319569"  ;;
            esac
            URL="${BASE}/${SAMPLE}/alignment/${SAMPLE}.chrom${CHR}.ILLUMINA.bwa.${POP}.low_coverage.${DATE}.bam"
            samtools view -h "$URL" "$REGION" | \
                awk '/^@SQ/{if(/\tSN:'"${CHR}"'\t/)print; next}
                     /^@/{print; next}
                     $7!="=" && $7!="'"${CHR}"'" && $7!="*"{next}
                     ++n<='"${READS}" | \
                samtools view -b -o "${SAMPLE}.chr${CHR}.tmp.bam"
        done
        samtools merge -f "${SAMPLE}.merged.tmp.bam" \
            "${SAMPLE}.chr11.tmp.bam" "${SAMPLE}.chr20.tmp.bam"
        mv "${SAMPLE}.merged.tmp.bam" "$OUT"
        rm "${SAMPLE}.chr11.tmp.bam" "${SAMPLE}.chr20.tmp.bam"
    else
        echo "[${STEP}/6] ${SAMPLE} BAM already exists, skipping."
    fi
    [[ -f "${OUT}.bai" ]] || samtools index "$OUT"
done

# ── Manifest files ─────────────────────────────────────────────────────────────
# PED: YRI trio (NA19238=mother, NA19239=father, NA19240=daughter)
#    + NA12878 singleton (CEU, unrelated)
cat > test.ped <<EOF
YRI	NA19238	0	0	2	1
YRI	NA19239	0	0	1	1
YRI	NA19240	NA19239	NA19238	2	2
CEU	NA12878	0	0	2	1
EOF

# BAMs manifest with absolute paths (pipeline resolves via Nextflow file())
cat > test.bams <<EOF
NA19238	${FIXTURES}/NA19238.hg19.bam
NA19239	${FIXTURES}/NA19239.hg19.bam
NA19240	${FIXTURES}/NA19240.hg19.bam
NA12878	${FIXTURES}/NA12878.hg19.bam
EOF

echo ""
echo "Done. Fixtures ready in ${FIXTURES}"
echo "Run the test profile with:"
echo "    nextflow run . -profile test,conda"
