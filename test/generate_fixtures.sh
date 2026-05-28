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
# UCSC hg19 chr20, sequence renamed to GRCh37 naming (no "chr" prefix).
REF="hg19.chr20.fa.gz"
if [[ ! -f $REF ]]; then
    echo "[1/3] Downloading chr20 reference..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz
    bgzip -cd chr20.fa.gz | sed 's/>chr20\b/>20/' | bgzip -c > "$REF"
    rm chr20.fa.gz
else
    echo "[1/3] Reference already exists, skipping."
fi
[[ -f ${REF}.fai ]] || samtools faidx "$REF"

# ── BAMs ───────────────────────────────────────────────────────────────────────
# Stream full chr20 directly from 1000G via HTTPS.
# Strip the full b37 header to SN:20 only, drop reads whose mate maps elsewhere.
#
# Samples: NA19238 (YRI) + NA12878 (CEU singleton).
# NA19238 exercises the family joint-calling path (via -profile test,fam).
# NA12878 exercises the singleton path (via -profile test).
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
STEP=1
for ENTRY in "NA19238:YRI:20130415" "NA12878:CEU:20121211"; do
    SAMPLE="${ENTRY%%:*}"; REST="${ENTRY#*:}"; POP="${REST%%:*}"; DATE="${REST##*:}"
    STEP=$((STEP + 1))
    OUT="${SAMPLE}.hg19.bam"
    if [[ ! -f $OUT ]]; then
        echo "[${STEP}/3] Streaming ${SAMPLE} full chr20..."
        URL="${BASE}/${SAMPLE}/alignment/${SAMPLE}.chrom20.ILLUMINA.bwa.${POP}.low_coverage.${DATE}.bam"
        samtools view -h "$URL" | \
            awk '/^@SQ/{if(/\tSN:20\t/)print; next}
                 /^@/{print; next}
                 $7!="=" && $7!="20" && $7!="*"{next}
                 1' | \
            samtools view -b -o "$OUT"
    else
        echo "[${STEP}/3] ${SAMPLE} BAM already exists, skipping."
    fi
    [[ -f "${OUT}.bai" ]] || samtools index "$OUT"
done

# ── Manifest files ─────────────────────────────────────────────────────────────
# singleton: NA12878 only (-profile test)
cat > test.ped <<EOF
CEU	NA12878	0	0	2	1
EOF
cat > test.bams <<EOF
NA12878	${FIXTURES}/NA12878.hg19.bam
EOF

# family: both samples in one family (-profile test,fam)
cat > test_fam.ped <<EOF
FAM	NA12878	0	0	2	1
FAM	NA19238	0	0	2	1
EOF
cat > test_fam.bams <<EOF
NA12878	${FIXTURES}/NA12878.hg19.bam
NA19238	${FIXTURES}/NA19238.hg19.bam
EOF

echo ""
echo "Done. Fixtures ready in ${FIXTURES}"
echo "Run with: nextflow run . -profile test,conda"
echo "     or:  nextflow run . -profile test,fam,conda"
