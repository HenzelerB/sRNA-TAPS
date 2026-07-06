#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 12 ]; then
    echo "usage: three_letter_align.sh FASTQ OUTPUT_BAM C2T_INDEX G2A_INDEX PYTHON THREE_LETTER THREADS BOWTIE_LOG FLAGSTAT_LOG PLUS_LOG MINUS_LOG TMPROOT" >&2
    exit 2
fi

FASTQ=$1
OUTPUT_BAM=$2
C2T_INDEX=$3
G2A_INDEX=$4
PYTHON=$5
THREE_LETTER=$6
THREADS=$7
BOWTIE_LOG=$8
FLAGSTAT_LOG=$9
PLUS_LOG=${10}
MINUS_LOG=${11}
TMPROOT=${12}

mkdir -p "$(dirname "$OUTPUT_BAM")" "$(dirname "$BOWTIE_LOG")" "$TMPROOT"
WORK=$(mktemp -d "$TMPROOT/srnataps.3letter.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

"$PYTHON" "$THREE_LETTER" convert-fastq \
    --input "$FASTQ" \
    --output "$WORK/converted.fastq" \
    > "$WORK/convert.log"

WORKERS=$(( THREADS / 2 ))
if [ "$WORKERS" -lt 1 ]; then WORKERS=1; fi

awk -v out="$WORK" -v n="$WORKERS" '
{
    chunk = int((NR - 1) / 4) % n
    file = sprintf("%s/chunk_%03d.fastq", out, chunk)
    print >> file
}
' "$WORK/converted.fastq"

run_branch() {
    local branch=$1
    local index=$2
    local strand=$3
    local branch_log=$4
    local pids=()

    : > "$branch_log"
    for fq in "$WORK"/chunk_*.fastq; do
        (
            bowtie \
                -n 1 -l 10 -e 100 \
                "$strand" \
                -k 1 \
                -p 1 \
                -q --sam \
                -x "$index" "$fq" \
                2>> "$branch_log" \
            | samtools view -bS -F 4 -o "$fq.$branch.bam" -
        ) &
        pids+=("$!")
    done

    local status=0
    for pid in "${pids[@]}"; do
        wait "$pid" || status=1
    done
    if [ "$status" -ne 0 ]; then
        echo "Three-letter $branch alignment failed" >&2
        exit "$status"
    fi

    samtools cat -o "$WORK/$branch.bam" "$WORK"/chunk_*.fastq."$branch".bam
}

run_branch plus "$C2T_INDEX" --norc "$PLUS_LOG"
run_branch minus "$G2A_INDEX" --nofw "$MINUS_LOG"

"$PYTHON" "$THREE_LETTER" restore \
    --original-fastq "$FASTQ" \
    --plus-bam "$WORK/plus.bam" \
    --minus-bam "$WORK/minus.bam" \
    --output-bam "$WORK/restored.bam" \
    > "$WORK/restore.log"

samtools sort -@ "$WORKERS" -o "$OUTPUT_BAM" "$WORK/restored.bam"
samtools index "$OUTPUT_BAM" "$OUTPUT_BAM.bai"
samtools flagstat "$OUTPUT_BAM" > "$FLAGSTAT_LOG"

TOTAL=$(awk 'END { print NR / 4 }' "$WORK/converted.fastq")
MAPPED=$(samtools view -c -F 4 "$OUTPUT_BAM")
FAILED=$(( TOTAL - MAPPED ))

awk -v total="$TOTAL" -v mapped="$MAPPED" -v failed="$FAILED" '
BEGIN {
    mapped_pct = total ? 100 * mapped / total : 0
    failed_pct = total ? 100 * failed / total : 0
    print "# reads processed: " total
    printf "# reads with at least one reported alignment: %d (%.2f%%)\n", mapped, mapped_pct
    printf "# reads that failed to align: %d (%.2f%%)\n", failed, failed_pct
    print "Reported " mapped " alignments"
}
' > "$BOWTIE_LOG"
