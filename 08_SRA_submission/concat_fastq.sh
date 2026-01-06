#!/usr/bin/env bash
set -euo pipefail

OUTDIR="concat"
LOG="${OUTDIR}/concat.log"

mkdir -p "$OUTDIR"
echo "Concatenation started: $(date)" > "$LOG"
echo "----------------------------------------" >> "$LOG"

for d in A*/; do
    sample_dir="${d%/}"

    l7_file=$(ls "${sample_dir}"/*_L7_1.fq.gz 2>/dev/null | head -n 1)

    if [[ -z "$l7_file" ]]; then
        echo "WARNING: No L7 files found in ${sample_dir}, skipping" >> "$LOG"
        continue
    fi

    base=$(basename "$l7_file")
    prefix=${base%_L7_1.fq.gz}

    r1_l7="${sample_dir}/${prefix}_L7_1.fq.gz"
    r1_l8="${sample_dir}/${prefix}_L8_1.fq.gz"
    r2_l7="${sample_dir}/${prefix}_L7_2.fq.gz"
    r2_l8="${sample_dir}/${prefix}_L8_2.fq.gz"

    out_r1="${OUTDIR}/${prefix}_1.fq.gz"
    out_r2="${OUTDIR}/${prefix}_2.fq.gz"

    echo "Sample: ${sample_dir}" >> "$LOG"

    echo "  R1 inputs:" >> "$LOG"
    echo "    ${r1_l7}" >> "$LOG"
    echo "    ${r1_l8}" >> "$LOG"
    echo "  R1 output: ${out_r1}" >> "$LOG"
    echo "  Command: cat ${r1_l7} ${r1_l8} > ${out_r1}" >> "$LOG"

    echo "  R2 inputs:" >> "$LOG"
    echo "    ${r2_l7}" >> "$LOG"
    echo "    ${r2_l8}" >> "$LOG"
    echo "  R2 output: ${out_r2}" >> "$LOG"
    echo "  Command: cat ${r2_l7} ${r2_l8} > ${out_r2}" >> "$LOG"

    cat "$r1_l7" "$r1_l8" > "$out_r1"
    cat "$r2_l7" "$r2_l8" > "$out_r2"

    echo "  Status: OK" >> "$LOG"
    echo "" >> "$LOG"
done
