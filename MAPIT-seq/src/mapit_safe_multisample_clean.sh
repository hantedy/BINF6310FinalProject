#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# --- Configuration ---
GENOME_NAME="GRCh38"
REPLICATE=1  # Must be an integer for MAPIT
FASTQ_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/example/fastq_samples"
CONFIG_PATH="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/conf/GRCh38.json"
MAPIT_BIN="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/Mapit"
OUTPUT_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT_output_safe"
LOGFILE="${OUTPUT_DIR}/mapit_run.log"

mkdir -p "$OUTPUT_DIR"
chmod u+rwx "$OUTPUT_DIR"

# Ensure annotation BED file is writable
ANNOTATION_BED="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT_output/gencode.v40.chr_patch_hapl_scaff.annotation_exon.bed"
if [[ -f "$ANNOTATION_BED" ]]; then
    chmod u+w "$ANNOTATION_BED"
fi

# --- Logging function ---
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

# ============= MAPIT-seq Pipeline ===================================================================
log "Starting MAPIT-safe full multi-sample pipeline..."

# --- Step 1: Validate MAPIT binary ---
if [[ ! -x "$MAPIT_BIN" ]]; then
    log "ERROR: MAPIT binary not found or not executable: $MAPIT_BIN"
    exit 1
fi

# --- Step 2: Validate config ---
if [[ ! -f "$CONFIG_PATH" ]]; then
    log "ERROR: Config file does not exist: $CONFIG_PATH"
    exit 1
fi

# --- Step 3: Prepare & Map Each Sample ---
for fq in "$FASTQ_DIR"/*.fastq; do
    if [[ ! -f "$fq" ]]; then
        log "No FASTQ files found in $FASTQ_DIR, skipping."
        continue
    fi

    SAMPLE_NAME=$(basename "$fq" .fastq)
    log "Preparing sample: $SAMPLE_NAME (Replicate: $REPLICATE)"

    "$MAPIT_BIN" prepare \
        -v "$GENOME_NAME" \
        -n "$SAMPLE_NAME" \
        -r "$REPLICATE" \
        -o "$OUTPUT_DIR" 2>&1 | tee -a "$LOGFILE"

    log "Mapping sample: $SAMPLE_NAME"
    "$MAPIT_BIN" mapping \
        -v "$GENOME_NAME" \
        --fq "$fq" \
        -n "$SAMPLE_NAME" \
        -r "$REPLICATE" \
        -o "$OUTPUT_DIR" 2>&1 | tee -a "$LOGFILE"
done

# --- Step 4: RNA editing calling ---
log "Running MAPIT editing call step..."
"$MAPIT_BIN" callediting \
    -v "$GENOME_NAME" \
    --sampleList "$FASTQ_DIR/sample_list.txt" \
    -o "$OUTPUT_DIR" \
    --prefix "MAPIT" \
    -e Both 2>&1 | tee -a "$LOGFILE"

# --- Step 5: Target calling ---
log "Running MAPIT target calling step..."
"$MAPIT_BIN" calltargets \
    -v "$GENOME_NAME" \
    -i "$OUTPUT_DIR/MAPIT_edits.tsv" \
    -l Both \
    --treatName TREAT \
    --controlName CONTROL \
    -o "$OUTPUT_DIR" 2>&1 | tee -a "$LOGFILE"

log "MAPIT-safe pipeline complete. Output in $OUTPUT_DIR"

