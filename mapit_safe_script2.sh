#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# --- Configuration ---
GENOME_NAME="GRCh38"
FASTQ_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/example/fastq_samples"
MAPIT_BIN="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/Mapit"
OUTPUT_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT_output_safe"
CONFIG_PATH="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/conf/GRCh38.json"
REPLICATE="R1"
RNA_STRAND="FR"
LOGFILE="${OUTPUT_DIR}/mapit_run.log"

mkdir -p "$OUTPUT_DIR"

# Ensure write permissions to any pre-existing BED files
find "$OUTPUT_DIR" -type f -name "*.bed" -exec chmod u+w {} \;

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

log "Starting MAPIT-safe full multi-sample pipeline..."

# --- Validate MAPIT binary ---
if [[ ! -x "$MAPIT_BIN" ]]; then
    log "ERROR: MAPIT binary not found or not executable: $MAPIT_BIN"
    exit 1
fi

# --- Validate config ---
if [[ ! -f "$CONFIG_PATH" ]]; then
    log "ERROR: Config file does not exist: $CONFIG_PATH"
    exit 1
fi

# --- Step 1: Prepare genome & annotation for each sample ---
for fq in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$fq" .fastq)
    SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE_NAME"
    mkdir -p "$SAMPLE_OUT"
    log "Preparing sample: $SAMPLE_NAME"
    "$MAPIT_BIN" prepare \
        -v "$GENOME_NAME" \
        -n "$SAMPLE_NAME" \
        -r "$REPLICATE" \
        -o "$SAMPLE_OUT" \
        2>&1 | tee -a "$LOGFILE"
done

# --- Step 2: Map sequencing reads ---
for fq in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$fq" .fastq)
    SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE_NAME"
    log "Mapping sample: $SAMPLE_NAME"
    "$MAPIT_BIN" mapping \
        -v "$GENOME_NAME" \
        -n "$SAMPLE_NAME" \
        -r "$REPLICATE" \
        --fq "$fq" \
        --rna-strandness "$RNA_STRAND" \
        -o "$SAMPLE_OUT" \
        2>&1 | tee -a "$LOGFILE"
done

# --- Step 3: Call RNA editing (multi-sample support) ---
SAMPLE_LIST="$OUTPUT_DIR/sample_list.txt"
ls "$FASTQ_DIR"/*.fastq | xargs -n1 basename | sed 's/.fastq//' > "$SAMPLE_LIST"
log "Calling RNA editing on all samples..."
"$MAPIT_BIN" callediting \
    -v "$GENOME_NAME" \
    --sampleList "$SAMPLE_LIST" \
    -o "$OUTPUT_DIR" \
    --prefix "MAPIT" \
    -e Both \
    2>&1 | tee -a "$LOGFILE"

# --- Step 4: Call targets ---
log "Calling targets..."
"$MAPIT_BIN" calltargets \
    -v "$GENOME_NAME" \
    -i "$OUTPUT_DIR/MAPIT_edits.tsv" \
    -l Both \
    --treatName TREAT \
    --controlName CONTROL \
    -o "$OUTPUT_DIR" \
    2>&1 | tee -a "$LOGFILE"

# --- Step 5: Run FLARE analysis per sample ---
for fq in "$FASTQ_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$fq" .fastq)
    SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE_NAME"
    log "Running FLARE for sample: $SAMPLE_NAME"
    "$MAPIT_BIN" FLARE \
        -v "$GENOME_NAME" \
        -n "$SAMPLE_NAME" \
        -r "$REPLICATE" \
        --regions "$OUTPUT_DIR/MAPIT_targets.bed" \
        -o "$SAMPLE_OUT" \
        -e AG \
        2>&1 | tee -a "$LOGFILE"
done

# --- Step 6: Ensure BED permissions permanently ---
find "$OUTPUT_DIR" -type f -name "*.bed" -exec chmod u+rw {} \;

log "MAPIT-safe full multi-sample pipeline complete. Output in $OUTPUT_DIR"

