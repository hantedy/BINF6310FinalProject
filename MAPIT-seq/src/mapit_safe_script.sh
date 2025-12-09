#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# --- Configuration ---
GENOME_NAME="GRCh38"
# Auto-detect sample name from FASTQ files
FASTQ_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/example/fastq_samples"
FIRST_FASTQ=$(ls "$FASTQ_DIR"/*.fastq 2>/dev/null | head -n 1)
if [[ -z "$FIRST_FASTQ" ]]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR" >&2
    exit 1
fi
SAMPLE_NAME="$(basename "$FIRST_FASTQ" .fastq)"
CONFIG_PATH="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/conf/GRCh38.json"
MAPIT_BIN="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT-seq/MAPIT-seq/Mapit"
OUTPUT_DIR="/Users/hanchoi/PycharmProjects/BINF6310FinalProject/MAPIT_output_safe"
LOGFILE="${OUTPUT_DIR}/mapit_run.log"

mkdir -p "$OUTPUT_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

log "Starting MAPIT-safe pipeline..."

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

log "Running MAPIT prepare step..."
"$MAPIT_BIN" prepare --config "$CONFIG_PATH" 2>&1 | tee -a "$LOGFILE"

log "Running MAPIT mapping step..."
"$MAPIT_BIN" mapping --config "$CONFIG_PATH" 2>&1 | tee -a "$LOGFILE"

log "Running MAPIT editing call step..."
"$MAPIT_BIN" callediting --config "$CONFIG_PATH" 2>&1 | tee -a "$LOGFILE"

log "Running MAPIT target calling step..."
"$MAPIT_BIN" calltargets --config "$CONFIG_PATH" 2>&1 | tee -a "$LOGFILE"

log "MAPIT-safe pipeline complete. Output in $OUTPUT_DIR"
