#!/bin/bash

# ==========================================
# DESCRIPTION
# ==========================================

# This script takes BAM file as an input and expects txt files with CNV simulation instruction in "sim" subdirectory.
# It will look for indexed reference in a "reference" subdirectory, if it doesn't find it, it downloads and indexes it.
# Parameters from CONFIGURATION must be double-checked before launching.

# Official github with theoretical info:
# https://github.com/adamewing/bamsurgeon

# ==========================================
# INSTALLATION (MUST BE DONE ONLY ONCE)
# ==========================================

# set correct priority of conda channels, this must be done only once on a system:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# create "bamsurgeon_env" environment for bamsurgeon and its dependancies only:
# conda create -n bamsurgeon_env bamsurgeon -y

# ==========================================
# CONFIGURATION
# ==========================================

THREADS=4

INPUT_BAM="wgs_30x_hs37d5.bam"

# Reference genome: "hg19" or "hg38" or "hs37d5"
REF="hs37d5"

# Set CHR to "all" for WGS, or a specific number inside quotes (e.g., "20") for single chromosome mode
CHR="all"

# ==========================================
# 0. DEPENDENCY CHECK
# ==========================================

echo "=========================================="
echo "CHECKING DEPENDENCIES"
echo "=========================================="

# Check for BamSurgeon (addsv.py)
if ! command -v addsv.py &> /dev/null; then
    echo "❌ ERROR: 'addsv.py' (BamSurgeon) not found."
    echo "   Please set up proper conda environment (e.g., 'conda activate bamsurgeon_env')."
    exit 1
else
    echo "✅ Found BamSurgeon (addsv.py)"
fi

# Check for BWA
if ! command -v bwa &> /dev/null; then
    echo "❌ ERROR: 'bwa' not found."
    echo "   Please install BWA (e.g., 'conda install bwa')."
    exit 1
else
    # bwa doesn't have a simple --version flag that outputs clean text, 
    # but it prints version info to stderr when run with no args or invalid args.
    # We grep it out safely.
    BWA_VER=$(bwa 2>&1 | grep -i "Version" | head -n 1)
    echo "✅ Found BWA: $BWA_VER"
fi

# Check for Samtools
if ! command -v samtools &> /dev/null; then
    echo "❌ ERROR: 'samtools' not found."
    echo "   Please install Samtools (e.g., 'conda install samtools')."
    exit 1
else
    SAM_VER=$(samtools --version | head -n 1)
    echo "✅ Found Samtools: $SAM_VER"
fi

# ==========================================
# 1. PRE-FLIGHT CHECKS
# ==========================================

SIM_DIR="./sim"

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM '$INPUT_BAM' not found."
    exit 1
fi

if [ ! -d "$SIM_DIR" ]; then
    echo "Error: Directory '$SIM_DIR' does not exist."
    exit 1
fi

shopt -s nullglob
TXT_FILES=("$SIM_DIR"/*.txt)
shopt -u nullglob

if [ ${#TXT_FILES[@]} -eq 0 ]; then
    echo "Error: No .txt files found in '$SIM_DIR'. Exiting."
    exit 1
fi

echo "Found ${#TXT_FILES[@]} variant list(s)."

# ==========================================
# 2. UNIVERSAL REFERENCE SETUP
# ==========================================

REF_DIR="./reference"
mkdir -p "$REF_DIR"

# 2a. Detect BAM Naming Convention (For logging only)
if [ -f "$INPUT_BAM" ]; then
    BAM_FIRST_CHR=$(samtools view -H "$INPUT_BAM" | grep "^@SQ" | head -n 1 | sed 's/.*SN:\([^\t]*\).*/\1/')
    echo "Detected BAM Header First Chromosome: '$BAM_FIRST_CHR'"
fi

# 2b. Define Reference Paths & URLs
if [ "$REF" == "hs37d5" ]; then
    # --- hs37d5 Logic ---
    REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    REF_CLEAN="${REF_DIR}/hs37d5.fa"
    REF_RAW="${REF_DIR}/hs37d5.fa.gz"

    if [ "$CHR" != "all" ]; then
        echo "⚠️  WARNING: Single chromosome mode ($CHR) requested for 'hs37d5'."
        echo "    'hs37d5' is not available as split chromosomes."
        echo "    Downloading and indexing the FULL genome reference instead."
    fi

else
    # --- hg19 / hg38 Logic ---
    if [ "$CHR" == "all" ]; then
        # WGS
        REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF}/bigZips/${REF}.fa.gz"
        REF_CLEAN="${REF_DIR}/${REF}.fa"
        REF_RAW="${REF_DIR}/${REF}.fa.gz"
    else
        # Single Chromosome
        REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF}/chromosomes/chr${CHR}.fa.gz"
        REF_CLEAN="${REF_DIR}/chr${CHR}.fa"
        REF_RAW="${REF_DIR}/chr${CHR}.fa.gz"
    fi
fi

# 2c. Download & Decompress
if [ ! -f "$REF_CLEAN" ]; then
    echo "Downloading Reference from: $REF_URL"
    
    # Clear old indices to prevent mismatches
    rm -f "${REF_CLEAN}.bwt" "${REF_CLEAN}.fai" "${REF_CLEAN}.pac" "${REF_CLEAN}.ann" "${REF_CLEAN}.amb" "${REF_CLEAN}.sa"
    
    wget -q -O "$REF_RAW" "$REF_URL"
    
    echo "Decompressing Reference..."
    # Simply decompress the file, treating the source headers as the source of truth
    zcat "$REF_RAW" > "$REF_CLEAN"
    
    rm "$REF_RAW"
else
    echo "Reference already exists: $REF_CLEAN"
fi

# 2d. Indexing (Only if missing)
if [ ! -f "${REF_CLEAN}.fai" ]; then
    echo "Indexing for Samtools (fai)..."
    samtools faidx "$REF_CLEAN"
fi

if [ ! -f "${REF_CLEAN}.bwt" ]; then
    echo "Indexing for BWA..."
    bwa index "$REF_CLEAN"
fi

# ==========================================
# 3. BATCH EXECUTION
# ==========================================

for VAR_LIST in "${TXT_FILES[@]}"; do
    
    BASE_INPUT=$(basename "$INPUT_BAM" .bam)
    BASE_VAR=$(basename "$VAR_LIST" .txt)
    OUTPUT_BAM="${BASE_INPUT}_${BASE_VAR}.bam"
    
    # Create isolated temp directory for each run
    TEMP_DIR="bamsurgeon_tmp_${BASE_VAR}"
    rm -rf "$TEMP_DIR" 
    mkdir -p "$TEMP_DIR"
    TEMP_UNSORTED="${TEMP_DIR}/unsorted.bam"

    echo "------------------------------------------------"
    echo "Processing: $VAR_LIST"
    echo "------------------------------------------------"

    # Clean legacy temp file if present
    rm -rf addsv.tmp 

    addsv.py \
      -v "$VAR_LIST" \
      -f "$INPUT_BAM" \
      -r "$REF_CLEAN" \
      -o "$TEMP_UNSORTED" \
      --aligner mem \
      -p "$THREADS" \
      --alignerthreads "$THREADS" \
      --mindepth 1

    if [ ! -f "$TEMP_UNSORTED" ]; then
        echo "Error: BamSurgeon failed for $VAR_LIST"
        continue
    fi

    echo "Sorting and indexing..."
    samtools sort -@ "$THREADS" -o "$OUTPUT_BAM" "$TEMP_UNSORTED"
    samtools index "$OUTPUT_BAM"

    # Cleanup
    rm -rf "$TEMP_DIR"
    rm -rf addsv.tmp
    rm "$TEMP_UNSORTED" 2>/dev/null
    
    echo "Finished: $OUTPUT_BAM"

done

echo "Batch processing complete."