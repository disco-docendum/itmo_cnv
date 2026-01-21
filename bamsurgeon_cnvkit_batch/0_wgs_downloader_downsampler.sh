# ==========================================
# DESCRIPTION
# ==========================================

# This script downloads and downsamples HG002 WGS data with 60x depth from NCBI server of Genome in a Bottle project.
# Streaming is not viable because single connection break will mean starting all over.
# Parameters from CONFIGURATION must be double-checked before launching.

# ==========================================
# INSTALLATION (MUST BE DONE ONLY ONCE)
# ==========================================

# set correct priority of conda channels, this must be done only once on a system:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# create "samtools_env" environment for samtools:
# conda create -n samtools_env samtools -y

# ==========================================
# CONFIGURATION
# ==========================================

THREADS=4

#Downloading BAM and BAI:
URL_BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai"
DOWNLOAD_BAM="wgs_60x_hs37d5.bam.bai"

URL_BAI="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai"
DOWNLOAD_BAI="wgs_60x_hs37d5.bam.bai"

#Name of downsampled BAM file:
OUTPUT_BAM="wgs_30x_hs37d5.bam"

#Downsampling variables:
SEED=42
DOWNSAMPING_RATE=0.5

# ==========================================
# 0. DEPENDENCY CHECK
# ==========================================

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
# 1. DOWNLOADING
# ==========================================

# The robust download commands with retry
echo "Starting resilient download..."
until curl -L -C - --retry 20 --retry-delay 5 --limit-rate 6m -o "$DOWNLOAD_BAM" "$URL_BAM"; do
    echo "Check: Download interrupted. Sleeping 5s before resuming..."
    sleep 5
done
echo "BAM download complete!"

echo "Starting resilient download..."
until curl -L -C - --retry 20 --retry-delay 5 --limit-rate 6m -o "$DOWNLOAD_BAI" "$URL_BAI"; do
    echo "Check: Download interrupted. Sleeping 5s before resuming..."
    sleep 5
done
echo "BAI download complete!"

# ==========================================
# 2. DOWNSAMPLING AND FILTERING
# ==========================================

SUBSAMPLE_ARG=$(awk "BEGIN {print $SEED + $DOWNSAMPLE_FRACTION}")

echo "--------------------------------------"
echo "Starting downsampling..."
echo "  Seed:       $SEED"
echo "  Fraction:   $DOWNSAMPING_RATE"
echo "--------------------------------------"

samtools view -h -f 2 -F 256 "$DOWNLOAD_BAM" | \
samtools view -b -s "$SUBSAMPLE_ARG" -@ "$THREADS" - > "$OUTPUT_BAM"

echo "Starting building index..."
samtools index "$OUTPUT_BAM"

echo "Building stats report file..."
samtools stats -@ "$THREADS" "$OUTPUT_BAM" > "$OUTPUT_BAM.stats"