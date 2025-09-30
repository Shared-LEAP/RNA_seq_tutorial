#!/bin/bash
#SBATCH --job-name=rna_expression
#SBATCH --output="slurmLog-rna_expression.txt"
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --ntasks-per-node=50

# RNA-seq Pipeline SLURM Script
# This script performs: quality control, trimming, and gene expression quantification

# Load required modules (modify based on your cluster's available modules)
module load miniconda3

# Set variables
RAW_READS_DIR="Raw_reads"
TRIMMED_DIR="Trim_out"
RSEM_REF_DIR="rsem_reference"
QC_DIR="fastqc_results"
SAMPLE1="ENCFF489WGO"
SAMPLE2="ENCFF522FSF"

# Create directories
mkdir -p ${TRIMMED_DIR} ${QC_DIR} ${RSEM_REF_DIR}

echo "Starting RNA-seq pipeline at: $(date)"

# Step 1: Unzip gz files
echo "Step 1: Unzipping files..."
conda activate gzip

if [ ! -f "${RAW_READS_DIR}/${SAMPLE1}.fastq" ]; then
    gzip -dk "${RAW_READS_DIR}/${SAMPLE1}.fastq.gz"
else
    echo "${SAMPLE1}.fastq already exists, skipping unzip"
fi

if [ ! -f "${RAW_READS_DIR}/${SAMPLE2}.fastq" ]; then
    gzip -dk "${RAW_READS_DIR}/${SAMPLE2}.fastq.gz"
else
    echo "${SAMPLE2}.fastq already exists, skipping unzip"
fi

conda deactivate

# Step 2: Initial Quality Control with FastQC
echo "Step 2: Initial Quality Control..."
conda activate fastqc

fastqc "${RAW_READS_DIR}/${SAMPLE1}.fastq" "${RAW_READS_DIR}/${SAMPLE2}.fastq" -o ${QC_DIR}

conda deactivate

# Step 3: Trim Reads using Trimmomatic
echo "Step 3: Trimming reads..."
conda activate trimmomatic

# Check if Trimmomatic adapter file exists, if not download it
if [ ! -f "TruSeq3-PE.fa" ]; then
    echo "Downloading Trimmomatic adapter file..."
    wget https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE.fa
fi

trimmomatic PE \
    -threads ${SLURM_CPUS_PER_TASK} \
    "${RAW_READS_DIR}/${SAMPLE1}.fastq" \
    "${RAW_READS_DIR}/${SAMPLE2}.fastq" \
    "${TRIMMED_DIR}/output_forward_paired.fastq" \
    "${TRIMMED_DIR}/output_forward_unpaired.fastq" \
    "${TRIMMED_DIR}/output_reverse_paired.fastq" \
    "${TRIMMED_DIR}/output_reverse_unpaired.fastq" \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

conda deactivate

# Step 4: Post-trimming Quality Control
echo "Step 4: Post-trimming Quality Control..."
conda activate fastqc

fastqc \
    "${TRIMMED_DIR}/output_forward_paired.fastq" \
    "${TRIMMED_DIR}/output_forward_unpaired.fastq" \
    "${TRIMMED_DIR}/output_reverse_paired.fastq" \
    "${TRIMMED_DIR}/output_reverse_unpaired.fastq" \
    -o ${QC_DIR}

conda deactivate

# Step 5: RSEM Reference Preparation (only if not already done)
echo "Step 5: Preparing RSEM reference..."
cd ${RSEM_REF_DIR}

# Check if reference files exist, if not download them
if [ ! -f "GRCh38.primary_assembly.genome.fa" ]; then
    echo "Downloading reference genome..."
    curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz
fi

if [ ! -f "gencode.v48.primary_assembly.annotation.gtf" ]; then
    echo "Downloading annotation file..."
    curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
    gunzip gencode.v48.primary_assembly.annotation.gtf.gz
fi

# Build RSEM reference if it doesn't exist
if [ ! -d "hg38_v48.transcripts.fa" ]; then
    echo "Building RSEM reference..."
    conda activate rsem
    
    rsem-prepare-reference \
        --gtf gencode.v48.primary_assembly.annotation.gtf \
        --star \
        -p ${SLURM_CPUS_PER_TASK} \
        GRCh38.primary_assembly.genome.fa \
        hg38_v48
    
    conda deactivate
else
    echo "RSEM reference already exists, skipping preparation"
fi

cd ..

# Step 6: Gene Expression Quantification with RSEM
echo "Step 6: Gene expression quantification..."
conda activate rsem

# Create output directory for RSEM results
mkdir -p rsem_results

rsem-calculate-expression \
    --paired-end \
    --star \
    --star-gzipped-read-file \
    -p ${SLURM_CPUS_PER_TASK} \
    --output-genome-bam \
    "${TRIMMED_DIR}/output_forward_paired.fastq" \
    "${TRIMMED_DIR}/output_reverse_paired.fastq" \
    "${RSEM_REF_DIR}/hg38_v48" \
    "rsem_results/${SAMPLE1}_${SAMPLE2}"

conda deactivate

echo "RNA-seq pipeline completed at: $(date)"
echo "Results are in: rsem_results/"
echo "Quality control reports are in: ${QC_DIR}/"