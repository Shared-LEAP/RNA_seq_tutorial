#!/bin/bash
#SBATCH --job-name=setup_conda
#SBATCH --output=conda_setup_%j.out
#SBATCH --error=conda_setup_%j.err
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=normal

# Script to set up conda environments for RNA-seq pipeline

module load miniconda3

# Create conda environments
echo "Setting up conda environments..."

# gzip environment
conda create -n gzip gzip -y

# fastqc environment
conda create -n fastqc fastqc -y

# trimmomatic environment
conda create -n trimmomatic trimmomatic -y

# rsem environment (using conda-forge as priority to avoid conflicts)
conda create -n rsem -c conda-forge -c bioconda rsem star -y

echo "Conda environments created successfully!"
echo "Available environments:"
conda env list