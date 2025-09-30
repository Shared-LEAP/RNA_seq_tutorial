# RNA-seq Data Analysis Start to Finish

Tool: RNA-seq Analysis

Author: Josue Navarrete

Date: 09/25/2025

Status: In-progress

Summary: The following are steps to obtain, process, and analyze RNAseq data for differential gene expression analysis. 

**Note:
Make sure to change directories to relative paths in your system. Make appropriate folder for project**

_________________________________________________________________________________________________________________________________________________________

Option 1:

Run on cluster (personal computer not sufficient) 

        1. conda_setup_script.sh 
                This will set up the required environment in your project folder

                Use the following in terminal:
                    sbatch conda_setup_script.sh

                Verify queue:
                    squeue
        
        2. rna_expression_pipeline.sh 
                This will run each step mentioned below for proprocess RNA raw reads, alignment, quantification, and differential analysis

                Use the following in terminal:
                    sbatch rna_expression_pipeline.sh

                Verify queue:
                    squeue


Option 2:
Follow the steps below as a tutorial.

_________________________________________________________________________________________________________________________________________________________

Installation Requirements
1. gzip
2. fastqc 
3. trimmomatic 
4. RSEM + STAR (i believe if u install rsem it comes with star)
5. Deseq2 for differential gene expression


Create conda environments 
```   
        conda create -n <library name> <library name>
        
        ex: 
            
            conda create -n fastqc fastqc
```
Sample/Data: grab the paired fastq files 


ENCFF489WGO.fastq.gz
ENCFF522FSF.fastq.gz

From:

https://www.encodeproject.org/experiments/ENCSR002JOW/#:~:text=Isogenic%20replicateLibrary-,Accession,-File%20typeRun
<br/>

Biosample Summary - Homo sapiens K562 genetically modified (deletion) using CRISPR targeting H. sapiens RNASEH2A




Code:
RNA data files:
ENCFF489WGO.fastq.gz
ENCFF522FSF.fastq.gz
-----------------------------------------------------------------------------
Step 1: Unzip gz 
    conda activate gzip

    gzip -dk /Users/jnavarrete/Desktop/multi_omic-practice/rna_seq/Raw_reads/ENCFF489WGO.fastq.gz.download/ENCFF489WGO.fastq.gz

    gzip -dk /Users/jnavarrete/Desktop/multi_omic-practice/rna_seq/Raw_reads/ENCFF522FSF.fastq.gz.download/ENCFF522FSF.fastq.gz 

    conda deactivate
-----------------------------------------------------------------------------
Step 2: Perform Quality Control using fastqc
    conda activate fastqc

    fastqc Raw_reads/ENCFF489WGO.fastq Raw_reads/ENCFF522FSF.fastq

    conda deactivate 
-----------------------------------------------------------------------------
Step 3: Trim Reads using trimmomatic
    conda activate trimmomatic

    trimmomatic PE \
    /Users/jnavarrete/Desktop/multi_omic-practice/rna_seq/Raw_reads/ENCFF489WGO.fastq \
    /Users/jnavarrete/Desktop/multi_omic-practice/rna_seq/Raw_reads/ENCFF522FSF.fastq \
    output_forward_paired.fastq output_forward_unpaired.fastq \
    output_reverse_paired.fastq output_reverse_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

    conda deactivate
-----------------------------------------------------------------------------
Step 4: Recheck Quality using fastqc 
    conda activate fastqc

    fastqc Trim_out/output_forward_unpaired.fastq Trim_out/output_reverse_unpaired.fastq

    conda deactivate 
-----------------------------------------------------------------------------
Step 5: Run rsem for gene expression from RNA-Seq dataset

Summary:
    rsem quantifies gene and isoform expression levels from RNA-Seq data 

    rsem-prepare-reference: generates a set of reference transcript sequences 
    rsem-calculate-expression: aligns RNA-seq reads to the reference transcripts and uses the alignments to estimate
    abundances.

    rsem output: transcript-level and a gene-level count estimate


### Conda installation trobleshooting
    (I had dependency confilct issue with conda)
    Conda environment a to install rsem
    - Delete rsem env 
    - Create new rsem conda env 
    - Use conda-forge as priority channel
        conda install -c conda-forge -c bioconda rsem
    - Worked!, 
    - Check using 
        rsem-calculate-expression --version

### Create and enter reference directory
    mkdir -p rsem_reference
    cd rsem_reference

### Download Refernece Genome
    curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
    curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz

### Uncompress
    gunzip *.gz

### Build RSEM reference (this will take a while)
```
rsem-prepare-reference \
    --gtf gencode.v48.primary_assembly.annotation.gtf \
    --star \
    -p 8 \
    GRCh38.primary_assembly.genome.fa \
    hg38_v48
```
### star wasnt install with rsem
```
conda install -c bioconda star
```
### Then re-run the same rsem-prepare-reference command

-----------------------------------------------------------------------------

