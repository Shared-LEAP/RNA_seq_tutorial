# RNA-seq Data Analysis Start to Finish

Tool: RNA-seq Analysis

Author: Josue Navarrete

Date: 09/25/2025

Status: In-progress

Summary: The following are steps to obtain, process, and analyze RNAseq data for differential gene expression analysis. 

**Note:
Make sure to change directories to relative paths in your system. Make appropriate folder for project**

_________________________________________________________________________________________________________________________________________________________

Run on cluster (personal computer not sufficient) 

_________________________________________________________________________________________________________________________________________________________

Installation Requirements
1. gzip
2. fastqc 
3. trimmomatic 
4. RSEM + STAR (i believe if u install rsem it comes with star)
5. Deseq2 for differential gene expression


Create conda environments 
```
# 1. Create gzip envirnoment
        conda create -n gzip -c conda-forge
# 1A. Install gzip in the gzip environment
        conda activate gzip
        conda install -c conda-forge gzip -y
        conda deactivate

# 2. Create fastqc envirnoment
        conda create -n fastqc fastqc
# 2A. Install fastqc in the fastqc environment
        conda activate fastqc
        conda install -c bioconda fastqc -y
        conda deactivate

# 3. Create trimmomatic enviroment 
        conda create -n trimmomatic -c bioconda -y
# 3A. Install trimmomatic in the fastqc environment
        conda activate trimmomatic
        conda install -c bioconda trimmomatic -y
        conda deactivate

# 4. Create rsem environment
        conda create -n rsem -c conda-forge -c bioconda rsem star -y
# 4A. Install rsem in the rsem environment
        conda activate rsem
        conda install -c conda-forge -c bioconda rsem star -y
        conda deactivate
```


# Use Any long read RNA Samples
## The following are general steps for processing Any lond read RNA Data:
-----------------------------------------------------------------------------
Step 1: Unzip gz 
    conda activate gzip
    
    gzip -dk <path to your folder>/Raw_reads/ENCFF489WGO.fastq.gz.download/ENCFF489WGO.fastq.gz

    gzip -dk <path to your folder>/Raw_reads/ENCFF522FSF.fastq.gz.download/ENCFF522FSF.fastq.gz 

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
    <path to your folder>/Raw_reads/ENCFF489WGO.fastq \
    <path to your folder>/Raw_reads/ENCFF522FSF.fastq \
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

### Create directories for RSEM outputs
        mkdir -p RSEM_reference
        mkdir -p RSEM_results
        
        source activate ../anaconda3/envs/rsem

# Define reference files using your specific paths
        GENOME_FASTA="/home/jnavarrete/anagha_mouse/rsem_reference/GRCm39.primary_assembly.genome.fa"
        GTF_FILE="/home/jnavarrete/anagha_mouse/rsem_reference/gencode.vM31.primary_assembly.annotation.gtf"

# Step 5.1: Prepare reference
        echo "Step 5.1: Preparing RSEM reference..."
        rsem-prepare-reference \
            --gtf "$GTF_FILE" \
            --star \
            -p 10 \
            "$GENOME_FASTA" \
            /home/jnavarrete/anagha_mouse/RSEM_reference/mouse_reference

# Check if reference preparation was successful
        if [ $? -eq 0 ]; then
            echo "Reference preparation completed successfully!"
        else
            echo "ERROR: Reference preparation failed! Please check the error messages above."
            exit 1
        fi
        
        conda deactivate
        
        
        source activate ../anaconda3/envs/rsem

# Step 5.2: Calculate expression for each sample
        echo "Step 5.2: Calculating expression for each sample..."
        
        trimmed_files=(/home/jnavarrete/anagha_mouse/Trimmed_reads/*_trimmed.fastq)
        
        for trimmed_file in "${trimmed_files[@]}"; do
            sample_name=$(basename "$trimmed_file" _trimmed.fastq)
            
            echo "Processing sample: $sample_name"
            
            rsem-calculate-expression \
                --star \
                -p 10 \
                --estimate-rspd \
                --output-genome-bam \
                "$trimmed_file" \
                /home/jnavarrete/anagha_mouse/RSEM_reference/mouse_reference \
                /home/jnavarrete/anagha_mouse/RSEM_results/"$sample_name"
        done
        
        conda deactivate

# Generate summary report
        echo "Generating RSEM summary report..."
        
        source activate ../anaconda3/envs/rsem
        rsem-generate-data-matrix RSEM_results/*.genes.results > RSEM_results/gene_expression_matrix.txt
        rsem-generate-data-matrix RSEM_results/*.isoforms.results > RSEM_results/isoform_expression_matrix.txt
        conda deactivate
        
        echo "RSEM analysis completed!"
        echo "Gene-level results: RSEM_results/gene_expression_matrix.txt"
        echo "Isoform-level results: RSEM_results/isoform_expression_matrix.txt"
-----------------------------------------------------------------------------

