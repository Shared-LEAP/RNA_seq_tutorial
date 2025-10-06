#!/bin/bash
#SBATCH --job-name=rna_mouse_unzip
#SBATCH --output="slurmLog-unzip.log"
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G


# Step 1: Unzip gz files
    echo "Step 1: Decompressing files..."
    source activate ../anaconda3/envs/gzip
    
    gzip -dk <Your Sample RNA Files>.fastq.gz
    gzip -dk <Your Sample RNA Files>.fastq.gz
    gzip -dk /home/jnavarrete/anagha_mouse/raw_reads/Z671DMSOR3ns_combined.fastq.gz
    gzip -dk /home/jnavarrete/anagha_mouse/raw_reads/Z671DOXR1ns_combined.fastq.gz
    gzip -dk /home/jnavarrete/anagha_mouse/raw_reads/Z671DOXR2ns_combined.fastq.gz
    gzip -dk /home/jnavarrete/anagha_mouse/raw_reads/Z671DOXR3ns_combined.fastq.gz
    
    conda deactivate
    echo "Decompressed files"


# Step 2: Initial Quality Control with FastQC
    echo "Step 2: Initial Quality Control..."
    # Create the directory first
    mkdir -p FastQC_reads
    
    source activate ../anaconda3/envs/fastqc
    
    fastqc /home/jnavarrete/anagha_mouse/unzipped_reads/*.fastq -o FastQC_reads
    
    
    conda deactivate 
    echo "FastQC completed"
    
    mkdir -p Trimmed_reads


# Step 3: Trim Reads using trimmomatic

    fastq_files=(/home/jnavarrete/anagha_mouse/unzipped_reads/*.fastq)
    source activate ../anaconda3/envs/trimmomatic
    
    for fastq_file in "${fastq_files[@]}"; do
        sample_name=$(basename "$fastq_file" _combined.fastq)
        
        echo "Processing: $sample_name"
        
         trimmomatic SE \
            -threads 4 \
            "$fastq_file" \
            "Trimmed_reads/${sample_name}_trimmed.fastq" \
            ILLUMINACLIP:/home/jnavarrete/anagha_mouse/TruSeq3-SE.fa:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36
    done
    conda deactivate
    echo "Trimming completed"


# Step 5: Run rsem for gene expression from RNA-Seq dataset
# Summary:
    #rsem quantifies gene and isoform expression levels from RNA-Seq data 

    #rsem-prepare-reference: generates a set of reference transcript sequences 
    #rsem-calculate-expression: aligns RNA-seq reads to the reference transcripts and uses the alignments to estimate
    #abundances.

    #rsem output: transcript-level and a gene-level count estimate

    echo "Step 5: Running RSEM for gene expression quantification..."
    echo "Getting Sweaty.."
    # Create directories for RSEM outputs
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
