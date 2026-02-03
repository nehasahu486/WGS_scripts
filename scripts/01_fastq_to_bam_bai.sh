#!/bin/bash
#SBATCH -p tsl-medium                           # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,fail              # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send email

# Load packages
source package /tgac/software/production/bin/bowtie-2.2.5 #Bowtie-2.2.5
source package c92263ec-95e5-43eb-a527-8f1496d56f1a #Samtools - 1.18

# DEFINE PATHS - Make sure to EDIT this accordingly!!!

# Reference genome path
reference_genome="/hpc-home/nep23rer/njt_lab_sequences/MGGv8_genome.fasta"
index_prefix="MGGv8_index"

# Input directory containing FASTQ files
input_dir="/tsl/data/reads/ntalbot/my_deletion_mutants_directory/"

# STEP 1: Create index for reference genome
bowtie2-build $reference_genome $index_prefix

# STEP 2: Loop through all forward reads in input_dir (and subdirectories if any)
find "$input_dir" -type f -name "*_1.fq.gz" | while read -r forward_read; do
    # Extract sample name from filename
    sample_name=$(basename "$forward_read" _1.fq.gz)

    # Define reverse read path (same folder as forward read)
    reverse_read="$(dirname "$forward_read")/${sample_name}_2.fq.gz"
    sam_file="${sample_name}.sam"
    bam_file="${sample_name}.bam"
    sorted_bam_file="${sample_name}_sorted.bam"

    # Align reads
    bowtie2 -x $index_prefix -1 <(zcat "$forward_read") -2 <(zcat "$reverse_read") -S "$sam_file"

    # Convert SAM to BAM
    samtools view -@ 16 -Sb -o "$bam_file" "$sam_file"

    # Sort BAM
    samtools sort -@ 16 "$bam_file" -o "$sorted_bam_file"

    # Index BAM
    samtools index "$sorted_bam_file"
done