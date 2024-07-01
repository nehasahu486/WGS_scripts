#!/bin/bash
#SBATCH -p tsl-long                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,fail              # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address


source package /tgac/software/production/bin/samtools-0.1.19
source package /tgac/software/production/bin/bowtie-2.2.5


# Specify the path to the reference genome
# Define variables
reference_genome="/location/of/reference_genome_fasta/MGGv8_genome.fasta"
index_prefix="MGGv8_index"
input_dir="/location/of/input/directory"
# Create index for reference genome
bowtie2-build $reference_genome $index_prefix

# Loop through all pairs of reads in the specified input directory
for forward_read in ${input_dir}/*_1.fq.gz; do
    # Extract sample name from filename
    sample_name=$(basename $forward_read _1.fq.gz)

    # Define filenames for forward and reverse reads
    reverse_read="${input_dir}/${sample_name}_2.fq.gz"
    sam_file="${sample_name}.sam"
    bam_file="${sample_name}.bam"
    sorted_bam_file="${sample_name}_sorted.bam"
    indexed_sorted_bam_file="${sample_name}_sorted.bam.bai"

    # Align reads to reference genome using zcat to decompress on-the-fly
    bowtie2 -x $index_prefix -1 <(zcat $forward_read) -2 <(zcat $reverse_read) -S $sam_file

    # Convert SAM to BAM
    samtools view -@ 16 -Sb -o $bam_file $sam_file

    # Sort BAM by genomic coordinates
    samtools sort -@ 16 $bam_file -o $sorted_bam_file
    
    # Index sorted BAM file
    samtools index $sorted_bam_file
done


