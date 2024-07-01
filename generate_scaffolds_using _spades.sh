#!/bin/bash
#SBATCH -p tsl-long                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,fail              # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address


 
# Load spades (SPAdes - 3.13.1)
source package 3a579940-1cba-4d60-86ef-43b8705935fb

# Create output directory
mkdir -p output

# Run SPAdes assembly using the .fq files from Novogene or any other service provider
spades.py -1 *1.fq -2 *2.fq -o output/

#The output directory will have a file called scaffolds.fasta. This is the assembly from .fq files.
