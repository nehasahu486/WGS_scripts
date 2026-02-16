#!/bin/bash
#SBATCH -p tsl-medium
#SBATCH --wckey=talbot_core
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --mem 100000
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk

source package /tgac/software/production/bin/bowtie-2.2.5
source package c92263ec-95e5-43eb-a527-8f1496d56f1a  # Samtools 1.18

# Setup
reference_genome="MGGv8_genome.fasta"
index_prefix="MGGv8_index"
input_dir1="/tsl/data/reads//ntalbot/wgs_bip1_ko_mgg_08118_guy11/wgs_bip1_ko_mgg_08118_guy11/_bip1_ko_t24_guy11/" # transformant_24
input_dir2="/tsl/data/reads//ntalbot/wgs_bip1_ko_mgg_08118_guy11/wgs_bip1_ko_mgg_08118_guy11/bip1_ko_t13_guy11/" # transformant_13
output_dir="./alignments_wgs"

mkdir -p $output_dir

echo "========================================="
echo "Starting alignment pipeline with FLAG filtering"
echo "========================================="

# Build index if needed
if [ ! -f "${index_prefix}.1.bt2" ]; then
    echo "Building bowtie2 index..."
    bowtie2-build $reference_genome $index_prefix
fi

# Knockout gene regions
declare -A GENES
GENES[mst12]="3:4220800-4225271"
GENES[bip1]="2:5750363-5753268"



find "$input_dir1" "$input_dir2" -type f \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz" \) | while read -r forward_read; do
    current=$((current + 1))
    
    # Determine file extension and sample name
    if [[ "$forward_read" == *"_1.fq.gz" ]]; then
        sample_name=$(basename "$forward_read" _1.fq.gz)
        reverse_read="$(dirname "$forward_read")/${sample_name}_2.fq.gz"
    else
        sample_name=$(basename "$forward_read" _1.fastq.gz)
        reverse_read="$(dirname "$forward_read")/${sample_name}_2.fastq.gz"
    fi
    
    echo "========================================="
    echo "[$current/$total] Processing: $sample_name"
    echo "========================================="
    
    if [ ! -f "$reverse_read" ]; then
        echo "ERROR: Reverse read not found"
        continue
    fi
    
    sam_file="${output_dir}/${sample_name}.sam"
    
    # Define output BAM files for different flag categories
    bam_all="${output_dir}/${sample_name}_all.bam"
    bam_unique="${output_dir}/${sample_name}_unique.bam"
    bam_paired="${output_dir}/${sample_name}_properly_paired.bam"
    bam_unique_paired="${output_dir}/${sample_name}_unique_paired.bam"
    bam_first="${output_dir}/${sample_name}_first_in_pair.bam"
    bam_second="${output_dir}/${sample_name}_second_in_pair.bam"
    bam_unmapped="${output_dir}/${sample_name}_unmapped.bam"
    bam_secondary="${output_dir}/${sample_name}_secondary.bam"
    bam_supplementary="${output_dir}/${sample_name}_supplementary.bam"
    
    # Align
    echo "[$(date +%H:%M:%S)] Running bowtie2..."
    bowtie2 -p 28 -x $index_prefix \
            -1 <(zcat "$forward_read") \
            -2 <(zcat "$reverse_read") \
            -S "$sam_file"
    
    if [ ! -s "$sam_file" ]; then
        echo "ERROR: SAM file not created"
        continue
    fi
    echo "SAM created: $(du -h $sam_file | cut -f1)"
    
    echo "[$(date +%H:%M:%S)] Creating filtered BAM files..."
    
    # 1. ALL READS (no filter)
    echo "  - All reads"
    samtools view -@ 28 -bS "$sam_file" | samtools sort -@ 28 -o "$bam_all" -
    samtools index "$bam_all"
    samtools flagstat "$bam_all" > "${output_dir}/${sample_name}_all_stats.txt"
    
    # 2. UNIQUE READS (MAPQ >= 10, primary alignments only)
    # Excludes: multimappers (MAPQ<10), unmapped (0x4), secondary (0x100), supplementary (0x800)
    echo "  - Unique reads (MAPQ>=10, primary only)"
    samtools view -@ 28 -bS -q 10 -F 4 -F 256 -F 2048 "$sam_file" | samtools sort -@ 28 -o "$bam_unique" -
    if [ -s "$bam_unique" ]; then
        samtools index "$bam_unique"
        samtools flagstat "$bam_unique" > "${output_dir}/${sample_name}_unique_stats.txt"
    fi
    
    # 3. PROPERLY PAIRED READS (0x2 = read mapped in proper pair)
    echo "  - Properly paired reads"
    samtools view -@ 28 -bS -f 2 "$sam_file" | samtools sort -@ 28 -o "$bam_paired" -
    if [ -s "$bam_paired" ]; then
        samtools index "$bam_paired"
        samtools flagstat "$bam_paired" > "${output_dir}/${sample_name}_paired_stats.txt"
    fi
    
    # 4. UNIQUE + PROPERLY PAIRED (most stringent)
    # Combines: uniquely mapped (MAPQ>=10) AND properly paired AND primary only
    echo "  - Unique AND properly paired (MOST STRINGENT)"
    samtools view -@ 28 -bS -q 10 -f 2 -F 4 -F 256 -F 2048 "$sam_file" | samtools sort -@ 28 -o "$bam_unique_paired" -
    if [ -s "$bam_unique_paired" ]; then
        samtools index "$bam_unique_paired"
        samtools flagstat "$bam_unique_paired" > "${output_dir}/${sample_name}_unique_paired_stats.txt"
    fi
    
    # 5. FIRST IN PAIR (0x40 = first read in pair)
    echo "  - First in pair"
    samtools view -@ 28 -bS -f 64 "$sam_file" | samtools sort -@ 28 -o "$bam_first" -
    if [ -s "$bam_first" ]; then
        samtools index "$bam_first"
    fi
    
    # 6. SECOND IN PAIR (0x80 = second read in pair)
    echo "  - Second in pair"
    samtools view -@ 28 -bS -f 128 "$sam_file" | samtools sort -@ 28 -o "$bam_second" -
    if [ -s "$bam_second" ]; then
        samtools index "$bam_second"
    fi
    
    # 7. UNMAPPED READS (0x4 = read unmapped)
    echo "  - Unmapped reads"
    samtools view -@ 28 -bS -f 4 "$sam_file" | samtools sort -@ 28 -o "$bam_unmapped" -
    if [ -s "$bam_unmapped" ]; then
        samtools index "$bam_unmapped"
    fi
    
    # 8. SECONDARY ALIGNMENTS (0x100 = not primary alignment)
    echo "  - Secondary alignments"
    samtools view -@ 28 -bS -f 256 "$sam_file" | samtools sort -@ 28 -o "$bam_secondary" -
    if [ -s "$bam_secondary" ]; then
        samtools index "$bam_secondary"
    fi
    
    # 9. SUPPLEMENTARY ALIGNMENTS (0x800 = supplementary alignment)
    echo "  - Supplementary alignments"
    samtools view -@ 28 -bS -f 2048 "$sam_file" | samtools sort -@ 28 -o "$bam_supplementary" -
    if [ -s "$bam_supplementary" ]; then
        samtools index "$bam_supplementary"
    fi
    
    # Knockout analysis on key BAM files
    echo ""
    echo "Knockout region analysis:"
    for gene in "${!GENES[@]}"; do
        region="${GENES[$gene]}"
        all_count=$(samtools view -c "$bam_all" "$region" 2>/dev/null || echo 0)
        unique_count=$(samtools view -c "$bam_unique" "$region" 2>/dev/null || echo 0)
        paired_count=$(samtools view -c "$bam_paired" "$region" 2>/dev/null || echo 0)
        unique_paired_count=$(samtools view -c "$bam_unique_paired" "$region" 2>/dev/null || echo 0)
        secondary_count=$(samtools view -c "$bam_secondary" "$region" 2>/dev/null || echo 0)
        
        echo "  ${gene}:"
        echo "    ALL: ${all_count} | UNIQUE: ${unique_count} | PAIRED: ${paired_count} | UNIQUE+PAIRED: ${unique_paired_count} | SECONDARY: ${secondary_count}"
    done
    
    # Cleanup
    echo "[$(date +%H:%M:%S)] Cleaning up..."
    rm -f "$sam_file"
    
    echo "Done with $sample_name"
    echo ""
done