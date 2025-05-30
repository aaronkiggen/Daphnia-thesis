# # # # # # # # # # # # # # # # # # # # # # # # # # #
# This SLURM job script automates the RNA-seq data processing pipeline for a given species using a reference genome, annotation file and reads in SRA format. #
# it includes:

#SRA to FASTQ conversion using fasterq-dump.#

#Quality trimming and adapter removal with trim_galore.#

#HISAT2 genome indexing and alignment of paired-end reads.#

#BAM file processing, including filtering, sorting, and indexing.#

#(Optional) Quality control using Qualimap.#

#Gene-level quantification using featureCounts, with strandness based on species number.#
# # # # # # # # # # # # # # # # # # # # # # # # # # #



#!/bin/bash -l

#SBATCH--cluster = genius
#SBATCH--job - name = RNAseq_mapping_test
#SBATCH--nodes = 2
#SBATCH--ntasks - per - node = 8
#SBATCH--time = 24: 00: 00
#SBATCH - A lp_svbelleghem

# Function to print usage information
usage() {
  echo "Usage: $0 <species_number> <reference_genome.fa> <reference_annotation.gtf>"
  echo "Species numbers: 1-8 (strandness=2 for 1 & 3, 0 for others)"
  exit 1
}

# Validate arguments
if ["$#" - ne 3];
then
echo "Error: Incorrect number of arguments provided."
usage
fi

# Assign and validate species number
species = $1
if [
  [!"$species" = ~ ^ [1 - 9] $]
];
then
echo "Error: Invalid species number. Must be 1-9."
usage
fi

# Validate reference files
reference_fa = $2
reference_gtf = $3
if [!-f "$reference_fa"] || [!-f "$reference_gtf"];
then
echo "Error: Reference genome or GTF file not found."
exit 1
fi

# Input / output directories
input_dir = "/lustre1/scratch/354/vsc35429/metastudy_final/input/SRA"
output_dir = "/lustre1/scratch/354/vsc35429/metastudy_final/output"
species_dir = "$input_dir/$species"
species_output_dir = "$output_dir/$species"

# Create output structure
mkdir - p "$species_output_dir" / {
  fastq,
  featurecounts,
  sample_alignments_summary,
  qualimap_reports,
  hisat2_index
}

# Environment setup
module purge
module load HISAT2 / 2.1 .0 - intel - 2018 a
module load SAMtools / 1.9 - foss - 2018 a
module load Java
module load picard / 2.18 .23 - Java - 1.8 .0_171
module load pigz
source / vsc - hard - mounts / leuven - data / 354 / vsc35429 / miniconda3 / etc / profile.d / conda.sh

# Progress logging
PROGRESS_LOG = "$species_output_dir/Pipeline_Progress.log"
echo "$(date): Pipeline started." > "$PROGRESS_LOG"

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# SRA TO FASTQ conversion #
# # # # # # # # # # # # # # # # # # # # # # # # # # #

conda activate featurecounts
mkdir - p "$species_output_dir/tmp"

cd $species_dir

echo "$(date): SRA conversion starting" >> "$PROGRESS_LOG"
# Process only one SRA file using head - n1
find "$species_dir" - type f - name "*.sra" |
  while read - r sra_file;
do
  # Capture the full filename, which includes ".sra"
prefix = $(basename "$sra_file") # e.g., SRR11176917.sra

# Use temporary directory
for conversion
fasterq - dump "$sra_file"--split - 3 - e 8 - O "$species_output_dir/tmp/"

# After conversion, fasterq - dump will produce files:
  # $species_output_dir / tmp / $ {
    prefix
  }
_1.fastq and $ {
  prefix
}
_2.fastq
# Compress these files in parallel.
pigz - p 8 "$species_output_dir/tmp/${prefix}_1.fastq"
"$species_output_dir/tmp/${prefix}_2.fastq"

# Rename the compressed files to remove the ".sra"
portion.
mv "$species_output_dir/tmp/${prefix}_1.fastq.gz"
"$species_output_dir/tmp/${prefix%.sra}_1.fastq.gz"
mv "$species_output_dir/tmp/${prefix}_2.fastq.gz"
"$species_output_dir/tmp/${prefix%.sra}_2.fastq.gz"

# Move the renamed FASTQ files to the final fastq output directory.
mv "$species_output_dir/tmp/" * .fastq.gz "$species_output_dir/fastq/"
done
rm - rf "$species_output_dir/tmp"
echo "$(date): SRA conversion completed successfully" >> "$PROGRESS_LOG"

conda deactivate

# # # # # # # # # # # # # # # # # # # # #
# Trim Galore stage #
# # # # # # # # # # # # # # # # # # # # #
conda activate RNA
echo "$(date): Adapter trimming starting" >> "$PROGRESS_LOG"
TRIMMED_DIR = "$species_output_dir/trimmed"
mkdir - p "$TRIMMED_DIR"

find "$species_output_dir/fastq" - maxdepth 1 - name "*_1.fastq.gz" |
  while read - r FASTQ1;
do
  base = $(basename "$FASTQ1"
    _1.fastq.gz)
FASTQ2 = "${FASTQ1%_1.fastq.gz}_2.fastq.gz"

trim_galore--paired--fastqc--cores 8\
--output_dir "$TRIMMED_DIR"\
"$FASTQ1"
"$FASTQ2"

# Rename trimmed files to a consistent naming convention
mv "$TRIMMED_DIR/${base}_1_val_1.fq.gz"
"$TRIMMED_DIR/${base}.filtered.1.fastq.gz"
mv "$TRIMMED_DIR/${base}_2_val_2.fq.gz"
"$TRIMMED_DIR/${base}.filtered.2.fastq.gz"
done
conda deactivate

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# HISAT2 alignment
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

conda activate featurecounts

mkdir - p "$species_output_dir/sambam"

echo "$(date): HISAT2 indexing starting" >> "$PROGRESS_LOG"
hisat2_index = "$species_output_dir/hisat2_index/$species"
chmod - R 755 "$species_output_dir/hisat2_index"
hisat2 - build - p 16 "$reference_fa"
"$hisat2_index"

echo "$(date): Alignment starting" >> "$PROGRESS_LOG"
find "$TRIMMED_DIR" - name "*.filtered.1.fastq.gz" |
  while read - r fq1;
do
  base = $(basename "$fq1".filtered .1.fastq.gz)
fq2 = "${fq1%1.fastq.gz}2.fastq.gz"

hisat2 - p 8 - x "$hisat2_index"\ -
  1 "$fq1" - 2 "$fq2"\ -
  S "$species_output_dir/sambam/${base}.sam"\
  --summary - file "$species_output_dir/sample_alignments_summary/${base}.txt"
done

#if 1 and 3, RF
else no option
# hisat2 - p 8 - x "$hisat2_index"\
# - 1 "$fq1" - 2 "$fq2"\
# - S "$species_output_dir/sambam/${base}.sam"\
#--rna - strandness RF\
#--summary - file "$species_output_dir/sample_alignments_summary/${base}.txt"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BAM Processing
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
echo "$(date): BAM conversion starting" >> "$PROGRESS_LOG"
find "$species_output_dir/sambam" - name "*.sam" |
  while read - r sam_file;
do
  base = $(basename "$sam_file".sam)
bam_file = "$species_output_dir/sambam/${base}.filtered.sorted.bam"

# Convert, filter, sort, and index
samtools view - @ 8 - bS "$sam_file" | \
  samtools view - @ 8 - f 0x02 - q 20 - b | \
  samtools sort - @ 8 - o "$bam_file"

samtools index - @ 8 "$bam_file"
rm "$sam_file"
done

# # # # # # # # # # # # # # # # # # # # # # # # #
# BAM quality control(QualiMap)
# # # # # # # # # # # # # # # # # # # # # # # # #
#echo "$(date): Running QualiMap QC" >> "$PROGRESS_LOG"
#find "$species_output_dir/sambam" - name "*.filtered.sorted.bam" |
  while read - r bam_file;
do
  # base = $(basename "$bam_file".filtered.sorted.bam)
# qualimap_dir = "$species_output_dir/qualimap_reports/$base"
#
# qualimap rnaseq - bam "$bam_file" - gtf "$reference_gtf" - outdir "$qualimap_dir"\
#--java - mem - size = 100 G - pe
#done

# # # # # # # # # # # # # # # # # # # # # # # # #
# FeatureCounts
# # # # # # # # # # # # # # # # # # # # # # # # #
conda activate featurecounts
echo "$(date): Quantification starting" >> "$PROGRESS_LOG"

# Determine strandness
strandness = 0
if [
  ["$species" == "1"]
] || [
  ["$species" == "3"]
];
then
strandness = 2
fi

# Process all BAM files in one run
featureCounts - T 16 - p - s $strandness - a "$reference_gtf"\ -
  t exon - g gene_id\ -
  o "$species_output_dir/featurecounts/gene_counts.txt"\
"$species_output_dir" / sambam
/*.filtered.sorted.bam

# Create clean matrix
#cut -f1,7- "$species_output_dir/featurecounts/gene_counts.txt" | \
#  grep -v '^#' > "$species_output_dir/featurecounts/gene_counts_matrix.tsv"
conda deactivate

echo "$(date): Pipeline completed successfully" >> "$PROGRESS_LOG"
