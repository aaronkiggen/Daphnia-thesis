# FST Analysis Pipeline – Chaturvedi et al.

This repository contains a three-step pipeline for read mapping, BAM processing, and SNP calling, used in the analysis of genome-wide differentiation (F<sub>ST</sub>) in *Daphnia magna*. The scripts are designed for use on the VSC HPC cluster (Genius) and focus on reproducibility and parallel processing of multiple samples and chromosomes.

---

## 📁 Directory Structure
- `README.md` — Pipeline documentation
- `scripts/`
  - `1_mapping_filtering.sh` — Mapping reads, filtering, sorting, deduplication
  - `2_snp_calling.sh` — Per-chromosome SNP calling
  - `parseVCF.py` — Custom VCF filtering script
- `bams/` — BAM files (input/output)
- `genome/` — Reference genome FASTA
- `final/` — Final filtered SNP calls
- `logs/` — SLURM output logs

---

## 🚀 Pipeline Overview

### Step 1: Mapping & BAM Processing

**Script:** `1_mapping_filtering.sh`  
**Function:**  
- Maps paired-end reads to the reference genome using `BWA-MEM`  
- Filters for properly paired reads with mapping quality ≥ 20  
- Sorts BAM files and removes PCR duplicates using Picard  
- Indexes the final BAM file  

**Run with SLURM Array:**
```bash
sbatch --array=1-24 scripts/1_mapping_filtering.sh


### ✅ Step 2: Per-Chromosome SNP Calling

**Script:** `2_snp_calling.sh`  
**Function:**  
- Runs `bcftools mpileup` and `bcftools call` per chromosome across all samples  
- Normalizes VCFs to biallelic SNPs  
- Optionally filters out indels using `vcftools`  
- Applies additional filtering with `parseVCF.py`  
- Renames sample IDs by stripping BAM suffixes  

**Run with SLURM Array:**
```bash
sbatch --array=1-13 scripts/2_snp_calling.sh


### ✅ Step 3: Custom VCF Filtering

Script: parseVCF.py
Function:

Filters variants based on the following criteria:

- Genotype Quality (GQ ≥ 30)

- Read Depth (DP ≥ 10)

- Genotype Type: Retains only heterozygous and homozygous alternate genotypes

Skips indels to retain SNPs only

Outputs compressed, filtered VCF files


### 📦 Output Files
| Output Type              | Description                          | Location                             |
|--------------------------|------------------------------------|------------------------------------|
| BAM Files                | Final deduplicated BAMs             | `bams/*.filtered.sorted.nd.bam`    |
| Raw VCFs per chromosome  | Initial VCFs from SNP calling       | `final/chaturvedi_chr_*.vcf.gz`    |
| Filtered SNP calls       | High-quality SNPs after filtering   | `final/chaturvedi_magna_chr_*.H.calls.gz` |

