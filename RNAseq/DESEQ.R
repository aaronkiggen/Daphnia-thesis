#### LIBRARIES ####
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggrepel)


########## change depending on what species ############ 
dir_path <- "C:/Users/u0173289/OneDrive - KU Leuven/PROJECTS/thesis to paper/featurecounts/magna_8"
species_name <- "magna_8"
clone_of_interest <- "M8" #FOR MITSUKURI AND SINENSIS ITS DAY

##### PCA ############

counts_file <- file.path(dir_path, "gene_counts.txt")
metadata_file <- file.path(dir_path, "metadata.csv")
# ----- Process the FeatureCounts output -----
# Read featureCounts output; skip commented lines (starting with "#")
counts_raw <- read.table(counts_file, header = TRUE, sep = "\t", comment.char = "#")

# Process column names: extract sample IDs from file paths.
# This example assumes the file names end with ".filtered.sorted.bam" which will be removed.
#colnames(counts_raw) <- gsub(".*\\.sambam\\.(SRR[0-9]+)\\.filtered\\.sorted\\.bam", "\\1", colnames(counts_raw)) #3 4
#colnames(counts_raw) <- gsub(".*\\.sambam\\.(ERR[0-9]+)\\.filtered\\.sorted\\.bam", "\\1", colnames(counts_raw)) #5
#colnames(counts_raw) <- gsub(".*\\.sambam\\.(SRR[0-9]+).*\\.filtered\\.sorted\\.bam", "\\1", colnames(counts_raw)) #1 mitsukuri
colnames(counts_raw) <- gsub(".*\\.sambam\\.([^.]+)\\.filtered\\.sorted\\.bam", "\\1", colnames(counts_raw)) #magna

colnames(counts_raw)

background_genes <- counts_raw$Geneid
write.csv(background_genes, file = paste0(species_name, "background_genes"))


# Remove the first six annotation columns
counts <- counts_raw[, -c(1:6)]
# Set rownames from the Geneid column in counts_raw
rownames(counts) <- counts_raw$Geneid


str(counts)
# ----- Read metadata -----
# Assuming metadata is in CSV format with a header and a column "sample"
metadata <- read.csv(metadata_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
str(metadata)

# Set row names of metadata to sample names
rownames(metadata) <- metadata$sample



str(metadata)
str(counts)
samples_in_counts <- colnames(counts)
metadata <- metadata[metadata$sample %in% samples_in_counts, ]
# Reorder metadata to match the order of columns in counts
metadata <- metadata[colnames(counts), ]


#clone of interest
metadata <- metadata[metadata$clone == clone_of_interest, ] #DAY

counts <- counts[, rownames(metadata)]

# ----- PCA by Condition -----
dds_condition <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                                colData = metadata,
                                                design = ~ condition)
vsd_condition <- DESeq2::vst(dds_condition, blind = FALSE)
pca_data_condition <- DESeq2::plotPCA(vsd_condition, intgroup = "condition", returnData = TRUE)
percentVar_condition <- round(100 * attr(pca_data_condition, "percentVar"))

p_condition <- ggplot2::ggplot(pca_data_condition, 
                               ggplot2::aes(PC1, PC2, color = condition, label = name)) +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +  # Adds sample names to the dots
  ggplot2::ggtitle(paste("PCA for", species_name, "- Condition")) +
  ggplot2::xlab(paste0("PC1: ", percentVar_condition[1], "% variance")) +
  ggplot2::ylab(paste0("PC2: ", percentVar_condition[2], "% variance")) +
  ggplot2::theme_minimal()
p_condition
# Save objects with species name in their variable names
assign(paste0("dds_condition_", species_name), dds_condition)
assign(paste0("vsd_condition_", species_name), vsd_condition)
assign(paste0("pca_data_condition_", species_name), pca_data_condition)
assign(paste0("percentVar_condition_", species_name), percentVar_condition)
assign(paste0("p_condition_", species_name), p_condition)

# ----- PCA by Clone -----
#dds_clone <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
#                                            colData = metadata,
#                                            design = ~ clone)
#vsd_clone <- DESeq2::vst(dds_clone, blind = FALSE)
#pca_data_clone <- DESeq2::plotPCA(vsd_clone, intgroup = "clone", returnData = TRUE)
#percentVar_clone <- round(100 * attr(pca_data_clone, "percentVar"))

#p_clone <- ggplot2::ggplot(pca_data_clone, ggplot2::aes(PC1, PC2, color = clone)) +
#  ggplot2::geom_point(size = 3) +
#  ggplot2::ggtitle(paste("PCA for", species_name, "- Clone")) +
#  ggplot2::xlab(paste0("PC1: ", percentVar_clone[1], "% variance")) +
#  ggplot2::ylab(paste0("PC2: ", percentVar_clone[2], "% variance")) +
#  ggplot2::theme_minimal()

#assign(paste0("dds_clone_", species_name), dds_clone)
#assign(paste0("vsd_clone_", species_name), vsd_clone)
#assign(paste0("pca_data_clone_", species_name), pca_data_clone)
#assign(paste0("percentVar_clone_", species_name), percentVar_clone)
#assign(paste0("p_clone_", species_name), p_clone)

# ----- PCA by Clone:Condition Interaction -----
# Create an interaction factor combining clone and condition
#metadata$clone_condition <- factor(paste(metadata$clone, metadata$condition, sep = "_"))

#dds_interaction <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
#                                                  colData = metadata,
#                                                  design = ~ clone_condition)
#vsd_interaction <- DESeq2::vst(dds_interaction, blind = FALSE)
#pca_data_interaction <- DESeq2::plotPCA(vsd_interaction, intgroup = "clone_condition", returnData = TRUE)
#percentVar_interaction <- round(100 * attr(pca_data_interaction, "percentVar"))

#p_interaction <- ggplot2::ggplot(pca_data_interaction, ggplot2::aes(PC1, PC2, color = clone_condition)) +
#  ggplot2::geom_point(size = 3) +
#  ggplot2::ggtitle(paste("PCA for", species_name, "- Clone:Condition Interaction")) +
#  ggplot2::xlab(paste0("PC1: ", percentVar_interaction[1], "% variance")) +
#  ggplot2::ylab(paste0("PC2: ", percentVar_interaction[2], "% variance")) +
#  ggplot2::theme_minimal()

#assign(paste0("dds_interaction_", species_name), dds_interaction)
#assign(paste0("vsd_interaction_", species_name), vsd_interaction)
#assign(paste0("pca_data_interaction_", species_name), pca_data_interaction)
#assign(paste0("percentVar_interaction_", species_name), percentVar_interaction)
#assign(paste0("p_interaction_", species_name), p_interaction)

# Finally, return a list with all the plots using species-specific names
#result <- list(condition = p_condition, clone = p_clone, interaction = p_interaction)
#assign(paste0("pca_plots_", species_name), result)



######## DEGs ###########
padj_cutoff = 0.05
lfc_threshold = 1

#dd_interaction <- DESeq(get(paste0("dds_interaction_", species_name)))
#assign(paste0("dd_interaction_", species_name), dd_interaction)

#dd_clone <- DESeq(get(paste0("dds_clone_", species_name)))
#assign(paste0("dd_clone_", species_name), dd_clone)





##### pre-filtering #####
dim(dds_condition)
dim(dds_condition[rowSums(counts(dds_condition)) > 10, ])
dd_condition <- dds_condition[rowSums(counts(dds_condition)) > 10, ]





#res_interaction <- DESeq2::results(get(paste0("dd_interaction_", species_name)), 
#                       alpha = padj_cutoff, 
#                       lfcThreshold = lfc_threshold)
#assign(paste0("res_interaction_", species_name), res_interaction)

#res_clone <- DESeq2::results(get(paste0("dd_clone_", species_name)), 
#                       alpha = padj_cutoff, 
#                       lfcThreshold = lfc_threshold)
#assign(paste0("res_clone_", species_name), res_clone)


dds <- DESeq(dd_condition)

Names <- resultsNames(dds)
Names

res <- results(dds,name="condition_fish_vs_control")


summary(res)
res

resOrdered  <- res[order(res$padj),]
head(resOrdered)
mcols(res)$description
plotMA(res, ylim=c(-2,2))


# Filter significant genes based on adjusted p-values and/or log2 fold change
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# 3. Count the Number of DEGs
num_degs <- nrow(sig_genes)
print(num_degs)

sig_genes <- sig_genes[order(sig_genes$log2FoldChange),]
sig_genes <- sig_genes[order(sig_genes$padj),]
sig_genes

resSig <- subset(resOrdered, padj < 0.1)

genes_m8 <- data.frame(res)



#degs_interaction <- subset(get(paste0("res_interaction_", species_name)), padj < padj_cutoff & abs(log2FoldChange) > lfc_threshold)
#assign(paste0("degs_interaction_", species_name), degs_interaction)
#write.csv(get(paste0("degs_interaction_", species_name)), file = paste0(species_name, "_degs_interaction.csv"))

#degs_clone <- subset(get(paste0("res_clone_", species_name)), padj < padj_cutoff & abs(log2FoldChange) > lfc_threshold)
#assign(paste0("degs_clone_", species_name), degs_clone)
#write.csv(get(paste0("degs_clone_", species_name)), file = paste0(species_name, "_degs_clone.csv"))

assign(paste0("res_condition_", species_name), sig_genes)
write.csv(get(paste0("res_condition_", species_name)), file = paste0(species_name, "_res_condition.csv"))

species_name

