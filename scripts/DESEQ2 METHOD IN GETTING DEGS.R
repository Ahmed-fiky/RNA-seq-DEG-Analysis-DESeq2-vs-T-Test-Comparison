if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
library(tidyverse)
library(DESeq2)

counts <- read.csv("C:/Users/alfek/Downloads/Diabetes_infection_count.csv", row.names = 1)
meta <- read.csv("C:/Users/alfek/Downloads/Diabetes_infection_pheno.csv", row.names = 1)

dim(counts)
dim(meta)

meta = meta[, c("CL4"), drop = FALSE]
colnames(counts) == rownames(meta)
rownames(meta) <- gsub("-", ".", rownames(meta), fixed = TRUE)
colnames(meta) <- "sample_group"
meta$sample_group <- gsub("Infection", "Infected", meta$sample_group, fixed = TRUE)

all(colnames(counts) == rownames(meta))

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ sample_group
)
dds

# remove all rows having total counts less than 10
filtered_counts <- rowSums(counts(dds)) >= 10

dds <- dds[filtered_counts,]
dds
# set the factor level were reference is the control
dds$sample_group <- relevel(dds$sample_group, ref = "Healthy")
dds$sample_group

# collapse technical replicates
# run deseq
dds <- DESeq(dds)
res <- results(dds)
res


summary(res)
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)



# MA plot a scatter plot used to visualize the differences between two samples, commonly for gene expression data (RNA-Seq, microarrays)
png("MA_plot.png", width=2000, height=1200, res=300)
plotMA(res)
dev.off()
#blue dots are for deg and as we go to right either upper or lower parts are very good candidates to further work on

# Step 3: Run the analysis
dds <- DESeq(dds)

# Step 4: Get the results
res <- results(dds)
res
# Order results by adjusted p-value
res_ordered_padj <- res[order(res$padj), ]
res_ordered_padj


# Get top 20 DEGs
DEG20_DESEQ <- head(as.data.frame(res_ordered_padj), 20)
DEG20_DESEQ 

# Save to file
getwd()
setwd("C:/Users/alfek/Desktop/Bioinformatics projects/DESeq2 vs Python T-test Top 20 DEG Overlap/rna-seq-deg-analysis/rna-seq-deg-analysis/results")
write.csv(DEG20_DESEQ, "DEG20_DESEQ.csv")

#comparing DEGS between python ttest and R DESEQ2

install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

DEG_20_ttst <- read.csv("DEG_20_ttst.csv")

deseq2_genes <- as.character(rownames(DEG20_DESEQ))
ttst_genes <- as.character(DEG_20_ttst$X)

deseq2_genes
ttst_genes

intersect(deseq2_genes, ttst_genes)
print(intersect(deseq2_genes, ttst_genes))

shared <- intersect(deseq2_genes, ttst_genes)



venn_plot <- ggVennDiagram(list(DESeq2 = deseq2_genes, ttest = ttst_genes),
                           set_color = "cyan") +
  scale_fill_gradient(low = "lightgreen", high = "orange") +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "DESeq2 vs Python T-test: Top 20 DEG Overlap",
       caption = paste("Shared:", paste(shared, collapse = ", "))) +
  theme(plot.title = element_text(color = "cyan"),
        plot.caption = element_text(color = "cyan"))
venn_plot

# NOW SAVE THE RESULTS OF THE VENN DIAGRAM

png("DESeq2_vs_ttest_venn.png", width=2000, height=1200, res=300, bg="black")
venn_plot
dev.off()


# Conclusion:
# DESeq2 is better for RNA-seq specifically because:
# - T-test assumes normal distribution and treats each gene independently
#   it was designed for microarrays, not RNA-seq
# - DESeq2 was built specifically for RNA-seq count data, uses a negative
#   binomial model which fits count data naturally, and applies shrinkage
#   to stabilize fold change estimates for lowly expressed genes

# The 5 shared genes (NCAPH2, RARA-AS1, JUNB, DLG4, KLF13) are the most
# reliable DEGs — both methods agree on them despite their different
# statistical approaches, making them the strongest candidates for follow-up