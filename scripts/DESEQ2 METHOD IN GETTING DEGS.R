if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(tidyverse)
library(DESeq2)

counts <- read.csv("C:/Users/alfek/Downloads/Diabetes_infection_count.csv", row.names = 1)
pheno <- read.csv("C:/Users/alfek/Downloads/Diabetes_infection_pheno.csv", row.names = 1)

dim(counts)
dim(pheno)

pheno = pheno[, c("CL4"), drop = FALSE]
rownames(pheno) <- gsub("-", ".", rownames(pheno), fixed = TRUE)
colnames(pheno) <- "sample_group"
pheno$sample_group <- gsub("Infection", "Infected", pheno$sample_group, fixed = TRUE)
pheno
all(colnames(counts) == rownames(pheno))


# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = pheno,
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

#contrasts
resultsNames(dds)

# e.g: treated_4hrs, treated_8hrs, untreated

# for e.g, results(dds, contrast = c("insulin", "treated_4hrs", "Healthy"))

# MA plot a scatter plot used to visualize the differences between two samples, commonly for gene expression data (RNA-Seq, microarrays)
plotMA(res)

#blue dots are for deg and as we go to right either upper or lower parts are very good candidates to further work on

j# Step 3: Run the analysis
dds <- DESeq(dds)

# Step 4: Get the results
res <- results(dds)

# Order results by adjusted p-value
res_ordered_pval <- res[order(res$padj), ]

# Get top 20 DEGs
top20_R <- head(as.data.frame(res_ordered), 20)
top20_R

# Save to file
write.csv(top20, "top20_DEGs_R_pval.csv")
getwd()
setwd("C:/Users/alfek/Downloads")

#comparing DEGS between python ttest and R DESEQ2

install.packages("ggVennDiagram")
library(ggVennDiagram)

python_DEG20 <- read.csv("DEG_20.csv")

library(ggVennDiagram)

deseq2_genes <- as.character(rownames(top20_R))
python_genes <- as.character(python_DEG20$X)

intersect(deseq2_genes, python_genes)
print(intersect(deseq2_genes, python_genes))

shared <- intersect(deseq2_genes, python_genes)

ggVennDiagram(list(DESeq2 = deseq2_genes, Python_ttest = python_genes)) +
  scale_fill_gradient(low = "lightgreen", high = "orange") +
  scale_color_manual(values = c("blue", "red")) +
  labs(caption = paste("Shared:", paste(shared, collapse = ", ")))
