# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")
# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# read in sample info
colData <- read.csv('sample_info.csv')

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))

# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ dexamethasone)
dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds
# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res

# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)

#Volcano plot
library(ggplot2)

# Add log fold-change and significance
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$threshold <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")


#Enrichment Analysis
library(dplyr)

sig_genes <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# Extract gene names
deg_list <- rownames(sig_genes)
write.csv(deg_list, "DEGs.csv", row.names = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)  # Use appropriate organism database

# Convert gene names to Entrez IDs
gene_entrez <- mapIds(org.Hs.eg.db, deg_list, keytype = "SYMBOL", column = "ENTREZID")

# GO Enrichment Analysis
ego <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# KEGG Pathway Enrichment
kegg <- enrichKEGG(gene = gene_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Visualize results
dotplot(ego, showCategory = 10)
dotplot(kegg, showCategory = 10)
barplot(kegg, showCategory = 10, title = "KEGG Pathway Enrichment")
barplot(GO, showCategory = 10, title = "GO Pathway Enrichment")

#Enriched Pathways
library(enrichplot)

emapplot(ego)   # For GO analysis
emapplot(kegg)  # For KEGG pathways

#Save results
write.csv(as.data.frame(ego), "GO_enrichment.csv")
write.csv(as.data.frame(kegg), "KEGG_enrichment.csv")


