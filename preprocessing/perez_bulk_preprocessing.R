# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(MOFA2)
library(MOFAcellulaR)
library(qs)
library(dplyr)
library(tidyverse)
library(magrittr)
library(edgeR)
library(DESeq2)
library(nichenetr)

library(ggprism)
library(ggpmisc)

set.seed(43648)

source("utils/rna_utils.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bulk_perez_cts <- read.csv("data/raw/bulk_RNA/perez_bulk/bulk_counts.csv", row.names = "Symbol")

## organize metadata
bulk_perez_meta <- GEOquery::getGEO(filename = "data/raw/bulk_RNA/perez_bulk/GSE167186-GPL20301_series_matrix.txt.gz")
bulk_perez_meta <- bulk_perez_meta@phenoData@data
bulk_perez_meta <- bulk_perez_meta[, c(1, 52:59)]
colnames(bulk_perez_meta) <- c("sample", "walkTest", "ageNum", "biodex", "gripStrength", "group", "legPress", "sppb", "upAndGo")

rownames(bulk_perez_meta) <- bulk_perez_meta$sample
bulk_perez_meta$walkTest %<>% as.numeric()
bulk_perez_meta$ageNum %<>% as.numeric()
bulk_perez_meta$biodex %<>% as.numeric()
bulk_perez_meta$gripStrength %<>% as.numeric()
bulk_perez_meta$legPress %<>% as.numeric()
bulk_perez_meta$sppb %<>% as.numeric()
bulk_perez_meta$upAndGo %<>% as.numeric()

bulk_perez_meta$group[bulk_perez_meta$group == "UNCLASSIFIED"] <- "Old Healthy"
bulk_perez_meta$group %<>% as.factor()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Filter and Normalize                          ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bulk_perez_cts <- bulk_perez_cts[rowSums(bulk_perez_cts > 10) > 3, ]

dds_perez <- DESeqDataSetFromMatrix(countData = bulk_perez_cts, colData = bulk_perez_meta, design = ~group)
vsd_perez_dst <- vst(dds_perez)
vsd_in_perez <- assay(vsd_perez_dst)
# vsd_in_perez <- vsd_in_perez[rownames(vsd_in_perez) %in% rownames(factor_loadings_mat), ]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    save                                  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(vsd_in_perez, "data/processed/perez_bulk_normalized.csv", row.names = TRUE)
write.csv(bulk_perez_meta, "data/processed/perez_bulk_metadata.csv", row.names = TRUE)





