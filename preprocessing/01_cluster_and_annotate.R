# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(harmony)
set.seed(43648)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts <- Read10X("data/processed/converted_files/mtx_files/", gene.column = 1)

metadata <- read.csv("data/processed/converted_files/metadata.csv")
colnames(metadata)[1] <- "barcode"
rownames(metadata) <- metadata$barcode

load("data/util_data/cycle.rda") 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              organize metadata                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rename columns and data types so metadata is consistent across datasets
## add additional information for easier intergration with other datasets
metadata <- metadata %>% dplyr::rename(
  "manuscript_anno1" = "Annotation",
  "sex" = "Sex",
  "age_bin" = "age_pop",
  "age_numeric" = "age",
  "donor" = "sample",
  "sample" = "orig.ident"
)

metadata$tech <- plyr::revalue(metadata$tech,
                               c("scRNA" = "cells",
                                 "snRNA" = "nuclei"))

metadata$sex <- plyr::revalue(metadata$sex,
                              c("Male" = "M",
                                "Female" = "F"))
metadata$donor_tech <- paste0(metadata$donor, "_", metadata$tech)
metadata$study <- "lai_muscle_multiomics"

## create seurat object
seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     QC                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## add qc stats
seurat$mitoRatio <- PercentageFeatureSet(object = seurat, pattern = "^MT-") / 100
seurat$riboRatio <- PercentageFeatureSet(object = seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)

## filtering was already performed by authors of original manuscript
## no additional filtering necessary
metadata <- seurat@meta.data
metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")

metadata %>%
  ggplot(aes(color = sample, x = nFeature_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 400)

metadata %>%
  ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2) 

metadata %>%
  ggplot(aes(color = sample, x = riboRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.04)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  cluster                                 ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- seurat %>% NormalizeData()
seurat <- seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
seurat <- seurat %>% FindVariableFeatures(nfeatures = 3000)
seurat <- seurat %>% ScaleData(vars.to.regress = c("mitoRatio","nCount_RNA"))
seurat <- seurat %>% RunPCA()

seurat <- seurat %>% RunHarmony(group.by.vars = "donor_tech", dims.use = 1:30)
seurat <- seurat %>% RunUMAP(reduction = "harmony", dims = 1:30)
seurat <- seurat %>% RunUMAP(reduction = "pca", dims = 1:30, reduction.name = "uncorrected_umap")

DimPlot(seurat, group.by = c("donor","tech"))
DimPlot(seurat, group.by = c("donor"))
DimPlot(seurat, group.by = c("Phase", "sex","age_bin"))
DimPlot(seurat, group.by = c("manuscript_anno1"), label = TRUE)


# change names to be consistent across datasets
seurat$annotation_level1 <- plyr::revalue(seurat$manuscript_anno1, 
                                    c("Adipocyte" = "Adipocyte",
                                      "Endo" = "Endothelial",
                                      "Erythrocyte" = "Erythrocyte",
                                      "FAP" = "FAP",
                                      "Lymphocyte" = "Lymphoid",
                                      "Mast cell" = "Mast",
                                      "MuSC" = "MuSC",
                                      "Myeloid cell" = "Myeloid",
                                      "Pericyte" = "Pericyte",
                                      "Schwann cell" = "Schwann",
                                      "SMC" = "SMC",
                                      "Specialized MF" = "FastSKM",
                                      "Tenocyte" = "Tenocyte",
                                      "Type I" = "SlowSKM",
                                      "Type II" = "FastSKM"))

DimPlot(seurat, group.by = "annotation_level1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, "data/processed/lai_preprocessed.qs")



