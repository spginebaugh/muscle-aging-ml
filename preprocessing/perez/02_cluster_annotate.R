# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(harmony)
library(stringr)
set.seed(43648)

source("scripts/utils/rna_utils.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/filtered_seurat.qs")

load("data/metadata/cycle.rda") # cell cycle

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Normalize and Cluster                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(seurat) <- "B01" 
seurat <- seurat %>% NormalizeData()
seurat <- seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
seurat <- seurat %>% FindVariableFeatures()
seurat <- seurat %>% ScaleData()
seurat <- seurat %>% RunPCA()

## batch correct
seurat <- seurat %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40, assay.use = "B01")
seurat <- seurat %>% RunUMAP(reduction = "harmony", dims = 1:40)
DimPlot(seurat, group.by = c("sample", "age_group"))

## cluster
seurat <- seurat %>% FindNeighbors(reduction = "harmony", dims = 1:40)
seurat <- seurat %>% FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8))

## select resolution
DimPlot(seurat, group.by = c("B01_snn_res.0.2"), label = TRUE)

## QC
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("nCount_B01", "nFeature_B01", "mitoRatio_B01", "riboRatio_B01")
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Annotate                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## score with Panglao gene lists
panglao_scores <- score_panglao(seurat)

## info plots
DimPlot(seurat, group.by = c("sample", "age_group"))
DimPlot(seurat, group.by = c("B01_snn_res.0.2"), label = TRUE)

## Marker FeaturePlots
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("PDGFRB")) # smooth muscle
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("PAX7")) # muscle stem cell
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("SOX2","SOX10")) # muscle stem cell

## doublets
DimPlot(seurat, group.by = c("scDbl_class", "DF_classification"))
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("cxds_scores", "bcds_scores", "hybrid_scores", "DF_score")
)

## annotation
seurat$annotation_level1 <- plyr::revalue(seurat$B01_snn_res.0.2, 
                                      c('0' = "Fast_SKM",
                                        '1' = "Slow_SKM",
                                        '2' = "FAP",
                                        '3' = "Endothelial",
                                        '4' = "SKM_prog",
                                        '5' = "VSMC",
                                        '6' = "Myeloid",
                                        '7' = "MuSC",
                                        '8' = "NMJ_1",
                                        '9' = "NMJ_2",
                                        '10' = "Lymphoid",
                                        '11' = "Mast"))
seurat$annotation_level1 <- factor(seurat$annotation_level1)
## check annotation
DimPlot(seurat, group.by = "annotation_level1", label = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    Save                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Idents(seurat) <- seurat$annotation_level1
qsave(seurat, "data/processed/annotated_seurat_perez.qs")
