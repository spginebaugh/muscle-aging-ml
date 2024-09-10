# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                                                            ~~
#                   Clustering for Perez snRNAseq Data                      ----
#                                                                            ~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
DefaultAssay(seurat) <- "B01" # use data without ambient RNA removal
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            remove contamination                          ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TODO: compare cellbender outputs, soupx, and other method for decontamination
#' Even with strong QC and ambient RNA removal, there are still many non-muscle
#' cells expressing muscle genes
#' We will remove these by identifying top muscle cell markers, then removing
#' non-muscle cells that express these markers

seurat$ismuscle <- ifelse(seurat$B01_snn_res.0.2 %in% c(1,0), "muscle","other")
Idents(seurat) <- seurat$ismuscle

muscle_markers_1 <- FindMarkers(seurat[, seurat$B01_snn_res.0.2 != 0],
                                ident.1 = "muscle", 
                                ident.2 = "other", 
                                only.pos = TRUE)
muscle_markers_0 <- FindMarkers(seurat[, seurat$B01_snn_res.0.2 != 1],
                                ident.1 = "muscle", 
                                ident.2 = "other", 
                                only.pos = TRUE)

# muscle_marker_genes <- rownames(muscle_markers)[muscle_markers$p_val_adj == 0 & 
#                                                   muscle_markers$avg_log2FC > 2 & 
#                                                   muscle_markers$pct.1 > 0.25]
# 
# seurat <- AddModuleScore(seurat, features = list(muscle_marker_genes = muscle_marker_genes))
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("MYH1","MYH2","MYL3","MYH7B", "MYH7","MYL2"),
)
FeaturePlot(
  seurat, max.cutoff = "q10",
  features = c("MYH1","MYH2","MYL3","MYH7B", "MYH7","MYL2"),
  slot = "counts"
)

muscle_counts <- GetAssayData(seurat, assay = "B01", layer = "counts")[c("MYH1","MYH2","MYL3","MYH7B", "MYH7","MYL2"),] %>%
  data.frame()
seurat$contam <- ifelse(colSums(muscle_counts) > 5, "muscle", "other") %>% as.character()

DimPlot(seurat, group.by = "contam")

keep_cells <- colnames(seurat)[(seurat$contam == "other") | (seurat$B01_snn_res.0.2 %in% c(0,1,9))]

seurat <- seurat[,colnames(seurat) %in% keep_cells]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  recluster                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(seurat) <- "B01" # use data without ambient RNA removal
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
Idents(seurat) <- seurat$B01_snn_res.0.2
## QC
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("nCount_B01", "nFeature_B01", "mitoRatio_B01", "riboRatio_B01")
)

FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("MYH1","MYH2","MYL3","MYH7B", "MYH7","MYL2"),
)

FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("MYPN","MYOT","MYLPF","MYOM1","MYL1"),
)

## remove on last contam cluster
seurat <- seurat[,!(seurat$B01_snn_res.0.2 %in% 12)]
Idents(seurat) <- seurat$B01_snn_res.0.2
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

## check annotation
DimPlot(seurat, group.by = "annotation_level1", label = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    Save                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Idents(seurat) <- seurat$annotation_level1
qsave(seurat, "data/processed/annotated_seurat_perez_v2.qs")
