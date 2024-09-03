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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/filtered_seurat.qs")

load("data/processed/cycle.rda") # cell cycle

doublet_meta <- read.csv("data/processed/perez_doublet_metadata.csv") # doublet metadata from 00_doublet detection
doublet_meta <- doublet_meta[, 5:ncol(doublet_meta)] # remove redundant info

## add doublet metadata to seurat
metadata <- seurat@meta.data
metadata <- left_join(metadata, doublet_meta, by = "barcodes")
rownames(metadata) <- metadata$barcodes
seurat@meta.data <- metadata


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Normalize and Cluster                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(seurat) <- "RNA_ranger" # use data without ambient RNA removal
seurat <- seurat %>% NormalizeData()
seurat <- seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
seurat <- seurat %>% FindVariableFeatures()
seurat <- seurat %>% ScaleData()
seurat <- seurat %>% RunPCA()
seurat <- seurat %>% RunUMAP(reduction = "pca", dims = 1:40)

## view umap to decide if batch correction is necessary
DimPlot(seurat, group.by = c("sample", "age_group"))

## batch correct
seurat <- seurat %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40, assay.use = "RNA_ranger")
seurat <- seurat %>% RunUMAP(reduction = "harmony", dims = 1:40)
DimPlot(seurat, group.by = c("sample", "age_group"))

## cluster
seurat <- seurat %>% FindNeighbors(reduction = "harmony", dims = 1:40)
seurat <- seurat %>% FindClusters(resolution = c(0.2, 0.4))

## select resolution
DimPlot(seurat, group.by = c("RNA_ranger_snn_res.0.2", "RNA_ranger_snn_res.0.4"))

## select snn_res.0.4 as resolution
Idents(seurat) <- seurat$RNA_ranger_snn_res.0.4

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Annotate                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## score with Panglao gene lists
seurat_ann <- seurat

## TODO: clean this up
panglao <- readr::read_tsv("data/processed/PanglaoDB_markers_27_Mar_2020.tsv")
panglao <- panglao[panglao$species %in% c("Hs", "Mm Hs"), ]
panglao <- panglao[panglao$`official gene symbol` %in% rownames(seurat_ann), ]
panglao <- split(panglao$`official gene symbol`, panglao$`cell type`)

seurat_ann <- AddModuleScore(seurat_ann, panglao)
colnames(seurat_ann@meta.data)[(ncol(seurat_ann@meta.data) - length(panglao) + 1):ncol(seurat_ann@meta.data)] <- names(panglao)

pang_ann <- seurat_ann@meta.data[, (ncol(seurat_ann@meta.data) - length(panglao) + 1):ncol(seurat_ann@meta.data)]
pang_ann <- scale(pang_ann)
agg_pang <- aggregate(pang_ann, list(seurat_ann@active.ident), mean)
rownames(agg_pang) <- agg_pang[, 1]
agg_pang <- agg_pang[, -1]
agg_pang <- t(agg_pang)

## info plots
DimPlot(seurat, group.by = c("sample", "age_group"))
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("nCount_RNA", "nFeature_RNA", "mitoRatio", "riboRatio")
)
DimPlot(seurat, group.by = c("RNA_ranger_snn_res.0.4"), label = TRUE)

## Marker FeaturePlots
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("PDGFRB")) # smooth muscle
FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("PAX7")) # muscle stem cell

## doublets
DimPlot(seurat, group.by = c("scDbl_class", "DF_classification"))
FeaturePlot(
  seurat,
  min.cutoff = "q1", max.cutoff = "q99",
  features = c("cxds_scores", "bcds_scores", "hybrid_scores", "DF_score")
)

## annotation
seurat$annotation_level1 <- plyr::revalue(
  seurat$RNA_ranger_snn_res.0.4,
  c(
    "0" = "fastSKM",
    "1" = "slowSKM",
    "2" = "fastSKM",
    "3" = "slowSKM",
    "4" = "FAP",
    "5" = "endothelial",
    "6" = "slowSKM",
    "7" = "NMJpost",
    "8" = "Myeloid",
    "9" = "MuSC",
    "10" = "vsmc",
    "11" = "FAP",
    "12" = "slowSKM",
    "13" = "NMJpre",
    "14" = "Tcell",
    "15" = "endothelial",
    "16" = "mast",
    "17" = "vsmc",
    "18" = "vsmc",
    "19" = "fastSKM"
  )
)

## check annotation
DimPlot(seurat, group.by = "annotation_level1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    Save                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Idents(seurat) <- seurat$annotation_level1
qsave(seurat, "data/processed/annotated_seurat_perez.qs")
