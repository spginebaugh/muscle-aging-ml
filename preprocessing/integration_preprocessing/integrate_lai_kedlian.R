#' preprocessing of individual datasets must be run prior to this code
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
lai_seurat <- qread("data/processed/lai_preprocessed.qs.qs")
ked_seurat <- qread("data/processed/kedlian_human_relabel.qs")

keep_genes <- intersect(rownames(lai_seurat), rownames(ked_seurat))
lai_seurat <- lai_seurat[rownames(lai_seurat) %in% keep_genes,]
ked_seurat <- ked_seurat[rownames(ked_seurat) %in% keep_genes,]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               merge datasets                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_seurat <- merge(x = lai_seurat, y = ked_seurat)
merge_seurat <- JoinLayers(merge_seurat)

merge_seurat <- merge_seurat %>% NormalizeData()
merge_seurat <- merge_seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
merge_seurat <- merge_seurat %>% FindVariableFeatures(nfeatures = 2000)
merge_seurat <- merge_seurat %>% ScaleData()
merge_seurat <- merge_seurat %>% RunPCA()
merge_seurat <- merge_seurat %>% RunHarmony(group.by.vars = "donor_tech", dims.use = 1:30)
merge_seurat <- merge_seurat %>% RunUMAP(reduction = "harmony", dims = 1:30)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Prelim plotting                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DimPlot(merge_seurat, group.by = c("donor","tech"))
DimPlot(merge_seurat, group.by = c("donor"))
DimPlot(merge_seurat, group.by = c("sample"))
DimPlot(merge_seurat, group.by = c("Phase", "sex","age_bin"))
DimPlot(merge_seurat, group.by = c("manuscript_anno1"), label = TRUE)
DimPlot(merge_seurat, group.by = c("manuscript_anno2"))
DimPlot(merge_seurat, group.by = c("manuscript_anno3"))   
DimPlot(merge_seurat, group.by = c("annotation_level1")) 
DimPlot(merge_seurat, group.by = c("annotation_level1"), split.by = c("age_bin"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(merge_seurat, "data/processed/integrated_lai_kedlian.qs")
