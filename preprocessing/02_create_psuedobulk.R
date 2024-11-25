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

source("utils/mofa_prep_functions.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/lai_preprocessed.qs")
output_dir <- file.path("data/mofa_inputs/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                prep for MOFA                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_annotation <- "donor_tech"
cell_annotation <- "annotation_level1"

## create pseudobulk counts of each celltype for each sample (called "views")
pb_counts <- prep_pseudobulk_counts(seurat, sample_annotation, cell_annotation)
pb_meta <- prep_pseudobulk_metadata(seurat, sample_annotation, cell_annotation)

## remove unwanted celltypes with low cell counts
cell_types_remove <- c("Bcell", "Erythrocyte", "Mast", "Schwann", "Adipocyte")
cell_type <- unique(as.character(pb_meta$cell_type))
cell_type <- cell_type[!(cell_type %in% cell_types_remove)]

pb_counts <- pb_counts[, pb_meta$cell_type %in% cell_type]
pb_meta <- pb_meta[pb_meta$cell_type %in% cell_type, ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            muscle specific DEGs                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' preliminary runs show that RNA contamination from muscle-specific genes
#'   e.g. TPM2 is prevalent in the data, even after RNA contamination removal
#'
#' To prevent this from impacting our MOFAcell run, we will filter out
#'   muscle specific marker genes
#'
#' We have to be careful here though, we don't want to filter out important
#'   nonspecific genes that are upregulated in muscle
#'
#' muscle is only in snRNAseq data, but its contamination is seen in both
#'   cells and nuclei

## remove MuSC and SMC because those are similar to muscle cells at later stages of differentiation
nuclei <- seurat[, seurat$tech == "nuclei" & !(seurat$annotation_level1 %in% c(cell_types_remove))]

Idents(nuclei) <- nuclei$annotation_level1
sc_markers_list_fast <- single_cell_de_consensus(nuclei, "annotation_level1", "FastSKM", c("SlowSKM","MuSC", "SMC"))
sc_markers_list_slow <- single_cell_de_consensus(nuclei, "annotation_level1", "SlowSKM", c("FastSKM","MuSC", "SMC"))

sc_markers_fast <- Reduce(intersect, sc_markers_list_fast)
sc_markers_slow <- Reduce(intersect, sc_markers_list_slow)

pb_markers_list_fast <- pseudobulk_marker_genes(pb_counts, pb_meta, "FastSKM", c("SlowSKM","MuSC","SMC"))
pb_markers_list_slow <- pseudobulk_marker_genes(pb_counts, pb_meta, "SlowSKM", c("FastSKM","MuSC","SMC"))

pb_markers_fast <- Reduce(intersect, pb_markers_list_fast)
pb_markers_slow <- Reduce(intersect, pb_markers_list_slow)

muscle_background_genes <- unique(intersect(sc_markers_fast, pb_markers_fast),
                                  intersect(sc_markers_slow, pb_markers_slow))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Sample Metadata                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_meta <- unique(seurat@meta.data[, c("donor_tech", "age_bin", "age_numeric","tech","sex")])
colnames(sample_meta)[1] <- "sample"
sample_id_column <- "sample"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 save files                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(pb_counts, paste0(output_dir, "pb_counts.csv"), row.names = TRUE)
write.csv(pb_meta, paste0(output_dir, "pb_meta.csv"), row.names = TRUE)
write.csv(muscle_background_genes, paste0(output_dir, "muscle_background_genes.csv"), row.names = FALSE)
write.csv(sample_meta, paste0(output_dir, "sample_metadata.csv"), row.names = FALSE)
