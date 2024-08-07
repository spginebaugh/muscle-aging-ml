# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                                                            ~~
#          Standard Filtering for Perez snRNAseq Data                      ----
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Functions                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
importSeurat <- function(cellranger_folder_path, file_names_vec, file_h5_path) {
  seurat_list <- list()
  for (file_index in 1:length(file_names_vec)) {
    seurat_data <- Read10X_h5(filename = paste0(
      cellranger_folder_path,
      file_names_vec[file_index],
      file_h5_path
    ))
    seurat_list[file_index] <- CreateSeuratObject(
      counts = seurat_data,
      project = file_names_vec[file_index]
    )
  }
  names(seurat_list) <- file_names_vec

  seurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.id = file_names_vec
  )
  seurat <- JoinLayers(seurat) # merge samples in Seurat V5

  colnames(seurat@meta.data)[1] <- "sample"

  return(seurat)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellbender Outputs                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_names_vec <- list.files("data/raw/cellranger_outputs/")


## loop over cellbender outputs for each sample and place in seurat object
seurat <- importSeurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/cellbender/cellbender_out_filtered.h5"
)

## import uncorrected data from cellranger as alternative and add to seurat object
## this will allow us to use cellbender's empty droplet selection with cellranger outputs
seurat_ranger <- importSeurat(
  cellranger_folder_path = "data/raw/cellranger_outputs/",
  file_names_vec = file_names_vec,
  file_h5_path = "/outs/raw_feature_bc_matrix.h5"
)


## cellbender has better emptydrop removal than cellranger
## subset to only cells, as determined by cellbender
seurat_ranger <- seurat_ranger[, colnames(seurat_ranger) %in% colnames(seurat)]

## normalize, add to original seurat object as separate assay
## then remove unneccesary seurat_ranger object
seurat_ranger %<>% NormalizeData()
seurat[["RNA_ranger"]] <- seurat_ranger[["RNA"]]

rm(seurat_ranger)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Add in Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## get patient info from metadata csv
runinfo <- read.csv("fastqs/Perez_SraRunInfo.csv")
runinfo <- runinfo[, c("GEO_Accession..exp.", "Group", "Age")]
colnames(runinfo) <- c("sample", "age_group", "age_numeric")

## organoize seurat object metadata
seurat$barcodes <- colnames(seurat)
metadata <- seurat@meta.data
metadata <- left_join(metadata, runinfo, by = "sample")
rownames(metadata) <- metadata$barcodes
seurat@meta.data <- metadata

## add in cell quality information
seurat$mitoRatio <- PercentageFeatureSet(object = seurat, pattern = "^MT-") / 100
seurat$riboRatio <- PercentageFeatureSet(object = seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     QC                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  geom_vline(xintercept = 300)

metadata %>%
  ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.05)

metadata %>%
  ggplot(aes(color = sample, x = riboRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.04)

metadata %>%
  ggplot(aes(x = sample, log10GenesPerUMI, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = alpha("white", 0.7)) +
  theme_classic() +
  geom_hline(yintercept = c(0.80)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  NoLegend()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Filtering                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- subset(seurat, nFeature_RNA > 300 & nCount_RNA > 500 & mitoRatio < 0.05 & log10GenesPerUMI > 0.8)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, "data/processed/filtered_seurat.qs")
