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
library(stringr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Import Cellbender Outputs                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed/perez_doublets_detected.qs")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Add in Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fix mistake in metadata colnames for consistency
colnames(seurat@meta.data) <- str_replace_all(colnames(seurat@meta.data),"b0","B0")

## get patient info from metadata csv
runinfo <- read.csv("data/metadata/Perez_SraRunInfo.csv")
runinfo <- runinfo[, c("GEO_Accession..exp.", "Group", "Age")]
colnames(runinfo) <- c("sample", "age_group", "age_numeric")

## organize seurat object metadata
metadata <- seurat@meta.data
metadata <- left_join(metadata, runinfo, by = "sample")
rownames(metadata) <- metadata$barcode
seurat@meta.data <- metadata

## add in cell quality information
seurat$mitoRatio_ranger <- PercentageFeatureSet(object = seurat, pattern = "^MT-", assay = "RNA") / 100
seurat$riboRatio_ranger <- PercentageFeatureSet(object = seurat, 
                                                pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", 
                                                assay = "RNA") / 100
seurat$log10GenesPerUMI_ranger <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)

seurat$mitoRatio_B01 <- PercentageFeatureSet(object = seurat, pattern = "^MT-", assay = "B01") / 100
seurat$riboRatio_B01 <- PercentageFeatureSet(object = seurat, 
                                             pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                                             assay = "B01") / 100
seurat$log10GenesPerUMI_B01 <- log10(seurat$nFeature_B01) / log10(seurat$nCount_B01)

seurat$mitoRatio_B05 <- PercentageFeatureSet(object = seurat, pattern = "^MT-", assay = "B05") / 100
seurat$riboRatio_B05 <- PercentageFeatureSet(object = seurat, 
                                             pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                                             assay = "B05") / 100
seurat$log10GenesPerUMI_B05 <- log10(seurat$nFeature_B05) / log10(seurat$nCount_B05)
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
  ggplot(aes(color = sample, x = nFeature_B05, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)

metadata %>%
  ggplot(aes(color = sample, x = nCount_B05, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio_ranger, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.02)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio_B01, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.02)

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio_B05, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.02)

metadata %>%
  ggplot(aes(color = sample, x = riboRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.04)

metadata %>%
  ggplot(aes(x = sample, log10GenesPerUMI_B05, fill = sample)) +
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
## Lots of noise in this dataset. Want to have strict QC filtering here
seurat <- subset(seurat, nFeature_B05 > 300 & 
                   nCount_B05 > 500 & 
                   mitoRatio_ranger < 0.02 & 
                   mitoRatio_B01 < 0.02 &
                   mitoRatio_B05 < 0.02 &
                   log10GenesPerUMI_B05 > 0.8)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, "data/processed/filtered_seurat.qs")
