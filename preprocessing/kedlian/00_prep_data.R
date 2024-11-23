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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts <- Read10X("data/manuscript_processed_data/kedlian_counts/", gene.column = 1)
colnames(counts) <- paste0("g", colnames(counts)) # so it doesn't start with number

metadata <- read.csv("data/manuscript_processed_data/metadata.csv", row.names = "X")
metadata$barcode <- paste0("g",rownames(metadata))
rownames(metadata) <- metadata$barcode

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              organize metadata                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rename columns and data types so metadata is consistent across datasets
metadata <- metadata %>% dplyr::rename(
  "manuscript_anno1" = "annotation_level0",
  "manuscript_anno2" = "annotation_level1",
  "manuscript_anno3" = "annotation_level2",
  "sex" = "Sex",
  "tech" = "batch",
  "age_bin" = "Age_bin",
  "age_numeric" = "Age_group",
  "donor" = "DonorID",
  "sample" = "SampleID"
)

metadata$age_numeric <- word(metadata$age_numeric,1,1,"-") %>% as.numeric()
metadata$sample <- paste0("g", metadata$sample)

vardata <- read.csv("data/manuscript_processed_data/vardata.csv")

load("data/metadata/cycle.rda") # cell cycle

seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)

seurat <- seurat[,seurat$is_doublet == "False"]

seurat$donor_tech <- paste0(seurat$donor, "_", seurat$tech)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     QC                                   ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## filtering was already performed by authors of original manuscript
seurat$mitoRatio <- PercentageFeatureSet(object = seurat, pattern = "^MT-") / 100
seurat$riboRatio <- PercentageFeatureSet(object = seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)

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
  geom_vline(xintercept = 400) # chose a value slightly higher than original manuscript

metadata %>%
  ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000) # same as manuscript

metadata %>%
  ggplot(aes(color = sample, x = mitoRatio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2) # same as manuscript

## didnt end up using these for filtering
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

DimPlot(seurat, group.by = c("donor","tech"))
DimPlot(seurat, group.by = c("donor"))
DimPlot(seurat, group.by = c("sample"))
DimPlot(seurat, group.by = c("Phase", "sex","age_bin"))
DimPlot(seurat, group.by = c("manuscript_anno1"), label = TRUE)
DimPlot(seurat, group.by = c("manuscript_anno2"))
DimPlot(seurat, group.by = c("manuscript_anno3"))        

FeaturePlot(seurat, min.cutoff = "q1", max.cutoff = "q99", features = c("MYH1","MYH2","MYH7","MYH7B"))

## need broader annotation for input to MOFAcell
seurat$annotation_level1 <- plyr::revalue(seurat$manuscript_anno1,
                                          c("Adipocyte" = "Adipocyte",
                                            "VenEC" = "Endothelial",
                                            "CapEC" = "Endothelial",
                                            "ArtEC" = "Endothelial",
                                            "LymphEC" = "Endothelial",
                                            "MuSC" = "MuSC",
                                            "SMC" = "SMC",
                                            "Pericyte" = "Pericyte",
                                            "NK-cell" = "Lymphoid",
                                            "T-cell" = "Lymphoid",
                                            "B-cell" = "Bcell",
                                            "B-plasma" = "Bcell",
                                            "FB" = "FAP",
                                            "Macrophage" = "Myeloid",
                                            "cDC2" = "Myeloid",
                                            "cDC1" = "Myeloid",
                                            "Monocyte" = "Myeloid",
                                            "Neutrophil" = "Myeloid",
                                            "Eosinophil" = "Myeloid",
                                            "mSchwann" = "Schwann",
                                            "nmSchwann" = "Schwann",
                                            "Mast" = "Mast",
                                            "MF-II" = "FastSKM",
                                            "MF-IIsc(fg)" = "FastSKM",
                                            "MF-IIsn(fg)" = "FastSKM",
                                            "Hyb" = "FastSKM",
                                            "Specialised MF" = "SlowSKM",
                                            "MF-I" = "SlowSKM",
                                            "MF-Isc(fg)" = "SlowSKM",
                                            "MF-Isn(fg)" = "SlowSKM",
                                            "RBC" = "Erythrocyte",
                                            "Tenocyte" = "Tenocyte",
                                            "PnFB" = "Tenocyte",
                                            "Mesothelium" = "Tenocyte",
                                            "EnFB" = "Tenocyte",
                                            "pDC" = "Myeloid"))        
DimPlot(seurat, group.by = "annotation_level1", label = TRUE)        

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Save Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, "data/processed/kedlian_human_relabel.qs")


        
        
        
        
        
        
        
        