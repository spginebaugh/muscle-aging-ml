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
source("utils/mofa_analysis_functions.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perez_bulk_dat <- read.csv("data/processed/perez_bulk_normalized.csv", row.names = 1)
perez_bulk_meta <- read.csv("data/processed/perez_bulk_metadata.csv", row.names = 1)

input_dir <- file.path("data/mofa_inputs/")
model_dir <- file.path("models/")
model_file <- "mofa_model1.qs"

sample_id_column <- "sample"

sample_meta <- read.csv(paste0(input_dir, "sample_metadata.csv"))
mofa_model <- qread(paste0(model_dir, model_file))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          project factor onto bulk                        ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## get factor loadings 
test_variable <- "age_bin"
aging_factor <- select_important_factor(mofa_model, sample_meta, sample_id_column, test_variable)
factor_loadings <- get_geneweights(model = mofa_model, factor = aging_factor)

## project onto bulk data 
factor_score_perez <- project_factor_on_bulk(perez_bulk_dat, factor_loadings, perez_bulk_meta)

## organize
factor_score_perez <- factor_score_perez[order(rownames(factor_score_perez)), sort(colnames(factor_score_perez))]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    save                                  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(factor_score_perez, "data/processed/factor_score_perez.csv", row.names = TRUE)

