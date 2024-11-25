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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Import data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_dir <- file.path("data/mofa_inputs/")
output_dir <- file.path("models/")

multiview_dat <- read.csv(paste0(input_dir, "multiview_dat.csv"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              MOFA parameters                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
num_factors <- 7
spikeslab_weights <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               Fit MOFA model                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mofa_obj <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(mofa_obj)
train_opts <- MOFA2::get_default_training_options(mofa_obj)
model_opts <- MOFA2::get_default_model_options(mofa_obj)

## FALSE avoids less sparse gene weights
model_opts$spikeslab_weights <- spikeslab_weights

## Define the number of factors needed (based on trial and error, unfortunately)
model_opts$num_factors <- num_factors

## create MOFA model
mofa_obj <- MOFA2::prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

## run MOFA2
mofa_model <- MOFA2::run_mofa(mofa_obj, outfile = NULL, save_data = FALSE, use_basilisk = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 save files                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(mofa_model, paste0(output_dir, "mofa_model1.qs"))
