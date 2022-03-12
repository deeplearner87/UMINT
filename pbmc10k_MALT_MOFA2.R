library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

rna <- read.csv("pbmc10k_MALT_rna_scaled.csv", header=TRUE, row.names=1)
protein <- read.csv("pbmc10k_MALT_adt_scaled.csv", header=TRUE, row.names=1)
rna <- data.matrix(rna)
protein <- data.matrix(protein)
data<- list(rna,protein)
names(data) <- c("RNA", "Protein")


#Create the MOFA object
MOFAobject <- create_mofa(data)

#Plot data overview
plot_data_overview(MOFAobject)

#Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

#Model options
model_opts <- get_default_model_options(MOFAobject)
#model_opts$num_factors <- 64
#model_opts$num_factors <- 15
model_opts

#Training options
train_opts <- get_default_training_options(MOFAobject)
#train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$verbose <- TRUE
train_opts

#Training...
#Prepare MOFa object
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

outfile = file.path(getwd(),"pbmc10k_MALT_model.hdf5")
#model <- run_mofa(MOFAobject, outfile)
pbmc10k_MALT_model <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)

#Overview of the trained MOFA model
slotNames(MOFAobject)
names(MOFAobject@data)


#filepath <- getwd()
#pbmc10k_MALT_model <- load_model(filepath)


#Extract factors
# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(pbmc10k_MALT_model, factors = "all")
lapply(factors,dim)
# Concatenate groups 
#factors <-do.call("rbind",factors)

#Extract weights
# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(pbmc10k_MALT_model, views = "all", factors = "all")
lapply(weights,dim)

#Extract data
# "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
#data <- get_data(model_pbmc10K)
#lapply(data, function(x) lapply(x, dim))[[1]]

write.csv(factors, 'pbmc10k_MALT_MOFA2_embedding.csv')
