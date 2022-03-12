#Applying MOFA+ to the Chronic Lymphocytic Leukemia cohort

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

DNA <- read.csv("LIHC_preprocessed_DNAMeth.csv", header=TRUE)
CNA <- read.csv("LIHC_preprocessed_CNA.csv", header=TRUE)
RNA <- read.csv("LIHC_preprocessed_RNASeq.csv", header=TRUE)

DNA <- data.matrix(t(DNA))
CNA <- data.matrix(t(CNA))
RNA <- data.matrix(t(RNA))

data<- list(DNA,CNA,RNA)
names(data) <- c("DNA", "CNA", "RNA")


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

outfile = file.path(getwd(),"LIHC_bulkMultiOmics_model.hdf5")
LIHC_bulkMultiOmics_model <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)

#Overview of the trained MOFA model
slotNames(MOFAobject)
names(MOFAobject@data)


#filepath <- getwd()
#LIHC_bulkMultiOmics_model <- load_model(filepath)


#Extract factors
# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(LIHC_bulkMultiOmics_model, factors = "all")
lapply(factors,dim)
# Concatenate groups 
#factors <-do.call("rbind",factors)

#Extract weights
# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(LIHC_bulkMultiOmics_model, views = "all", factors = "all")
lapply(weights,dim)

#Extract data
# "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
#data <- get_data(LIHC_bulkMultiOmics_model)
#lapply(data, function(x) lapply(x, dim))[[1]]

write.csv(factors, 'LIHC_bulkMultiOmics_MOFA2_embedding.csv')
