#'---------------------------------------------------------------------------------------
#'
#' Title......: Use of EPA-based outputs to group environments
#' Goal.......: Identify the pivot locations from the environmental weights
#' Author.....: Germano Costa Neto
#' Created at.: 2021-11-11
#' Modified at: 2022-07-27
#' 
#'------------------------------------------------------------------------------------


rm(list=ls())

require(plsdepot)
require(tidyverse)
require(reshape2)
require(plyr)
require(FactoMineR)
require(STPGA)


environmental_weights <- readRDS('./data/environmental_weights_matrix.rds')

# matrix of relationship for locations based on the environmental_weights
phy_matrix = EnvRtype::env_kernel(env.data = environmental_weights ,gaussian = T)[[1]]

superheat::superheat(phy_matrix,pretty.order.rows = T,pretty.order.cols = T)


PCA_phy = FactoMineR::PCA(t(environmental_weights ))
HCPC_phy = FactoMineR::HCPC(PCA_phy,nb.clust = -1)


HCPC_phy$data.clust$env = rownames(HCPC_phy$data.clust)

TS <- GenAlgForSubsetSelectionNoTest(P = phy_matrix , ntoselect = 5, # 5 groups
                                     InitPop=NULL,
                                     npop=100, nelite=5, 
                                     mutprob=.5, mutintensity = 1,
                                     niterations=200,minitbefstop=20, tabu=F,tabumemsize = 0,plotiters = T, 
                                     lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = 4)[[1]]

(TS = as.character(TS) )# training set (super-optmized) consisted by some genotypes at some environments



