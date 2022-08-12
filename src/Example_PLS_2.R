#'---------------------------------------------------------------------------------------
#'
#' Title......: EPA Analysis using PLS 2 algorithm
#' Goal.......: Compute environmental weights and model similarities for future locations
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

pheno_data_analysis <- readRDS('./data/MET_data.rds') %>% droplevels()

idcov   = c(16:30,36:40,46:75,81:82,91:92,95,101,104:105,108) # covariables used in this study
my_covs = names(pheno_data_analysis)[idcov]

S_matrix = 
  pheno_data_analysis %>%
  melt(measure.vars = my_covs) %>%
  ddply(.(Loc_no,variable),summarise, # sample by location
        q10=quantile(value,.1,na.rm=T), # you can use the quantiles you want.
        q50=quantile(value,.5,na.rm=T), # here I used q10, q50 and q90
        q90=quantile(value,.9,na.rm=T)) %>%
  melt(measure.vars = c('q10','q50','q90'),variable.name='qt') %>%
  acast(Loc_no~variable+qt,mean,value.var = 'value') %>%
  scale(center = T,scale = T)

superheat::superheat(S_matrix,pretty.order.rows = T,pretty.order.cols = T)


dim(S_matrix)

S0_matrix = 
  pheno_data_analysis[,c(2:3,9,10)] %>%
  melt(measure.vars = 'BLUE') %>%
  ddply(.(Loc_no,Cycle),summarise,
        q10=quantile(value,.1,na.rm=T),
        q50=quantile(value,.5,na.rm=T),
        q90=quantile(value,.9,na.rm=T)) %>%
  melt(measure.vars = c('q10','q50','q90'),variable.name='qt') %>%
  acast(Loc_no~qt,mean,value.var = 'value') %>%
  scale(center = T,scale = T)

dim(S0_matrix)



phy0_matrix = EnvRtype::env_kernel(env.data = S0_matrix,gaussian = T)[[2]]


require(plsdepot)
min_var_exp = .85
PLS_2       = plsreg2(predictors = S_matrix,responses = phy0_matrix,comps = 25,crosval = F)
pls$expvar
plot(PLS_2)

comp_min  = which.max(PLS_2  $expvar[,4][round(PLS_2  $expvar[,4],3) <= min_var_exp ])

PLS_expvar = data.frame(PLS_2$expvar)
PLS_expvar$comp <- rownames(PLS_expvar)

PLS_expvar = 
  PLS_expvar %>% melt() %>% 
  filter(variable %in% c('R2Xcum','R2Ycum')) %>% droplevels()

PLS_expvar$ncomp <- as.numeric(gsub(PLS_expvar$comp,pattern = 't',replacement = ''))

PLS_expvar %>% ggplot(aes(x=ncomp,y=value,colour=variable,fill=variable))+
  geom_point()+geom_line()+geom_hline(yintercept = min_var_exp)+
  geom_vline(xintercept = comp_min)+
  scale_y_continuous(name = 'Explained Variance %',labels = scales::percent)


environmental_weights = PLS_2$std.coefs

saveRDS(file = './data/environmental_weights_matrix.rds',object = environmental_weights)

# matrix of relationship for locations based on the environmental_weights
phy_matrix = EnvRtype::env_kernel(env.data = environmental_weights ,gaussian = T)[[2]]

superheat::superheat(phy_matrix,pretty.order.rows = T,pretty.order.cols = T)


