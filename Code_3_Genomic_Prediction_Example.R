#'-------------------------------------------------------------------------------------------------#
#' Title          : Genomic Prediction of GxY across years (by TPE)
#' Subtitle       : Using two yr to predict a next
#' Created at     : 2021-12-08
#' Last update at : 2022-08-10 
#' Goal           : Demonstrate the use of PLS 1 and PLS 2 for training GP models using EnvRtype
#'                  
#' Author         : Germano Costa-Neto <gmc222@@cornell.edu>
#'-------------------------------------------------------------------------------------------------#
rm(list=ls())

# Source
require(tidyverse)
require(plyr)
require(reshape2)
require(EnvRtype)
require(AGHmatrix)

source('./src/additional_src.R')




home.dir = getwd()
setwd(home.dir)

#'-------------------------------------------------------------------------------------------------#
# H matrix paremeters (see Martini et al. 2018)
#'-------------------------------------------------------------------------------------------------#
tau_G   = 1 
omega_G = 1

tau_e   = 1
omega_e = 1

#'-------------------------------------------------------------------------------------------------#
# Parameters for BGGE
#'-------------------------------------------------------------------------------------------------#
ite    = 1E3    # 15E3 as used in the paper
burn   = 2E2    # 5E3 as used in the paper
thi    = 10
tol    = 1e-10
myseed = 1112

#'-------------------------------------------------------------------------------------------------#
## Data sets ####
#'-------------------------------------------------------------------------------------------------#

# phenotypic and envirotyping data

pheno_data_analysis <- readRDS('./data/MET_data.rds') %>% droplevels()


# G matrix

G_matrix <- readRDS('./data/G_matrix.rds')  # G-matrix

(TPEY = sort(unique(pheno_data_analysis$TPE_year))) # TPE x Year combinations



#saveRDS(object = gger,file = 'PLS1_Rmatrix.rds')

#'-------------------------------------------------------------------------------------------------#
## definition of the training set ####
#'-------------------------------------------------------------------------------------------------#


trn = TPEY[1:20]
tsn = TPEY[4:23]



# ids
names=c()
#for(i in 1:length(trn)) names[i]    = paste0(names(trn)[[i]],'_to_',names(tsn)[[i]])
names    = paste0(trn,'_to_',tsn)
names

myTPEs <- colsplit(colsplit(names,pattern = '_to_',c('y1','y2'))[,1],pattern = '_',c('y1','ye'))[,2]
#names <- names[6]

#'-------------------------------------------------------------------------------------------------#
### RUNNING MODELS IN PARALELL ########
#'-------------------------------------------------------------------------------------------------#

require(doParallel)
require(foreach)

cl <- makeCluster(4) # number of clusters
registerDoParallel(cl)

#'-------------------------------------------------------------------------------------------------#

my_genomic_prediction_results = 
         foreach(SET = 1:length(names), .combine = "rbind",.errorhandling="pass",.packages = c('tidyverse','reshape2','plyr','EnvRtype')) %:%
         foreach(MYMODEL = 1:6, .combine = "rbind",.errorhandling="pass",.packages = c('tidyverse','reshape2','plyr','EnvRtype')) %dopar%
  {
    
    
    #'-------------------------------------------------------------------------------------------------#
    # phenotypes
    #'-------------------------------------------------------------------------------------------------#
    
    idpheno = c(1:3,8:10,109)
    my_met  = names(pheno_data_analysis)[idpheno]
    
    myPheno = pheno_data_analysis[,my_met]       %>% 
      filter(TPE_year %in% c(trn[SET],tsn[SET])) %>%
      droplevels()
    
    myPheno = data.frame(env=myPheno$environment,gid=myPheno$GID,value=myPheno$BLUE,TPE=myPheno$TPE,year=myPheno$Cycle)
    
    #'-------------------------------------------------------------------------------------------------#
    # G-matrix
    #'-------------------------------------------------------------------------------------------------#
    
    gid = levels(myPheno$gid)
    
    GRM <- G_matrix[which(row.names(G_matrix) %in% gid),which(colnames(G_matrix) %in% gid)]
    
    #'-------------------------------------------------------------------------------------------------#
    # Phenotypic data 
    #'-------------------------------------------------------------------------------------------------#
    
    myPheno = pheno_data_analysis[,my_met]%>% 
      filter(TPE_year %in% c(trn[SET],tsn[SET])) %>%
      filter(GID %in% rownames(GRM)) %>%
      droplevels()
    
    # checking
    myPheno = data.frame(env=myPheno$environment,gid=myPheno$GID,value=myPheno$BLUE,year=myPheno$Cycle,TPE_year=myPheno$TPE_year)
    gid = levels(myPheno$gid)
    
    # matching G-matrix and phenotypic data
    GRM = G_matrix[which(row.names(G_matrix) %in% gid),which(colnames(G_matrix) %in% gid)]
    
    
    #'-------------------------------------------------------------------------------------------------#
    # Envirotyping data
    #'-------------------------------------------------------------------------------------------------#
    
    idcov   = c(16:30,36:40,46:75,81:82,91:92,95,101,104:105,108) # covariables used in this study
    my_covs = names(pheno_data_analysis)[idcov]
    
    
    environmental_covariables = 
      pheno_data_analysis %>%
      filter(TPE_year %in% c(trn[SET],tsn[SET])) %>%
      melt(measure.vars = my_covs) %>%
      acast(environment~variable,mean,value.var = 'value') %>%
      scale(center = T,scale = T)
    
    environmental_covariables <-  environmental_covariables[,!is.nan(apply(environmental_covariables,2,sum))]
    
    saveRDS(object =  environmental_covariables,file = paste0('W_matrix_',names[SET]))
    
    W.matrix <-  environmental_covariables
    
    #'-------------------------------------------------------------------------------------------------#
    # Including Reaction-Norm coefficients - R-matrix
    #'-------------------------------------------------------------------------------------------------#
    
    Rmatrix  =   
      readRDS('./data/PLS1_Rmatrix.rds')  %>% # importing R-matrix
      melt(id.vars=c('GID','r')) %>% # GID = genotype ID, r = PLS 1 model accuracy
      filter(GID %in% rownames(GRM)) %>% droplevels() %>%  # filtering the genotypes of interest
      acast(GID~variable,value.var = 'value')             # genotype x coefficient
    
    myPheno$year = as.numeric(myPheno$year)
    (myYear= sort(unique(myPheno$year)))
    
    
    (myENVs <- as.character(levels(myPheno$env)))
    (myLOC <- colsplit(myENVs,pattern = '_',c('loc','year')))
    myNA <- which(myLOC$year %in% max(myYear))
    
    gid <- myPheno %>% 
      filter(!env %in% myENVs[myNA]) %>% 
      droplevels()
    
    gid <- as.character(levels(gid$gid))
    
    #'-------------------------------------------------------------------------------------------------#
    # Genetic similarity due to the shared reaction-norms
    #'-------------------------------------------------------------------------------------------------#
    id_Rmatrix = which(row.names(Rmatrix) %in% gid)
    GRM_GSP    = GK_Kernel(X = list(x = Rmatrix[id_Rmatrix ,]))[[1]]+diag(0.01,nrow = nrow(Rmatrix[id_Rmatrix ,]))
    
    #'-------------------------------------------------------------------------------------------------#
    # Environmental Similarity Matrices (ERM)
    #'-------------------------------------------------------------------------------------------------#
    
    # ERM based on W-matrix (models M02 and M03)
    ERM_M02  = EnvRtype::env_kernel(env.data = W.matrix,gaussian = F )[[2]] # linear
    ERM_M03  = EnvRtype::env_kernel(env.data = W.matrix,gaussian = T )[[2]] # nonlinear
    
    ERM_M01 =  ERM_M02 *0+diag(1,nrow = nrow(ERM_M02 ))     # dummy  identitiy matrix
    
    # ERM based on S-Matrix  and PLS 2
    
    env_weights_matrix =  readRDS(paste0('./data/PLS_Weights_for_each_site_TPE_',myTPEs[SET]))
    
    # organizing the kernel for location x location
    K_st = GK_Kernel(X = list(env_weights_matrix ))[[1]]  
    
    K_st = melt(K_st,varname=c('siteL','siteC'),value.name = 'site')
    
    K_wt = melt(ERM_M02,varname=c('envL','envC'),value.name = 'w')
    K_wt = data.frame(K_wt,
                      colsplit(K_wt$envL,pattern = '_',c('siteL','yL')),
                      colsplit(K_wt$envC,pattern = '_',c('siteC','yC')))
    
    ERM_M04  = merge(K_st,K_wt,by=c('siteL','siteC')) %>% acast(envL~envC,value.var = 'site')
    
    # ERM based on the S-matrix (PLS 2 outcomes)
    ERM_M04 = ERM_M04+diag(0.01,nrow = nrow(ERM_M04))
    
    rm(K_st,K_wt) #
    
  #  saveRDS(object = ERM_M04,file =  paste0('ERM_M04_',names[SET]))
    
    #'-------------------------------------------------------------------------------------------------#
    # Computing Hmatrix -- Merging Rmatrix and Gmatrix 
    #'-------------------------------------------------------------------------------------------------#
    
    GRM_Hmatrix  = doHkernel(R = GRM_GSP,G = GRM)
    
    #'-------------------------------------------------------------------------------------------------#
    # Computing GxE kernel merging past and future GxE
    #'-------------------------------------------------------------------------------------------------#
    #'
    pastGE = kronecker(ERM_M03,GRM_GSP,make.dimnames = T) # W-matrix x R-matrix
    expGE  = kronecker(ERM_M04,GRM    ,make.dimnames = T) # S-matrix x G-matrix
    
    my_gid_env = paste0(myPheno$env,':',myPheno$gid)
    
    # attention: to reduce dimensionality, remove NAs gid-env
    id = which(colnames(pastGE) %in% my_gid_env)
    pastGE <- pastGE[id,id]
    
    id = which(colnames(expGE) %in% my_gid_env)
    expGE <- expGE[id,id]
    
    # weighted GE matrix (Genotype x Environment relationship matrix, GERM)
    GERM <- doGEkernel(N=expGE, M=pastGE,tau = tau_G ,omega = omega_G) # omega=0 geralmente melhor
    dim( GERM)

    #'-------------------------------------------------------------------------------------------------#
    # kernel for S-matrix (S-W matrix)
    #'-------------------------------------------------------------------------------------------------#
    # PLOS: a kernel merging W and S matrix
    
    KWS = doSkernel(W = ERM_M03,S = ERM_M04)

    #'-------------------------------------------------------------------------------------------------#
    # Statistical Models -- get_kernel() ####
    #'-------------------------------------------------------------------------------------------------#
    
    # M01: Conventional multi-environment GBLUP
    m1 <- EnvRtype::get_kernel(
      K_G = list(G = GRM),
      K_E = list(E = ERM_M01),
      model = 'RNMM',  # denotes kronecker between E and G
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value')
    
    # M02: Reaction-norm GBLUP with a linear kernel for W-matrix (Ω)
    m2 <- EnvRtype::get_kernel(
      K_G = list(G = GRM),
      K_E = list(E = ERM_M02),
      model = 'RNMM',
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value')
    

    # M03: Reaction-norm GBLUP with a nonlinear Gaussian kernel for W-matrix (γ)
    m3 <- EnvRtype::get_kernel(
      K_G = list(G = GRM),
      K_E = list(E = ERM_M03),
      model = 'RNMM',
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value')
    
    #'-------------------------------------------------------------------------------------------------#
    # Integrating EPA outcomes in Predictive Models
    #'-------------------------------------------------------------------------------------------------#
    
    # M04: Reaction-norm GBLUP with environmental weights (Φ) from EPA
    m4 <- EnvRtype::get_kernel(
      K_G = list(G = GRM),
      K_E = list(S = ERM_M04),
      model = 'RNMM',
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value')
    
    # M05: Reaction-norm GBLUP with genotype-specific factors (R-matrix) from EPA
    m5 <- EnvRtype::get_kernel(
      K_G = list(G = GRM,H = GRM_Hmatrix),
      K_E = list(S = ERM_M04),
      model = 'RNMM',
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value') 
    
    
    m5 <- m5[- which(names(m5) %in% c('KG_G_H','KGE_GS'))]
    
    # OBS: I tested not to remove this kernel and it seems the results not change too much
    # that is why for the final model I decided to remove.
    # I suggest you to try with your own data and see which works better
    
    
    # M06: Reaction-norm GBLUP with single-step G×E kernel from EPA
    m6 <- EnvRtype::get_kernel(
      K_G = list(G = GRM),
      K_E = list(S = ERM_M04),
      model = 'RNMM',
      data = myPheno,
      env = 'env',gid = 'gid',y = 'value') 
    
    # replace the GxE kernel for the GERM
    
    currentGE = m6$KGE_GS$Kernel
      
    ids = which(rownames(GERM) %in% rownames( currentGE ))
    m6$KGE_GS$Kernel =  GERM[ids,ids]
    
    
    #'-------------------------------------------------------------------------------------------------#
    # cross-validation
    
    ts          <- which(myPheno$TPE_year %in% tsn[SET])
    myPheno$y     <- myPheno$value
    myPheno$y[ts] <- NA

    # png(filename = paste0('plot_',names_l[SET],'.png'),width = 900,height = 900)
    #superheat::superheat(acast(myPheno,gid~env,value.var = 'y'),left.label.text.size = 3,bottom.label.text.size = 3,bottom.label.text.angle = 90)
    #dev.off()
    
    Model_list <- list(m1,m2,m3,m4,m5,m6)
    Models <- c(paste0('M0',1:6))
    
    myPheno <- droplevels(myPheno)

    yhat = data.frame(obs=myPheno$value,pred=NA,gid=myPheno$gid, env=myPheno$env,
                      TPE=myPheno$TPE,pop=myPheno$y,TPE_year=myPheno$TPE_year,set=names[SET],
                      n_loc=nrow(ERM_M01), n_predict=length(myNA),Model = Models[MYMODEL])
    
    yhat$pop[!is.na(yhat$pop)] = 'training'
    yhat$pop[is.na(yhat$pop)] = 'testing'

      set.seed(myseed)
      fit <- EnvRtype::kernel_model(
                        data = myPheno,env = 'env',gid = 'gid',y = 'y',tol = tol,
                        random = Model_list[[MYMODEL]],
                        iterations = ite,burnin = burn,thining = thi)
        

      
      yhat$pred = fit$yHat 
      

    return(yhat)
  }

stopCluster(cl)

head(my_genomic_prediction_results )

# remember to save tour results
saveRDS(object = my_genomic_prediction_results ,file = 'myResults')




## OBS: how to compute the predictive ability and other metrics


my_genomic_prediction_results  %>%  
  filter(pop %in% 'testing') %>% 
  ddply(.(Model,TPE_year),summarise,cor=cor(obs,pred)) %>%
  ddply(.(Model),summarise,pa=round(median(cor),3))


my_genomic_prediction_results %>%  
  filter(pop %in% 'testing') %>% 
  ddply(.(Model,env,TPE_year),summarise,cor=cor(obs,pred))



## by genotype
my_genomic_prediction_results  %>%  
  # filter the testing set
  filter(pop %in% 'testing') %>% 
  # computing the correlations at genotype level
  ddply(.(Model,gid),summarise,cor=cor(obs,pred,method = 'spearman',use = 'pairwise.complete.obs')) %>% 
  
  # now lets do the plot ;-)
  ggplot()+geom_tile(aes(y=reorder(gid,-cor),x=Model,fill=cor),colour='white',size=0.02)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colours = rainbow(10),na.value = 'white',
                       breaks=c(0.2,.4,0.6,0.8,1),labels=c('0.20                 ','0.40  ','0.60','0.80','1.00'),
                       limits=c(0.2,1))+ 
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5))+
  #facet_grid(Model2~.,scales = 'free',space = 'free')+
  ylab('Genotype-specific\n predictive ability')+
  xlab(' ')+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.title.x.bottom  = element_text(size=13),
        strip.text.y = element_text(size=13,angle=360,face='bold'),
        axis.title.x.top     = element_text(size=15,face = 'bold'),
        axis.title.x = element_text(size=13,face='bold'),
        plot.title = element_text(size=15,face='bold',hjust = 0.5),
        # title = element_text(size=13,face='bold',hjust = 0.8),
        legend.title = element_text(size=13,face='bold'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.text = element_text(size=9),#legend.position = 'left',
        panel.background = element_rect(fill='white'),
        axis.text = element_text(size=11))+coord_flip()+
  labs(fill='Spearman\n   Rank')+
  theme_bw()+
  theme(axis.title.x.bottom  = element_text(size=13),
        axis.title.x.top     = element_text(size=15,face = 'bold'),
        #  panel.border = element_rect(color = "black",
        #                                    fill = NA,
        #                                   size = 1),
        axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.title = element_text(size=14,face='bold'),
        #     title = element_text(colour='black',fill='gray'),
        #  strip.background = element_rect(fill='gray98',colour='black',size=1),
        strip.text = element_text(size=15),
        plot.title = element_text(size=15,face='bold',hjust = 0.5),
        # title = element_text(size=13,face='bold',hjust = 0.8),
        legend.title = element_text(size=13,face='bold'),
        #   axis.line = element_line(colour = "black"),
        #  axis.line.y.right   = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        #  panel.grid.major.x =   element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=9),legend.position = 'none',
        #  panel.background = element_rect(fill='gray93'),
        axis.text = element_text(size=12))



