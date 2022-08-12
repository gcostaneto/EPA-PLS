#'---------------------------------------------------------------------------------------
#'
#' Title......: EPA Analysis using PLS 1 algorithm
#' Goal.......: Compute genotype-specific reaction-norms
#' Author.....: Germano Costa Neto
#' Created at.: 2021-11-11
#' Modified at: 2022-07-27
#' 
#'---------------------------------------------------------------------------------------



rm(list=ls())
home.dir=getwd()

require(plsdepot)
require(tidyverse)
require(reshape2)
require(plyr)


#'-------------------------------------------------------------------------------------------------#
# Data sets
#'-------------------------------------------------------------------------------------------------#


pheno_data_analysis <- readRDS('./data/MET_data.rds') %>% droplevels()

idcov   = c(16:30,36:40,46:75,81,91:92,95,101,104:105,108) # covariables used as example in this study
my_covs = names(pheno_data_analysis)[idcov]


idpheno = c(1:3,8:10,109) # id of the environment, locations, genotypes, traits
my_met  = names(pheno_data_analysis)[idpheno]

Mj = 
  pheno_data_analysis[,c(my_met,my_covs)] %>% 
  droplevels() %>% 
  ddply(.(environment),mutate,Mj = round(mean(BLUE,na.rm=T),3))
Mj

# Environmental-centering
Mj$GGE = round(Mj$BLUE-Mj$Mj,4)
(mygids = levels(Mj$GID))



path_output = paste0(home.dir,'/PLS 1 output')
dir.create(path_output)


min_var_exp = .90
coef=c()   # genotype-specific coefficients
gid_VIPs = c()

for(i in 1:length(mygids))
{
  cat(paste0('....doing GID N',i,' = ',mygids[i],'\n'))
  myG = Mj %>% filter(GID %in% mygids[i])
  
  W_matrix = myG[,my_covs]
  W_matrix =  W_matrix[,!is.nan(apply(W_matrix,2,sum))]
  
  pls= plsdepot::plsreg1(predictors = W_matrix ,
                  response = myG$GGE,comps = 30,crosval = F)
  
  
  # Estimate the number of LV
  comp_min = which.max(cumsum(pls$R2)[cumsum(pls$R2) <= min_var_exp ])
  
  
  # run PLS 1 again, but now with LV = comp_min
  
  pls=plsdepot::plsreg1(predictors = myG[,my_covs],
                  response = myG$GGE,comps = comp_min,crosval = F)
  
  # save VIP
  VIP_sum  = rowSums(pls$R2Xy)
  VIP_sum  = VIP_sum[-which(names(VIP_sum) %in% 'Y')]
  VIP_sum = melt(VIP_sum/max(VIP_sum))
  VIP_sum$EC = rownames(VIP_sum)
  
  # save GSP (genotype specific coefficients)
  coef=rbind(coef,data.frame(t(pls$std.coefs),GID=mygids[i],
                             r=round(cor(pls$y,pls$y.pred),3),
                             LV=length(pls$R2)))
  
  gid_VIPs = rbind(gid_VIPs,data.frame(GID=mygids[i],VIP_sum))
  
}

names(coef)

Rmatrix  =   
  coef %>% 
  melt(id.vars=c('GID','r','LV')) %>% # GID = genotype ID, r = PLS 1 model accuracy
  acast(GID~variable,value.var = 'value')             # genotype x coefficient


require(superheat)
superheat(Rmatrix,
          row.title = "Genotypes - wheat lines",
          pretty.order.rows = F,
          pretty.order.cols = T,
          row.dendrogram = F,
          col.dendrogram = F,
          grid.vline.col = "white",
          grid.hline.col = "white",
          #  title = paste0('Optmization of all TPEs with ',SVD$Ne, ' divergent locations\n',100*round(1-SVD$Ne/ncol(K_W),1),'% reduction of the TPE size'),
          #row.dendrogram = T,
          legend.width = 1,
          left.label.size = 0.2,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.4,
          bottom.label.text.angle = 90,
          # heat.pal = viridis::cividis(5),
          heat.pal = viridis::mako(100),
          legend.text.size = 13,
          # n.clusters.rows = 100,
          # left.label = 'variable',
          # n.clusters.cols = 3,
          ##  bottom.label = 'variable',
          #   X.text = round(as.matrix(a),1),X.text.col="white",
          legend.height=0.1)


head(coef)

gid_VIPs$importance = '< 80%'
imp=ddply(gid_VIPs,.(EC),summarise,imp=median(value))
gid_VIPs$importance[which(gid_VIPs$EC %in% imp$EC[which(imp$imp >=0.80)])] = '>=80%'


head(gid_VIPs)

gid_VIPs = data.frame(gid_VIPs,colsplit(gid_VIPs$EC,pattern = '__',c('EV','Dev_stage')))
head(gid_VIPs)

gid_VIPs$Dev_stage2 = factor(gid_VIPs$Dev_stage,levels = unique(gid_VIPs$Dev_stage)[c(6,2,5,1,4,3)])
levels(gid_VIPs$Dev_stage2) =  c('Static',
                                 'Stages 0 to 3','Stages 3 to 5',
                                 'Stages 5 to 6','Stages 6 to 7',
                                 'Stages 7 to 9')

dim(coef)

gid_VIPs$EF = paste0(gid_VIPs$EV,' [',gid_VIPs$Dev_stage2,']')

VIP_Ecs_single = 
  gid_VIPs %>% 
  ggplot(aes(y=reorder(EF,value),x=value,colour=importance,fill=importance),size=1.5,alpha=0.1)+
  geom_boxplot(alpha=0.3,size=0.4)+
  #geom_violin(alpha=0.3)+
  theme_classic()+
  ylab('Envirotyping Covariable\n(n=62)')+
  xlab('Variable Importance in Projection (VIP,%)')+
  geom_vline(xintercept = 0.8,colour='red')+
  scale_color_manual(values = c('gray80','royalblue'))+
  scale_fill_manual(values = c('gray80','royalblue'))+
  labs(fill='       VIP\nClassification',colour='       VIP\nClassification')+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1),labels = scales::percent)+
  theme(#axis.title.x.bottom  = element_text(size=13),
    axis.title.x.top     = element_text(size=15,face = 'bold'),
    #  panel.border = element_rect(color = "black",
    #                                    fill = NA,
    #                                   size = 1),
    #  axis.text.x = element_blank(),axis.ticks = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12),
    #     title = element_text(colour='black',fill='gray'),
    #  strip.background = element_rect(fill='gray98',colour='black',size=1),
    strip.text = element_text(size=15),
    plot.title = element_text(size=15,face='bold',hjust = 0.5),
    # title = element_text(size=13,face='bold',hjust = 0.8),
    legend.title = element_text(size=13,face='bold'),
    #   axis.line = element_line(colour = "black"),
    #  axis.line.y.right   = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.major.x =   element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=12)) 


VIP_Ecs_single = 
  gid_VIPs %>% 
  ggplot(aes(x=reorder(EF,value),y=value,colour=importance,fill=importance),size=1.5,alpha=0.1)+
  geom_boxplot(alpha=0.3,size=0.4)+
  #geom_violin(alpha=0.3)+
  theme_classic()+
  xlab('Envirotyping Covariable\n(n=61)')+
  ylab('Variable Importance\n in Projection (VIP, %)')+
  geom_hline(yintercept = 0.8,colour='red')+
  scale_color_manual(values = c('gray80','royalblue'))+
  scale_fill_manual(values = c('gray80','royalblue'))+
  labs(fill='       VIP\nClassification',colour='       VIP\nClassification')+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1),labels = scales::percent)+
  theme(#axis.title.x.bottom  = element_text(size=13),
    axis.title.x.top     = element_text(size=15,face = 'bold'),
    #  panel.border = element_rect(color = "black",
    #                                    fill = NA,
    #                                   size = 1),
    axis.text.x = element_text(angle=90,hjust=1),#,axis.ticks = element_blank(),
    axis.title = element_text(size=14,face='bold'),
    axis.text = element_text(size=12),
    #     title = element_text(colour='black',fill='gray'),
    #  strip.background = element_rect(fill='gray98',colour='black',size=1),
    strip.text = element_text(size=15),
    plot.title = element_text(size=15,face='bold',hjust = 0.5),
    # title = element_text(size=13,face='bold',hjust = 0.8),
    legend.title = element_text(size=13,face='bold'),
    #   axis.line = element_line(colour = "black"),
    #  axis.line.y.right   = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.major.x =   element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=12)) 
