# Load data
###########
library(corrplot)
library(RColorBrewer)
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]

# Create correlation plot
###########################
vars<- c("TMB","neoAgB1","neoAgB2","neoAgBnorm","pMHC1","pMHC2","pMHCnorm","MGBS1","MGBS2","MGBSnorm","MHC1","MHC2","heterogeneity","ploidy","purity","MHC1_zyg","MHC2_zyg","PD1","B2M","CYT")
cor_df<- ICB_data[,vars]
colnames(cor_df)<- c("TMB","abs. MHC-I neoAgB","abs. MHC-II neoAgB","abs. diff. neoAgB", "norm. MHC-I neoAgB","norm. MHC-II neoAgB","norm. diff. neoAgB", "MGBS-I","MGBS-II","MGBS-d","MHC-I expr.","MHC-II expr.","heterogeneity","ploidy","purity","MHC-I zygosity","MHC-II zygosity","PD1 (CD274) expr.","B2M expr.","CYT")
  
cor_df<- na.omit(cor_df)
# cor_df<- cor_df[cor_df$correctBiopsy==T,] # IMPORTANT????
# cor_df$correctBiopsy<- NULL
cor_data<- cor(cor_df)

cor_data_p<- cor.mtest(na.omit(cor_df), conf.level = 0.95)

# Not ggplot, use cowplod & recordplot()
library(cowplot)
corrplot(cor_data, method="ellipse", type= "lower", order="hclust", p.mat = cor_data_p$p, sig.level = 0.05, addrect = 0, insig='blank',
         number.cex = 0.8, tl.col = 'black', tl.srt = 45, tl.cex = 0.75)
p1_recorded <- recordPlot()
p<- ggdraw(p1_recorded)

# Save
#########
save.image(file="results/data/manuscript_correlation_analysis.RData")
saveRDS(p, file="results/data/manuscript_correlation_analysis.rds")


