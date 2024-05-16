# Load functions
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra) 

##########
# Fig. 1
##########

# Load plots
p_ls<- readRDS("results/data/manuscript_surv.rds")
p_BR_ls<- readRDS("results/data/manuscript_RECIST.rds")

# Create fig
p1<- p_ls$p_TMB +
  ggtitle("Tumor Mutation Burden")
p2<- p_ls$p_strongMHC1 +
  ggtitle("MGBS-I")
p3<- p_ls$p_strongMHC2 +
  ggtitle("MGBS-II")
p4<- p_ls$p_strongMHCnorm +
  ggtitle("MGBS-d")
p5<- p_ls$fp_multi_MB_1 + scale_x_discrete(labels=c("High MGBS-I","High TMB")) + scale_y_continuous(limits = c(0.15,3), trans='log10') + theme(axis.title = element_blank())
p6<- p_ls$fp_multi_MB_2 + scale_x_discrete(labels=c("High MGBS-II","High TMB")) + scale_y_continuous(limits = c(0.15,3), trans='log10') + theme(axis.title = element_blank())
p7<- p_ls$fp_multi_MB_norm + scale_x_discrete(labels=c("High MGBS-d","High TMB")) + scale_y_continuous(limits = c(0.15,3), trans='log10')
p8<- p_BR_ls$p_RECIST
p9<- p_ls$fp_th + theme(legend.position = "none")

p10<- plot_grid(
  p_BR_ls$p_R_ls$TMB, p_BR_ls$p_R_ls$MGBS1,
  p_BR_ls$p_R_ls$MGBS2, p_BR_ls$p_R_ls$MGBSnorm,
  align="hv",
  axis="bt"
)

p_top<- plot_grid(
  NA,p1,p2,
  p9,p3,p4,
  ncol=3,
  labels = c("a", "b", "c", "d","e","f"),
  label_size = 8
)
p_top

p_fp<- plot_grid(
  p5,p6,p7,
  ncol=1,
  align="hv",
  axis="bt"
)

p_bottom<- plot_grid(
  p_fp,p8,p10,
  labels = c("g", "h","i"),
  ncol=3,
  rel_widths = c(4,5,4),
  label_size = 8
)

p<- plot_grid(
  p_top,
  p_bottom,
  ncol = 1,
  rel_heights = c(1.5,1)
)

ggplot2::ggsave(paste0("results/figs/manuscript_fig1.pdf"), p, width = 178, height = 265/1.5, units = "mm")

# Some numbers
p_ls$HR_df[p_ls$HR_df$class=="Hi"&p_ls$HR_df$th==0.5,]
# cond class  th        HR      ci_l      ci_h           p n1   n th_label class_label cond_label
# 1     TMBclass    Hi 0.5 0.4590435 0.2833313 0.7437262 0.001563616 72 144      50%      Strong
# 3    MHC1class    Hi 0.5 0.6339154 0.3962685 1.0140820 0.057219842 72 144      50%      Strong
# 5    MHC2class    Hi 0.5 1.7439020 1.0771923 2.8232602 0.023668192 70 141      50%      Strong
# 7 MHCnormclass    Hi 0.5 0.4600508 0.2816971 0.7513274 0.001919444 70 141      50%      Strong

#########
# Fig. 2
#########

p_ls<- readRDS("results/data/manuscript_surv.rds")

p1<- plot_grid(
  NA,p_ls$p_ipi_vs_noIpi,
  ncol = 2
)

p2<- plot_grid(
  p_ls$p_noIpi_ls$highTMB, p_ls$p_ipi_ls$highTMB, p_ls$p_noIpi_ls$strongMHC1, p_ls$p_ipi_ls$strongMHC1, 
  p_ls$p_noIpi_ls$strongMHC2, p_ls$p_ipi_ls$strongMHC2, p_ls$p_noIpi_ls$strongMHCnorm, p_ls$p_ipi_ls$strongMHCnorm, 
  ncol=4) 

p3<- p_ls$fp_ipi_multi_2 + scale_y_continuous(limits = c(0.15,7), trans='log10')

p13<- plot_grid(
  p1, p3,
  ncol=2,
  labels = c("a", "c"),
  label_size = 8
)

p<- plot_grid(
  p13,
  p2,
  rel_heights = c(1,2),
  ncol=1,
  labels = c(NA,"b"),
  label_size = 8
)

ggplot2::ggsave(paste0("results/figs/manuscript_fig2.pdf"), p, width = 178, height = 0.5*265, units = "mm")

#########
# Fig. 3
#########

# Load plots
p_ls<- readRDS("results/data/manuscript_multivariate.rds")

# ROC & surv
p_ROC_surv<- plot_grid(
  p_ls$LR_model_ls$'TRUE'$p_ROC, p_ls$LR_model_ls$'TRUE'$p_surv,NA,
  ncol=1,
  labels = c("a","c"),
  label_size = 8
)

# Barplot
p_bp<- plot_grid(
  p_ls$LR_model_ls$'TRUE'$AUC_bp,
  labels="b",
  label_size = 8
)

# Plot
p<- plot_grid(
  p_ROC_surv, p_bp,
  ncol = 2,
  rel_widths = c(1,2)
)

ggplot2::ggsave("results/figs/manuscript_fig3.pdf", p, width = 178, height = 265/2, units = "mm")

#########
# Fig. 4
#########

# Load plots
p_ls<- readRDS("results/data/manuscript_DGE_analysis.rds")

p_right<- plot_grid(
  p_ls$p_RS_ls$highTMB, p_ls$p_RS_ls$strongMHC2,
  ncol=1
)

p<- plot_grid(
  p_ls$p_fGSEA + xlim(0,23) + ylim(0,32), p_right, 
  ncol=2,
  rel_widths = c(2,1)
)

ggplot2::ggsave("results/figs/manuscript_fig4.pdf", p, width = 3/4*178, height = 265/3, units = "mm")

###########################
# Fig. S1: power analysis
###########################

# Load plots
p_power<- readRDS("results/data/manuscript_AF_power_analysis.rds")

# Create plot
p<- plot_grid(
  p_power$SS + geom_vline(xintercept = 1, linetype=2) + geom_hline(yintercept = 1e+04, linetype=2),
  p_power$AF,
  ncol = 1
)

# Create figure
ggplot2::ggsave(paste0("results/figs/manuscript_figS1_power.pdf"), p, width = 178/2, height = 265/3, units = "mm")

#######################################
# Fig. S2: KM for different thresholds
####################################### 

p_ls<- readRDS("results/data/manuscript_surv.rds")

p<- plot_grid(
  p_ls$p_th_ls$TMBclass$`0.75` + ggtitle("Tumor Mutation Burden") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$TMBclass$`0.9` + ggtitle("Tumor Mutation Burden") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$TMBclass$`0.95` + ggtitle("Tumor Mutation Burden") + theme(legend.background=element_blank()),
  p_ls$p_th_ls$MHC1class$`0.75` + ggtitle("MGBS-I") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHC1class$`0.9` + ggtitle("MGBS-I") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHC1class$`0.95` + ggtitle("MGBS-I") + theme(legend.background=element_blank()),
  p_ls$p_th_ls$MHC2class$`0.75` + ggtitle("MGBS-II") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHC2class$`0.9` + ggtitle("MGBS-II") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHC2class$`0.95` + ggtitle("MGBS-II") + theme(legend.background=element_blank()),
  p_ls$p_th_ls$MHCnormclass$`0.75` + ggtitle("MGBS-d") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHCnormclass$`0.9` + ggtitle("MGBS-d") + theme(legend.background=element_blank()), 
  p_ls$p_th_ls$MHCnormclass$`0.95` + ggtitle("MGBS-d") + theme(legend.background=element_blank()),
  ncol = 3
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 2",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

ggplot2::ggsave(paste0("results/figs/manuscript_figS2_surv_th.pdf"), p, width = 178, height = 265, units = "mm")

################################# 
# Fig. S3: Correlation analysis
################################# 

# Load plots
p<- readRDS("results/data/manuscript_correlation_analysis.rds")

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 3",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

# Create figure
ggplot2::ggsave(paste0("results/figs/manuscript_figS3_correlation.pdf"), p, width = 178, height = 265, units = "mm")

################################# 
# Fig. S4: Validation studies
################################# 

p_ls<- readRDS("results/data/manuscript_validation.rds")

# Forest plot
# Adjust sizes facet
gt<- ggplot_gtable(ggplot_build(p_ls$p_val_ls$fp))
# gtable_show_layout(gt)
gt$heights[10] = 1/2*gt$heights[10]
# fp<- grid.draw(gt)

# Forest plot preIpi
fp_preIpi<- p_ls$p_val_ls$fp_preIpi

# Survival plots MGBS2
p_surv<- plot_grid(
  p_ls$p_val_ls$Hugo_2016$OS$strongMHC2$all + ggtitle("Hugo et al., 2016"),
  p_ls$p_val_ls$Gide_2019$OS$strongMHC2$all + ggtitle("Gide et al., 2019"),
  p_ls$p_val_ls$Rizvi_2015$PFS$strongMHC2 + ggtitle("Rizvi et al., 2015"),
  ncol=1
)

# RECIST plots
p_RECIST<- plot_grid(
  p_ls$p_val_ls$Hugo_2016$RECIST$all + theme(legend.position = "left"),
  plot_grid(NULL, p_ls$p_val_ls$Gide_2019$RECIST$all + theme(legend.position = "left"),ncol=2,rel_widths = c(1,5)),
  p_ls$p_val_ls$Rizvi_2015$RECIST$all + theme(legend.position = "left"),
  ncol=1
)

# Create plot
p<- plot_grid(
  gt,
  fp_preIpi,
  plot_grid(p_surv, p_RECIST,ncol=3,labels=c("b","c",NA), rel_widths = c(2,3,1), label_size = 8),
  ncol = 1,
  rel_heights = c(1.5,1,3),
  labels = c("a","b",NA),
  label_size = 8
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 4",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

# Create figure
ggplot2::ggsave(paste0("results/figs/manuscript_figS4_validation.pdf"), p, width = 178, height = 265, units = "mm")

#####################################
# Fig. S5: NeoAgB & pMHC
#####################################

p_ls<- readRDS("results/data/manuscript_surv_neoAgB.rds")

p<- plot_grid(
  p_ls$fp_neoAgB + ggtitle("Absolute Neoantigen Burden") + theme(title = element_text(size=8)), p_ls$fp_pMHC+ggtitle("Normalized Neoantigen Burden") + theme(title = element_text(size=8)),
  p_ls$p_neoAgB1 + ggtitle("MHC-I NeoAgB"), p_ls$p_pMHC1 + ggtitle("norm. MHC-I NeoAgB"),
  p_ls$p_neoAgB2 + ggtitle("MHC-II NeoAgB"), p_ls$p_pMHC2 + ggtitle("norm. MHC-II NeoAgB"),
  p_ls$p_neoAgBnorm + ggtitle("diff. NeoAgB"), p_ls$p_pMHCnorm + ggtitle("norm. diff. NeoAgB"),
  ncol=2
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 5",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p_neoAgB<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

# Create figures
ggplot2::ggsave(paste0("results/figs/manuscript_figS5_neoAgB.pdf"), p_neoAgB, width = 178, height = 265, units = "mm")

#################################################
# Fig. S6-7: Multivariate analysis
###################################################

p_ls<- readRDS("results/data/manuscript_multivariate.rds")

# 1. COX
t_preIpi<- p_ls$Cox_model_ls$'TRUE'$t_uni
t_ipiNaive<- p_ls$Cox_model_ls$'FALSE'$t_uni

p1<- plot_grid(
  t_ipiNaive, t_preIpi,
  ncol=2
)

p_cox<- plot_grid(
  p1, p_ls$Cox_model_ls$fp_multi + scale_y_continuous(limits = c(0.15,12), trans='log10') , NA,
  ncol = 1,
  rel_heights = c(1,1/2,1)
)

# 2. Log reg
t_preIpi<- p_ls$LR_model_ls$'TRUE'$t_uni
p_surv_preIpi<- p_ls$LR_model_ls$'TRUE'$p_surv
p_ROC_preIpi<- p_ls$LR_model_ls$'TRUE'$p_ROC
t_ipiNaive<- p_ls$LR_model_ls$'FALSE'$t_uni
p_surv_ipiNaive<- p_ls$LR_model_ls$'FALSE'$p_surv
p_ROC_ipiNaive<- p_ls$LR_model_ls$'FALSE'$p_ROC

# ROC & surv for ipi naive
p_ROC_surv<- plot_grid(
  p_ROC_ipiNaive, p_surv_ipiNaive,NA,
  ncol=1,
  labels = c("c","e"),
  label_size = 8
)

# Barplot
p_bp<- plot_grid(
  p_ls$LR_model_ls$'FALSE'$AUC_bp,
  labels="d",
  label_size = 8
)

# Bottom plots
p_bottom<- plot_grid(
  p_ROC_surv, p_bp,
  ncol=2,
  rel_widths = c(1,2)
)

# Create plot
p_logReg<- plot_grid(
  t_ipiNaive, t_preIpi,
  ncol = 2
)

p_logReg<- plot_grid(
  p_logReg, p_bottom,
  ncol=1
  )

# Logreg
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 6",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p_logReg<- plot_grid(
  p_logReg, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

# Cox
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 7",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p_cox<- plot_grid(
  p_cox, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

# Create figures
ggplot2::ggsave(paste0("results/figs/manuscript_figS6_LogReg.pdf"), p_logReg, width = 178, height = 265, units = "mm")
ggplot2::ggsave(paste0("results/figs/manuscript_figS7_Cox.pdf"), p_cox, width = 178, height = 265, units = "mm")

# Models
summary(p_ls$LR_model_ls$'FALSE'$fw)

# Call:
#   glm(formula = R ~ heterogeneity + ploidy + MGBS1, family = binomial, 
#       data = ICB_data[ICB_data$preIpi == isPreIpi, ])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3600  -0.7905   0.4197   0.8466   1.5964  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)    -3.2056     1.6625  -1.928  0.05384 . 
# heterogeneity  -9.8270     3.7430  -2.625  0.00865 **
#   ploidy          1.1088     0.4363   2.541  0.01104 * 
#   MGBS1         152.6444    67.5380   2.260  0.02381 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 100.087  on 72  degrees of freedom
# Residual deviance:  76.092  on 69  degrees of freedom
# AIC: 84.092
# 
# Number of Fisher Scoring iterations: 5

summary(p_ls$LR_model_ls$'TRUE'$fw)

# Call:
#   glm(formula = R ~ MGBS2 + MHC2_noBatch + hasLOHB2M + heterogeneity, 
#       family = binomial, data = ICB_data[ICB_data$preIpi == isPreIpi, 
#       ])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.6740  -0.4442   0.1063   0.4512   1.5441  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.923e+00  3.072e+00  -1.277  0.20160   
# MGBS2         -4.648e+01  1.865e+01  -2.492  0.01270 * 
#   MHC2_noBatch   1.101e-03  4.081e-04   2.698  0.00697 **
#   hasLOHB2M     -8.306e+00  3.776e+00  -2.200  0.02782 * 
#   heterogeneity  1.053e+01  4.789e+00   2.199  0.02787 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 59.401  on 42  degrees of freedom
# Residual deviance: 28.949  on 38  degrees of freedom
# AIC: 38.949
# 
# Number of Fisher Scoring iterations: 7

######################################################### 
# Fig. S8: Validation studies multivariate analysis
#########################################################

p_ls<- readRDS("results/data/manuscript_validation_multivariate.rds")
p_compare_ls<- readRDS("results/data/manuscript_validation_compare_variables.rds")

p_surv<- plot_grid(
  plot_grid(p_ls$ROC + theme(plot.title = element_blank()),ggdraw(get_legend(p_ls$predict_surv$Liu_2019$all), xlim=c(1,1), ylim = c(1,1)),NULL,nrow = 1, rel_widths = c(2,2,3)),
  plot_grid(
    p_ls$predict_surv$Liu_2019$all + ggtitle("Liu et al., 2019") + theme(legend.position = "none"),
    p_ls$predict_surv$Hugo_2016$all + ggtitle("Hugo et al., 2016") + theme(legend.position = "none"),
    p_ls$predict_surv$Gide_2019$all + ggtitle("Gide et al., 2019") + theme(legend.position = "none"),
    p_ls$predict_surv$Riaz_2017$all + ggtitle("Riaz et al., 2017") + theme(legend.position = "none"),
    nrow=1
  ),
  plot_grid(
    p_ls$predict_surv$Liu_2019$'FALSE' + ggtitle("Treatment naive") + theme(legend.position = "none"),
    p_ls$predict_surv$Liu_2019$'TRUE' + ggtitle("Ipilimumab pretreated") + theme(legend.position = "none"),
    p_ls$predict_surv$Riaz_2017$'FALSE' + ggtitle("Treatment naive") + theme(legend.position = "none"),
    p_ls$predict_surv$Riaz_2017$'TRUE' + ggtitle("Ipilimumab pretreated") + theme(legend.position = "none"),
    nrow=1
  ),
  ncol=1,
  rel_heights = c(1,1,1)
)

# Get HM legend
legend<- get_legend(p_compare_ls$pred_surv + theme(legend.position = "right"))
p_legend<- ggdraw(legend)

# Plot bottom part
p_bottom<- plot_grid(
  p_compare_ls$pred_surv + theme(legend.position = "none", axis.text.x = element_text(angle=45), title = element_blank()),
  p_compare_ls$pred_R + theme(legend.position = "none", axis.text.x = element_text(angle=45), title = element_blank()),
  plot_grid(p_ls$bp_heterogeneity,p_legend,ncol=1,labels = c("c",NA)),
  ncol=3,
  rel_widths = c(2,2,2)
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 8",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p_surv, p_bottom, footer,
  ncol = 1,
  rel_heights = c(10,10,1),
  labels=c("a","b",NA)
)

# Create figure
ggplot2::ggsave(paste0("results/figs/manuscript_figS8_validation_multivariate.pdf"), p, width = 178, height = 265, units = "mm")

# Model
p_ls$model

# Call:  glm(formula = R ~ MHC2_z + MGBS2_z, family = binomial, data = ICB_data)
# 
# Coefficients:
#   (Intercept)       MHC2_z      MGBS2_z  
# 0.1395      -0.9765       0.8967  
# 
# Degrees of Freedom: 45 Total (i.e. Null);  43 Residual
# (14 observations deleted due to missingness)
# Null Deviance:	    63.68 
# Residual Deviance: 48.73 	AIC: 54.73

# Correlation Cox HRs
cor_df<- p_compare_ls$pred_surv$data
vars<- levels(cor_df$predictor)

cor_matrix<- matrix(NA, length(vars), 4, dimnames = list(vars, c("liu_naive", "liu_pre", "riaz_naive", "riaz_pre")))
for(v in vars){
  cor_matrix[v,"liu_naive"]<- as.numeric(cor_df[cor_df$predictor==v & cor_df$study=="Liu et al., 2019" & cor_df$preIpi=="Treatment naive", "HR"])
  cor_matrix[v,"liu_pre"]<- as.numeric(cor_df[cor_df$predictor==v & cor_df$study=="Liu et al., 2019" & cor_df$preIpi=="Pretreated", "HR"])
  cor_matrix[v,"riaz_naive"]<- as.numeric(cor_df[cor_df$predictor==v & cor_df$study=="Riaz et al., 2017" & cor_df$preIpi=="Treatment naive", "HR"])
  cor_matrix[v,"riaz_pre"]<- as.numeric(cor_df[cor_df$predictor==v & cor_df$study=="Riaz et al., 2017" & cor_df$preIpi=="Pretreated", "HR"])
}

cor.test(cor_matrix[,"liu_naive"], cor_matrix[,"riaz_naive"]) # r=0.87; p = 1.1e-4
cor.test(cor_matrix[,"liu_pre"], cor_matrix[,"riaz_pre"]) # r = 0.050; p = 0.87

################################
# Fig. S9: GSEA
##############################

p_ls<- readRDS("results/data/manuscript_DGE_analysis.rds")

p<- plot_grid(
  NA,p_ls$p_GSEA_preIpi,
  rel_widths = c(1,19)
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Claeys et al. 2024 - Supplementary Figure 9",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, NA, footer,
  ncol = 1,
  rel_heights = c(12,7,1)
)

ggplot2::ggsave(paste0("results/figs/manuscript_figS9_GSEA.pdf"), p, width = 178, height = 265, units = "mm")

# Some numbers
p_ls$p_GSEA_preIpi$data[
  p_ls$p_GSEA_preIpi$data$pathway%in%c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS")
  & p_ls$p_GSEA_preIpi$data$cond=="all"
  & p_ls$p_GSEA_preIpi$data$var%in%c("highTMB","strongMHC2"),c("pathway","var","padj","NES")]
# pathway        var         padj       NES
# 1: HALLMARK_INTERFERON_GAMMA_RESPONSE    highTMB 1.137001e-21  2.811693
# 2:            HALLMARK_G2M_CHECKPOINT    highTMB 4.613648e-19  2.724963
# 3:               HALLMARK_E2F_TARGETS    highTMB 1.939716e-18  2.694781
# 4: HALLMARK_INTERFERON_GAMMA_RESPONSE strongMHC2 6.417147e-30 -3.257081
# 5:            HALLMARK_G2M_CHECKPOINT strongMHC2 8.802238e-05  1.730252
# 6:               HALLMARK_E2F_TARGETS strongMHC2 2.600254e-02  1.426587

p_ls$p_GSEA_preIpi$data[
  p_ls$p_GSEA_preIpi$data$pathway%in%c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS")
  & p_ls$p_GSEA_preIpi$data$cond=="FALSE",c("pathway","var","padj","NES")]
# pathway           var         padj        NES
# 1: HALLMARK_INTERFERON_GAMMA_RESPONSE       highTMB 1.045590e-17  2.7000517
# 2:               HALLMARK_E2F_TARGETS       highTMB 1.123571e-06  2.0017115
# 3:            HALLMARK_G2M_CHECKPOINT       highTMB 9.436590e-06  1.9313404
# 4:               HALLMARK_E2F_TARGETS    strongMHC1 3.705309e-24  2.8289734
# 5:            HALLMARK_G2M_CHECKPOINT    strongMHC1 3.273479e-16  2.5187709
# 6: HALLMARK_INTERFERON_GAMMA_RESPONSE    strongMHC1 9.310148e-02 -1.2591407
# 7: HALLMARK_INTERFERON_GAMMA_RESPONSE    strongMHC2 8.700460e-30 -3.3075026
# 8:            HALLMARK_G2M_CHECKPOINT    strongMHC2 9.996884e-01  0.7103756
# 9:               HALLMARK_E2F_TARGETS    strongMHC2 9.996884e-01  0.5627072
# 10:               HALLMARK_E2F_TARGETS strongMHCnorm 2.957434e-43  3.4033589
# 11:            HALLMARK_G2M_CHECKPOINT strongMHCnorm 3.724650e-33  3.1897813
# 12: HALLMARK_INTERFERON_GAMMA_RESPONSE strongMHCnorm 8.251887e-01  0.9293197

p_ls$p_GSEA_preIpi$data[
  p_ls$p_GSEA_preIpi$data$pathway%in%c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS")
  & p_ls$p_GSEA_preIpi$data$cond=="TRUE",c("pathway","var","padj","NES")]

# pathway           var         padj       NES
# 1:            HALLMARK_G2M_CHECKPOINT       highTMB 1.363141e-18  2.615735
# 2:               HALLMARK_E2F_TARGETS       highTMB 5.491196e-11  2.265284
# 3: HALLMARK_INTERFERON_GAMMA_RESPONSE       highTMB 3.051936e-08  2.013291
# 4:               HALLMARK_E2F_TARGETS    strongMHC1 1.754815e-14 -2.518577
# 5:            HALLMARK_G2M_CHECKPOINT    strongMHC1 4.178407e-13 -2.451385
# 6: HALLMARK_INTERFERON_GAMMA_RESPONSE    strongMHC1 8.656149e-02  1.286962
# 7:            HALLMARK_G2M_CHECKPOINT    strongMHC2 6.732831e-11  2.310844
# 8:               HALLMARK_E2F_TARGETS    strongMHC2 2.400939e-10  2.277960
# 9: HALLMARK_INTERFERON_GAMMA_RESPONSE    strongMHC2 4.867595e-06 -1.881148
# 10: HALLMARK_INTERFERON_GAMMA_RESPONSE strongMHCnorm 7.038563e-05 -1.813337
# 11:            HALLMARK_G2M_CHECKPOINT strongMHCnorm 2.180549e-02 -1.436483
# 12:               HALLMARK_E2F_TARGETS strongMHCnorm 5.961633e-01 -1.013165
