# Load functions
library(cowplot)
library(ggplot2)

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
  p_ls$p_fGSEA, p_right, 
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

# Survival plots MGBS2
p_surv<- plot_grid(
  p_ls$p_val_ls$Hugo_2016$OS$strongMHC2 + ggtitle("Hugo et al., 2016"),
  p_ls$p_val_ls$Gide_2019$OS$strongMHC2 + ggtitle("Gide et al., 2019"),
  p_ls$p_val_ls$Rizvi_2015$PFS$strongMHC2 + ggtitle("Rizvi et al., 2015"),
  ncol=1
)

# RECIST plots
p_RECIST<- plot_grid(
  plot_grid(NA, p_ls$p_val_ls$Gide_2019$RECIST,ncol=2,rel_widths = c(1,3)),
  p_ls$p_val_ls$Rizvi_2015$RECIST,
  ncol=1
)

# Create plot
p<- plot_grid(
  gt,
  plot_grid(p_surv, p_RECIST,ncol=3,labels=c("b","c",NA), rel_widths = c(2,3,1), label_size = 8),
  ncol = 1,
  rel_heights = c(2,3),
  labels = c("a",NA),
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

# 
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
#   glm(formula = R ~ MGBS2 + MHC2 + hasLOHB2M + heterogeneity, family = binomial, 
#       data = ICB_data[ICB_data$preIpi == isPreIpi, ])
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -2.68395  -0.32969   0.02498   0.34800   1.81980  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -6.416e-04  3.061e+00   0.000  0.99983   
# MGBS2         -6.531e+01  2.659e+01  -2.456  0.01404 * 
#   MHC2           9.524e-04  3.655e-04   2.606  0.00917 **
#   hasLOHB2M     -8.341e+00  4.195e+00  -1.988  0.04680 * 
#   heterogeneity  1.142e+01  5.191e+00   2.199  0.02786 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 56.814  on 40  degrees of freedom
# Residual deviance: 23.973  on 36  degrees of freedom
# AIC: 33.973
# 
# Number of Fisher Scoring iterations: 7

################################
# Fig. S8: GSEA
##############################

p_ls<- readRDS("results/data/manuscript_DGE_analysis.rds")

p<- plot_grid(
  NA,p_ls$p_GSEA_preIpi,
  rel_widths = c(1,19)
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
  p, NA, footer,
  ncol = 1,
  rel_heights = c(12,7,1)
)

ggplot2::ggsave(paste0("results/figs/manuscript_figS8_GSEA.pdf"), p, width = 178, height = 265, units = "mm")

