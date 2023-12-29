# Load data
############
library(tidyverse)
library(ggrepel)
library(fgsea)
load("data/MHC_immunotherapy.RData")
DGE_ls<- readRDS("temp/DGE_ls.rds")
  
# plot ipi versus no ipi
#########################

DGE_Ha_df<- NULL
for(v in names(DGE_ls)){
  for(preIpi in c("all","TRUE","FALSE")){
    DGE_Ha_df_tmp<- DGE_ls[[v]][[preIpi]]$GSEA
    DGE_Ha_df_tmp$var<- v
    DGE_Ha_df_tmp$cond<- preIpi
    DGE_Ha_df<- rbind(DGE_Ha_df,DGE_Ha_df_tmp)
  }
}

# Sort based on overal MGBS2
DGE_Ha_df_tmp<- DGE_ls$strongMHC2$all$GSEA
pws<- DGE_Ha_df_tmp$pathway[order(DGE_Ha_df_tmp$pval)]
DGE_Ha_df$pathway<- factor(DGE_Ha_df$pathway, pws)

# Name variable
DGE_Ha_df$var_name<- "TMB"
DGE_Ha_df$var_name[DGE_Ha_df$var=="strongMHC1"]<- "MGBS-I"
DGE_Ha_df$var_name[DGE_Ha_df$var=="strongMHC2"]<- "MGBS-II"
DGE_Ha_df$var_name[DGE_Ha_df$var=="strongMHCnorm"]<- "MGBS-d"
DGE_Ha_df$var_name<- factor(DGE_Ha_df$var_name, levels=unique(DGE_Ha_df$var_name))

# Name condition
DGE_Ha_df$cond_name<- "All patients (n=144)"
DGE_Ha_df$cond_name[DGE_Ha_df$cond=="TRUE"]<- "Ipilimumab pretreated patients (n=60)"
DGE_Ha_df$cond_name[DGE_Ha_df$cond=="FALSE"]<- "Treatment naive patients (n=84)"
DGE_Ha_df$cond_name<- factor(DGE_Ha_df$cond_name, levels=c("All patients (n=144)", "Treatment naive patients (n=84)", "Ipilimumab pretreated patients (n=60)"))

p_GSEA_preIpi<- ggplot(DGE_Ha_df, aes(x = pathway, y = var_name, colour = NES, size = -log10(pval))) +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust=1, vjust=1, size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # legend.position = "bottom",
    legend.key.size = unit(8, 'points'),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.position = "top",
    axis.text.y =  element_text(size = 7),
  ) +
  scale_colour_steps2(low = scales::muted("blue"), mid = "white", high = scales::muted("red")) +
  facet_wrap(~cond_name, ncol = 1) +
  scale_y_discrete(limits=rev) +
  labs(x = "pathway", y = "metric", colour = "NES")
# p_GSEA_preIpi

# Plot different TMB Vs MGBS2
DGE_Ha_df_sel<- merge(DGE_Ha_df[DGE_Ha_df$var=="highTMB"&DGE_Ha_df$cond=="all",],DGE_Ha_df[DGE_Ha_df$var=="strongMHC2"&DGE_Ha_df$cond=="all",],by="pathway")
DGE_Ha_df_sel$padj_mean<- rowMeans(DGE_Ha_df_sel[,c("padj.x", "padj.y")])
fGSEA_plot<- ggplot(DGE_Ha_df_sel, aes(x = -log10(padj.x), y = -log10(padj.y), fill = NES.y, colour = NES.x, size= -log10(padj_mean))) +
  geom_point(shape=21, stroke=3) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(6, 'points'),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7)  
    ) +
  scale_colour_steps2(low = scales::muted("blue"), mid = "white", high = scales::muted("red")) +
  scale_fill_steps2(low = scales::muted("blue"), mid = "white", high = scales::muted("red")) +
  labs(x = "-log10(Padj) TMB", y = "-log10(Padj) TMB MGBS-II", colour = "NES TMB", fill = "NES MGBS-II", size = "-log10(Padj)") +
  geom_label_repel(data = subset(DGE_Ha_df_sel, padj.x<1e-4|padj.y<1e-4),
                   aes(label = pathway, fontface=3),
                   size = 2,
                   colour="black",
                   box.padding   = 0.5,
                   point.padding = 0.5,
                   max.overlaps = 25,
                   segment.color = 'grey50',
                   fill=NA,
                   label.size=NA)

# fGSEA_plot

# Running score plots
#####################
p_RS_ls<- list()
pw<- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
for(v in c("highTMB","strongMHC2")){
  
  res_tmp<- DGE_ls[[v]][["all"]]
  p_pw<- signif(res_tmp$GSEA$padj[res_tmp$GSEA$pathway==pw],3)
  ES_pw<- res_tmp$GSEA$ES[res_tmp$GSEA$pathway==pw]
  plot_lims<- c(-0.1,0.6) 
  if(ES_pw<0) plot_lims<- c(-0.6,0.1) 
  
  if(v=="highTMB") varname<- "TMB"
  if(v=="strongMHC2") varname<- "MGBS-II"
  
  stat<- res_tmp$DGE$t
  names(stat)<- rownames(res_tmp$DGE)
  p_RS_tmp<- plotEnrichment(Ha_ls[[pw]],stat,ticksSize=.1,ticksLength = 0.2) + ggtitle("TMB - ALL") +
    ggtitle(paste0(pw,"\n",varname,"(Padj=",p_pw,")")) + 
    geom_line(size=0.5, col="green") +
    theme(
      plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
    ) +
    scale_x_continuous(name = "Rank") +
    scale_y_continuous(name = "Enrichment Score", limits = plot_lims)
  p_RS_ls[[v]]<- p_RS_tmp
}
p_RS_ls$highTMB
p_RS_ls$strongMHC2

# Save
#########
p_DGE_ls<- list(p_fGSEA=fGSEA_plot, p_RS_ls=p_RS_ls, p_GSEA_preIpi=p_GSEA_preIpi)
save.image(file="results/data/manuscript_DGE_analysis.RData")
saveRDS(p_DGE_ls, file="results/data/manuscript_DGE_analysis.rds")

