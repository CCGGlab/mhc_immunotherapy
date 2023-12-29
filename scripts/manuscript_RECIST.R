###########################
# manuscript_ICB_mhc2_RECIST
###########################

library(survminer)
library(survival)
library(reshape2)

# Load data
############
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]

# MHC binders
#############
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)

# Create data frame
###################
vars<- c("highTMB","strongMHC1", "strongMHC2","strongMHCnorm")
for(v in vars) ICB_data[,v]<- as.numeric(as.logical(as.character(ICB_data[,v])))
resp_t_df<- melt(ICB_data[,c("patient","RECIST2","R","preIpi",vars)])
  
# Make high/low discrete
resp_t_df$value[resp_t_df$value==1]<- "Hi"
resp_t_df$value[resp_t_df$value==0]<- "Lo"
resp_t_df$variable<- gsub("high ", "",resp_t_df$variable)

resp_t_df<- resp_t_df[!is.na(resp_t_df$value),]
  
# Fisher test responding
p<- NULL
for(v in vars) p<- c(p,fisher.test(table(ICB_data$R, ICB_data[,v]))$p.value)
names(p)<- vars
signif(p,2)
# highTMB    strongMHC1    strongMHC2 strongMHCnorm 
# 0.1800        0.7400        0.0430        0.004 

# Chisq test
resp_t_all<- table(resp_t_df$variable,resp_t_df$RECIST2,resp_t_df$value)
chisq.test(resp_t_all["highTMB",,]) # p-value = 0.02754
chisq.test(resp_t_all["strongMHC1",,]) # p-value = 0.01041
chisq.test(resp_t_all["strongMHC2",,]) # p-value = 0.1134
chisq.test(resp_t_all["strongMHCnorm",,]) # p-value = 0.0105

# Factor for plot
resp_t_df$variable<- factor(resp_t_df$variable)
resp_t_df$value<- factor(resp_t_df$value, levels=c("Lo","Hi"))

# Plot
resp_t_df$variable_label<- NA
resp_t_df$variable_label[resp_t_df$variable=="highTMB"]<- "TMB"
resp_t_df$variable_label[resp_t_df$variable=="strongMHC1"]<- "MGBS-I"
resp_t_df$variable_label[resp_t_df$variable=="strongMHC2"]<- "MGBS-II"
resp_t_df$variable_label[resp_t_df$variable=="strongMHCnorm"]<- "MGBS-d"
resp_t_df$variable_label<- factor(resp_t_df$variable_label, levels = c("TMB","MGBS-I","MGBS-II","MGBS-d"))

p_RECIST<- ggplot(resp_t_df, aes(fill=RECIST2, y=1, x=value)) + 
  geom_bar(position="fill", stat="identity", width = 0.95) +
  # geom_text(aes(x=value, y=R, label=R, vjust=-0.2), size=6/3) +
  facet_grid(.~variable_label) +
  scale_fill_brewer(palette = "RdBu") +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(fill="RECIS2") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size=6),
    legend.key.size =  unit(4,"mm"),
    legend.title = element_text(size=7),
    text = element_text(size=6),
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    panel.spacing.x = unit(0, "lines")
  )
p_RECIST
  
p_R_ls<- list()
for(v in c("TMB","MGBS1","MGBS2","MGBSnorm")){
  resp_t_df2<- melt(ICB_data[,c("patient","R","preIpi",v)])
  resp_t_df2$R2<- NA
  resp_t_df2$R2[resp_t_df2$R==TRUE]<- "R"
  resp_t_df2$R2[resp_t_df2$R==FALSE]<- "NR"
  p <- ggboxplot(data = resp_t_df2, x = "R2", y = "value",
                 color = "R2",
                 add = "jitter", xlab = F, ylab=v, add.params = list(size=0.5))+ 
    stat_compare_means(label = "p.format", label.x=1.5, hjust=0.5, size=7*0.35)
  
  # p<- ggpar(p = p, legend = "none", font.tickslab = 6, font.y = 7)
  if(v=="TMB") p<- ggpar(p = p, yscale = "log10")
  p<- p +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7),
    )
  p_R_ls[[v]]<- p
}
p

# RECIST ~ Ipilimumab
#######################
# resp_t_df$preIpi<- resp_t_df$patient%in%ICB_data$patient[ICB_data$preIpi==T]

p_RECIST_preIpi<- ggplot(resp_t_df, aes(fill=RECIST2, y=1, x=value)) + 
  geom_bar(position="fill", stat="identity", width = 0.95) +
  # geom_text(aes(x=value, y=R, label=R, vjust=-0.2), size=6/3) +
  facet_grid(preIpi~variable_label) +
  scale_fill_brewer(palette = "RdBu") +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(fill="RECIST") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size=6),
    legend.key.size =  unit(4,"mm"),
    legend.title = element_text(size=7),
    text = element_text(size=6),
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    panel.spacing.x = unit(0, "lines")
  )
p_RECIST_preIpi

resp_t_all<- table(resp_t_df$variable,resp_t_df$RECIST2,resp_t_df$value,resp_t_df$preIpi)
chisq.test(resp_t_all["highTMB",,,"TRUE"]) # p-value = 0.29
chisq.test(resp_t_all["strongMHC1",,,"TRUE"]) # p-value = 0.57
chisq.test(resp_t_all["strongMHC2",,,"TRUE"]) # p-value = 0.071
chisq.test(resp_t_all["strongMHCnorm",,,"TRUE"]) # p-value = 0.058

chisq.test(resp_t_all["highTMB",,,"FALSE"]) # p-value = 0.048
chisq.test(resp_t_all["strongMHC1",,,"FALSE"]) # p-value = 0.031
chisq.test(resp_t_all["strongMHC2",,,"FALSE"]) # p-value = 0.42
chisq.test(resp_t_all["strongMHCnorm",,,"FALSE"]) # p-value = 0.017

# Fisher test responding
resp_t_all<- table(resp_t_df$variable,resp_t_df$R,resp_t_df$value,resp_t_df$preIpi)
fisher.test(resp_t_all["highTMB",,,"TRUE"])$p.value # p-value = 0.12
fisher.test(resp_t_all["strongMHC1",,,"TRUE"])$p.value # p-value = 0.80
fisher.test(resp_t_all["strongMHC2",,,"TRUE"])$p.value # p-value = 0.017
fisher.test(resp_t_all["strongMHCnorm",,,"TRUE"])$p.value # p-value = 0.040

fisher.test(resp_t_all["highTMB",,,"FALSE"])$p.value # p-value = 0.66
fisher.test(resp_t_all["strongMHC1",,,"FALSE"])$p.value # p-value = 0.83
fisher.test(resp_t_all["strongMHC2",,,"FALSE"])$p.value # p-value = 0.66
fisher.test(resp_t_all["strongMHCnorm",,,"FALSE"])$p.value # p-value = 0.048

# Risk groups
#############
  
# Define groups
ICB_data$riskGroup<- "Int"
ICB_data$riskGroup[ICB_data$highTMB==T&ICB_data$strongMHCnorm==T]<- "Lo"
ICB_data$riskGroup[ICB_data$highTMB==F&ICB_data$strongMHCnorm==F]<- "Hi"

# Data frame
resp_t_df<- melt(ICB_data[,c("patient","RECIST2","riskGroup","R")])
  
# Plot
p_RECIST_RG<- ggplot(resp_t_df, aes(fill=RECIST2, y=1, x=riskGroup)) + 
  geom_bar(position="fill", stat="identity", width = 0.95) +
  # geom_text(aes(x=riskGroup, y=isR, label=isR, vjust=-0.2), size=6/3) +
  scale_fill_brewer(palette = "RdBu") +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(fill="RECIST") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size=6),
    legend.key.size =  unit(4,"mm"),
    legend.title = element_text(size=7),
    text = element_text(size=6),
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    panel.spacing.x = unit(0, "lines")
  )
p_RECIST_RG
  
# P value?
chisq.test(table(resp_t_df$R,resp_t_df$riskGroup)) # p-value = 0.01574
chisq.test(table(resp_t_df$RECIST2,resp_t_df$riskGroup)) # p-value = 0.060

# Add to list
##############
p_BR_ls<- list(p_RECIST=p_RECIST, p_R_ls=p_R_ls, p_RECIST_preIpi=p_RECIST_preIpi, p_RECIST_RG=p_RECIST_RG, p=p)

# Save
######
save.image(file="results/data/manuscript_RECIST.RData")
saveRDS(p_BR_ls, file="results/data/manuscript_RECIST.rds")

