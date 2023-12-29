###########################
# manuscript_ICB_mhc2_surv
###########################

library(survminer)
library(survival)
library(reshape2)

# Load data
############
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)

# Some number
#############
quantile(ICB_data$TMB, 0.5, na.rm=T) # 250.5
quantile(ICB_data$MGBS1, 0.5, na.rm=T) # 0.01594183 
quantile(ICB_data$MGBS2, 0.5, na.rm=T) # 0.1176544
quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)  # 0.009604019 

# MHC binders
#############
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)

# Univariate analysis
############################
vars<- c("highTMB","strongMHC1", "strongMHC2","strongMHCnorm")
surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))

# KM
for(v in vars){
  var_name<- gsub("high", "",v)
  fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", v)), data = ICB_data)
  legend_text<- paste0(c("Low", "High"),paste0(" (n=",fit$n,")"))
  p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9), title=var_name)$plot
  assign(paste0("p_", var_name), p)
  cat(v,surv_median(fit)[,"median"],"\n")
}

# Median surv
# highTMB 1.360712 NA 
# strongMHC1 1.713895 2.784394 
# strongMHC2 NA 1.566051 
# strongMHCnorm 1.374401  NA 

p_TMB
p_strongMHC1
p_strongMHC2
p_strongMHCnorm

# EF
for(v in vars){
  var_name<- gsub("high", "",v)
  fit <- survfit(as.formula(paste0("Surv(time = progression_free, event = progression) ~ ", v)), data = ICB_data)
  legend_text<- paste0(c("Low", "High"),paste0(" (n=",fit$n,")"))
  p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9), title=var_name)$plot
  cat(v,surv_median(fit)[,"median"],"\n")
  assign(paste0("p_EF_", var_name), p)
}

# highTMB 0.2642026 0.8199863 
# strongMHC1 0.3134839 0.5010267 
# strongMHC2 0.6570842 0.2642026 
# strongMHCnorm 0.2491444 0.7665982 

p_EF_TMB
p_EF_strongMHC1
p_EF_strongMHC2
p_EF_strongMHCnorm # P=0.01

# univariate hazard ratios
HR_df<- data.frame(
  cond=vars,
  HR=NA,
  ci_l=NA,
  ci_h=NA,
  p=NA,
  n1=NA,
  n=NA
)

for(i in 1:nrow(HR_df)){
  ICB_data$cond<- ICB_data[,HR_df$cond[i]]
  fit.coxph <- coxph(Surv(time = overall_survival, event = dead) ~ cond, data = ICB_data)
  HR_df[i,c("HR","ci_l","ci_h")]<- summary(fit.coxph)$conf.int["condTRUE", c(1,3,4)]
  HR_df[i,"p"]<- summary(fit.coxph)$coefficients["condTRUE",5]
  HR_df[i,c("n1","n")]<- c(summary(fit.coxph)$nevent, summary(fit.coxph)$n)
}

HR_df$cond<- c("high TMB", "strong MHC1", "strong MHC2", "strong MHCnorm")
HR_df$cond<- factor(HR_df$cond, levels=rev(HR_df$cond))
fp <- ggplot(data=HR_df, aes(x=cond, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.5, fatten = 10, shape = 18) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-2, size=7/3) +
  geom_text(aes(y=3, label=paste0("P=",signif(p,2))), fontface="italic", vjust=2, hjust=1, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.1,10), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
fp

#  Multivariate analysis with MB
################################
fp_multi_MB_1<- plot_CPH_forest(ICB_data, c("highTMB", "strongMHC1")) # P=0.0017 P=0.065
fp_multi_MB_2<- plot_CPH_forest(ICB_data, c("highTMB", "strongMHC2")) # P=0.0069 P=0.032
fp_multi_12<- plot_CPH_forest(ICB_data, c("strongMHC1", "strongMHC2")) # P=0.15 P=0.015
fp_multi_MB_12<- plot_CPH_forest(ICB_data, c("highTMB", "strongMHC1", "strongMHC2")) # P=0.0065; P=0.16; P=0.059  
fp_multi_MB_norm<- plot_CPH_forest(ICB_data, c("highTMB", "strongMHCnorm")) # P=0.015; P=0.0054;  
fp_multi_MB_12norm<- plot_CPH_forest(ICB_data, c("highTMB", "strongMHC1", "strongMHC2","strongMHCnorm")) # P=0.016; P=0.6; P=0.46; P=0.26  

# Use different thresholds
##########################

HR_df_all<- NULL
p_th_ls<- list()
vars<- c("TMBclass","MHC1class","MHC2class","MHCnormclass")
for(th in c(0.5, 0.75, 0.90,0.95)){
  ICB_data$TMBclass<- NA
  ICB_data$TMBclass[ICB_data$TMB < quantile(ICB_data$TMB, 1-th, na.rm=T)]<- "Lo" 
  ICB_data$TMBclass[ICB_data$TMB >= quantile(ICB_data$TMB, 1-th, na.rm=T)]<- "Im" 
  ICB_data$TMBclass[ICB_data$TMB >= quantile(ICB_data$TMB, th, na.rm=T)]<- "Hi" 
  
  ICB_data$MHC1class<- NA
  ICB_data$MHC1class[ICB_data$MGBS1 < quantile(ICB_data$MGBS1, 1-th, na.rm=T)]<- "Lo" 
  ICB_data$MHC1class[ICB_data$MGBS1 >= quantile(ICB_data$MGBS1, 1-th, na.rm=T)]<- "Im" 
  ICB_data$MHC1class[ICB_data$MGBS1 >= quantile(ICB_data$MGBS1, th, na.rm=T)]<- "Hi" 
  
  ICB_data$MHC2class<- NA
  ICB_data$MHC2class[ICB_data$MGBS2 < quantile(ICB_data$MGBS2, 1-th, na.rm=T)]<- "Lo" 
  ICB_data$MHC2class[ICB_data$MGBS2 >= quantile(ICB_data$MGBS2, 1-th, na.rm=T)]<- "Im" 
  ICB_data$MHC2class[ICB_data$MGBS2 >= quantile(ICB_data$MGBS2, th, na.rm=T)]<- "Hi" 
  
  ICB_data$MHCnormclass<- NA
  ICB_data$MHCnormclass[ICB_data$MGBSnorm < quantile(ICB_data$MGBSnorm, 1-th, na.rm=T)]<- "Lo" 
  ICB_data$MHCnormclass[ICB_data$MGBSnorm >= quantile(ICB_data$MGBSnorm, 1-th, na.rm=T)]<- "Im" 
  ICB_data$MHCnormclass[ICB_data$MGBSnorm >= quantile(ICB_data$MGBSnorm, th, na.rm=T)]<- "Hi" 
  
  # Plot survival
  surv_cols<- c(rgb(231/256, 184/256, 0/256),"black",rgb(46/256, 159/256, 223/256))
  if(th!=0.5){
    for(v in vars){
      fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", v)), data = ICB_data)
      legend_text<- paste0(c(paste0(">P",100*th), paste0("P",100*(1-th),"-","P",100*(th)), paste0("<P",100*(1-th))),paste0(" (n=",fit$n,")"))
      p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=6, font.title=8, legend=c(0.8,0.9), title=v)$plot
      p_th_ls[[v]][[as.character(th)]]<- p
    }
  }
  
  # univariate hazard ratios
  HR_df<- data.frame(
    cond=rep(vars,each=2),
    class=rep(c("Hi","Lo"),length(vars)),
    th=as.character(th),
    HR=NA,
    ci_l=NA,
    ci_h=NA,
    p=NA,
    n1=NA,
    n=NA
  )
  
  for(i in 1:nrow(HR_df)){
    ICB_data$cond<- ICB_data[,HR_df[i,"cond"]]
    ICB_data$cond<- factor(ICB_data$cond, levels=c("Im","Hi","Lo"))
    fit.coxph <- coxph(Surv(time = overall_survival, event = dead) ~ cond, data = ICB_data)
    HR_df[i,c("HR","ci_l","ci_h")]<- summary(fit.coxph)$conf.int[paste0("cond", HR_df[i,"class"]), c(1,3,4)]
    HR_df[i,"p"]<- summary(fit.coxph)$coefficients[paste0("cond",HR_df[i,"class"]), 5]
    HR_df[i,c("n1","n")]<- c(summary(fit.coxph)$nevent, summary(fit.coxph)$n)
  }
  HR_df_all<- rbind(HR_df_all,HR_df)
}  
# View(HR_df_all)

# Plot

# Change labelname threshold
HR_df_all$th_label<- paste0(100*(1-as.numeric(HR_df_all$th)),"%")
HR_df_all$th_label<- factor(HR_df_all$th_label, levels=c("5%","10%","25%","50%"))

# Change labelname class
HR_df_all$class_label<- NA
HR_df_all$class_label[HR_df_all$class=="Hi"]<- "Strong"
HR_df_all$class_label[HR_df_all$class=="Lo"]<- "Weak"

# Change labelname condition
HR_df_all$cond_label<- NA
HR_df_all$cond_label[HR_df_all$cond=="MHC1class"]<- "MGBS-I"
HR_df_all$cond_label[HR_df_all$cond=="MHC2class"]<- "MGBS-II"

# Restrict CI to legend (log) scales
HR_df_all$ci_l[HR_df_all$ci_l<0.1]<- 0.1

# plot
surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))
fp_th<- ggplot(data=subset(HR_df_all,cond%in%c("MHC1class","MHC2class")), aes(x=th_label, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)),col=class_label)) +
  # geom_pointrange(size=.5, fatten = 5) +
  geom_pointrange(size=0.5, fatten = 10, shape = 18) + 
  # scale_color_brewer(palette="Dark2") +
  scale_color_manual(values = rev(surv_cols)) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  # geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_wrap(cond_label~., nrow = 2) +
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.1,10), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    axis.text = element_text(size=7),
    axis.title = element_text(size=7),
    legend.position = "bottom",
    legend.title = element_blank()  
  )
fp_th

# Survival analysis on Risk Groups
##################################
ICB_data$riskGroup<- "Intermediary"
ICB_data$riskGroup[ICB_data$highTMB==F&ICB_data$strongMHCnorm==F]<- "High"
ICB_data$riskGroup[ICB_data$highTMB==T&ICB_data$strongMHCnorm==T]<- "Low"
ICB_data$riskGroup<- factor(ICB_data$riskGroup, levels=c("Low","Intermediary","High"))

# Observed
fit <- survfit(Surv(time = overall_survival, event = dead) ~ riskGroup, data = ICB_data) 
legend_text<- paste0(levels(ICB_data$riskGroup),paste0(" (n=",fit$n,")"))
p_RG_obs<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="Risk group", font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8))$plot
p_RG_obs
surv_median(fit)

# strata   median     lower    upper
# 1          riskGroup=Low       NA        NA       NA
# 2 riskGroup=Intermediary 1.878166 1.4264203       NA
# 3         riskGroup=High 1.077344 0.7118412 2.568104

# Progression free
fit <- survfit(Surv(time = progression_free, event = progression) ~ riskGroup, data = ICB_data) 
legend_text<- paste0(levels(ICB_data$riskGroup),paste0(" (n=",fit$n,")"))
p_RG_EF<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="Risk group", font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8))$plot
p_RG_EF
surv_median(fit)
# strata    median     lower     upper
# 1          riskGroup=Low 1.0650240 0.5201916        NA
# 2 riskGroup=Intermediary 0.3641342 0.2354552 1.1115674
# 3         riskGroup=High 0.2532512 0.2299795 0.5749487

# HRs?
ICB_data$riskGroup<- factor(ICB_data$riskGroup, levels=c("Intermediary","Low","High"))
fit.coxph <- coxph(Surv(time = overall_survival, event = dead) ~ riskGroup, data = ICB_data)
summary(fit.coxph)$conf.int
# exp(coef) exp(-coef) lower .95 upper .95
# riskGroupLow  0.3816499  2.6202022 0.1911448 0.7620228
# riskGroupHigh 1.5305939  0.6533412 0.9253718 2.5316501
summary(fit.coxph)$coefficients
# coef exp(coef)  se(coef)         z    Pr(>|z|)
# riskGroupLow  -0.9632515 0.3816499 0.3527987 -2.730315 0.006327388
# riskGroupHigh  0.4256558 1.5305939 0.2567473  1.657878 0.097342010

# Multivariate analysis with ipilimumab
########################################

vars<- c("highTMB","strongMHC1", "strongMHC2", "strongMHCnorm")

# KM all
fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ preIpi")), data = ICB_data)
legend_text<- paste0(c("No","Yes"), " (n=",fit$n,")")
p_ipi_vs_noIpi<- ggsurvplot(fit, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", linetype = c("solid","dashed"), palette =  c("black","black"), font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8), censor=F)$plot + theme(legend.background = element_blank())

# KM stratified preIpi
p_ipi_ls<- list()
p_noIpi_ls<- list()
for(v in vars){
  surv_cols<- c(rgb(231/256, 184/256, 0/256,.7),rgb(46/256, 159/256, 223/256,.7))
  
  fit_noPreIpi <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ preIpi + ", v)), data = subset(ICB_data,preIpi==F))
  legend_text<- paste0(c("Low", "High"), " (n=",fit_noPreIpi$n,")")
  p_noIpi_ls[[v]]<- ggsurvplot(fit_noPreIpi, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", linetype = "solid", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8), censor=F)$plot + theme(legend.background = element_blank())
  
  fit_preIpi <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ preIpi + ", v)), data = subset(ICB_data,preIpi==T))
  legend_text<- paste0(c("Low", "High"), " (n=",fit_preIpi$n,")")
  p_ipi_ls[[v]]<- ggsurvplot(fit_preIpi, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", linetype = "dashed", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8), censor=F)$plot + theme(legend.background = element_blank())
}    
p_ipi_ls$highTMB
p_ipi_ls$strongMHC1
p_ipi_ls$strongMHC2
p_ipi_ls$strongMHCnorm

p_noIpi_ls$highTMB
p_noIpi_ls$strongMHC1
p_noIpi_ls$strongMHC2
p_noIpi_ls$strongMHCnorm

# Cox HR
HR_multi_df<- data.frame(
  cond=rep(vars,2),
  HR=NA,
  ci_l=NA,
  ci_h=NA,
  p=NA,
  n1=NA,
  n=NA,
  preIpi=rep(c(T,F),each=length(vars))
)
HR_multi_df$ipiLabs<- "Ipilimumab pretreated"
HR_multi_df$ipiLabs[HR_multi_df$preIpi==F]<- "Treatment naive"
HR_multi_df$ipiLabs<- factor(HR_multi_df$ipiLabs, levels=c("Treatment naive","Ipilimumab pretreated"))

HR_multi_df_1<- HR_multi_df[HR_multi_df$cond!="strongMHC2"&HR_multi_df$cond!="strongMHCnorm",]
HR_multi_df_2<- HR_multi_df[HR_multi_df$cond!="strongMHC1"&HR_multi_df$cond!="strongMHCnorm",]
HR_multi_df_n<- HR_multi_df[HR_multi_df$cond!="strongMHC1"&HR_multi_df$cond!="strongMHC2",]

for(isPreIpi in c(T,F)){
  fit.coxph <- coxph(Surv(time = overall_survival, event = dead) ~ highTMB + strongMHC1 + strongMHC2 + strongMHCnorm, data = subset(ICB_data,preIpi==isPreIpi))
  for(var in vars){
    HR_multi_df[HR_multi_df$cond==var&HR_multi_df$preIpi==isPreIpi,c("HR","ci_l","ci_h")]<- summary(fit.coxph)$conf.int[paste0(var,"TRUE"), c(1,3,4)]
    HR_multi_df[HR_multi_df$cond==var&HR_multi_df$preIpi==isPreIpi,"p"]<- summary(fit.coxph)$coefficients[paste0(var,"TRUE"),5]
    HR_multi_df[HR_multi_df$cond==var&HR_multi_df$preIpi==isPreIpi,c("n1","n")]<- c(summary(fit.coxph)$nevent, summary(fit.coxph)$n)
  }
  fit.coxph_1 <- coxph(Surv(time = overall_survival, event = dead) ~ highTMB + strongMHC1, data = subset(ICB_data,preIpi==isPreIpi))
  for(var in c("highTMB","strongMHC1")){
    HR_multi_df_1[HR_multi_df_1$cond==var&HR_multi_df_1$preIpi==isPreIpi,c("HR","ci_l","ci_h")]<- summary(fit.coxph_1)$conf.int[paste0(var,"TRUE"), c(1,3,4)]
    HR_multi_df_1[HR_multi_df_1$cond==var&HR_multi_df_1$preIpi==isPreIpi,"p"]<- summary(fit.coxph_1)$coefficients[paste0(var,"TRUE"),5]
    HR_multi_df_1[HR_multi_df_1$cond==var&HR_multi_df_1$preIpi==isPreIpi,c("n1","n")]<- c(summary(fit.coxph_1)$nevent, summary(fit.coxph)$n)
  }
  fit.coxph_2 <- coxph(Surv(time = overall_survival, event = dead) ~ highTMB + strongMHC2, data = subset(ICB_data,preIpi==isPreIpi))
  for(var in c("highTMB","strongMHC2")){
    HR_multi_df_2[HR_multi_df_2$cond==var&HR_multi_df_2$preIpi==isPreIpi,c("HR","ci_l","ci_h")]<- summary(fit.coxph_2)$conf.int[paste0(var,"TRUE"), c(1,3,4)]
    HR_multi_df_2[HR_multi_df_2$cond==var&HR_multi_df_2$preIpi==isPreIpi,"p"]<- summary(fit.coxph_2)$coefficients[paste0(var,"TRUE"),5]
    HR_multi_df_2[HR_multi_df_2$cond==var&HR_multi_df_2$preIpi==isPreIpi,c("n1","n")]<- c(summary(fit.coxph_2)$nevent, summary(fit.coxph)$n)
  }
  fit.coxph_n <- coxph(Surv(time = overall_survival, event = dead) ~ highTMB + strongMHCnorm, data = subset(ICB_data,preIpi==isPreIpi))
  for(var in c("highTMB","strongMHCnorm")){
    HR_multi_df_n[HR_multi_df_n$cond==var&HR_multi_df_n$preIpi==isPreIpi,c("HR","ci_l","ci_h")]<- summary(fit.coxph_n)$conf.int[paste0(var,"TRUE"), c(1,3,4)]
    HR_multi_df_n[HR_multi_df_n$cond==var&HR_multi_df_n$preIpi==isPreIpi,"p"]<- summary(fit.coxph_n)$coefficients[paste0(var,"TRUE"),5]
    HR_multi_df_n[HR_multi_df_n$cond==var&HR_multi_df_n$preIpi==isPreIpi,c("n1","n")]<- c(summary(fit.coxph_n)$nevent, summary(fit.coxph)$n)
  }
}

for(i in 0:3){
  HR_multi_df_tmp<- list(HR_multi_df,HR_multi_df_1,HR_multi_df_2, HR_multi_df_n)[[i+1]]
  HR_multi_df_tmp$cond<- factor(HR_multi_df_tmp$cond, levels=rev(unique(HR_multi_df_tmp$cond)))
  fp_multi <- ggplot(data=HR_multi_df_tmp, aes(x=cond, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
    # geom_pointrange(size=.3, fatten = 5) + 
    geom_pointrange(size=0.5, fatten = 10, shape = 18) + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
    geom_text(aes(y=6.5, label=paste0("P=",signif(p,2))), fontface="italic", hjust=0, size=6/3) +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    facet_grid(.~ipiLabs) +  
    xlab("") +
    ylab("Hazard Ratio (Mean +/- 95% CI)") +
    scale_y_continuous(limits = c(0.15,30), trans='log10') +
    theme_bw() +  # use a white background
    theme(
      axis.text = element_text(size=7),
      axis.title = element_text(size=7)
    )
  assign(paste0("fp_multi_",i),fp_multi)
}

# Plots to list 
p_surv_ls<- list(
  p_TMB=p_TMB, p_strongMHC1=p_strongMHC1, p_strongMHC2=p_strongMHC2, p_strongMHCnorm=p_strongMHCnorm,
  p_EF_TMB=p_EF_TMB, p_EF_strongMHC1=p_EF_strongMHC1, p_EF_strongMHC2=p_EF_strongMHC2, p_EF_strongMHCnorm=p_EF_strongMHCnorm,
  p_ipi_ls=p_ipi_ls, p_noIpi_ls=p_noIpi_ls, p_ipi_vs_noIpi=p_ipi_vs_noIpi,
  HR_df=HR_df_all, HR_multi_df=HR_multi_df, 
  fp_multi_MB_1=fp_multi_MB_1, fp_multi_MB_2=fp_multi_MB_2, fp_multi_MB_12=fp_multi_MB_12,fp_multi_MB_norm=fp_multi_MB_norm,  
  fp_th=fp_th,
  p_RG_obs=p_RG_obs, p_RG_EF=p_RG_EF,
  fp_ipi_multi=fp_multi_0, fp_ipi_multi_1=fp_multi_1, fp_ipi_multi_2=fp_multi_2, fp_ipi_multi_n=fp_multi_3,
  p_th_ls=p_th_ls
  )

# Save
#########
save.image(file="results/data/manuscript_surv.RData")
saveRDS(p_surv_ls, file="results/data/manuscript_surv.rds")
# load(file="results/data/manuscript_surv.RData")
