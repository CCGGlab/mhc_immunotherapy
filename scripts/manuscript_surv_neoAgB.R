###########################
# manuscript_surv_neoAgB.R
###########################

library(survminer)
library(survival)
library(reshape2)

# Load data
############
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)

# Some numbers
################
median(ICB_data$TMB) # 250.5
median(ICB_data$neoAgB1) # 46
median(ICB_data$neoAgB2) # 93.4


for(param in c("neoAgB","pMHC")){
  
  # MHC binders
  #############
  ICB_data[[paste0("high",param,"1")]]<- ICB_data[[paste0(param,"1")]] > quantile(ICB_data[[paste0(param,"1")]], 0.5, na.rm=T)
  ICB_data[[paste0("high",param,"2")]]<- ICB_data[[paste0(param,"2")]] > quantile(ICB_data[[paste0(param,"2")]], 0.5, na.rm=T)
  ICB_data[[paste0("high",param,"norm")]]<- ICB_data[[paste0(param,"norm")]] > quantile(ICB_data[[paste0(param,"norm")]], 0.5, na.rm=T)
  
  # Univariate analysis
  ############################
  vars<- c("highTMB",paste0("high",param,"1"),paste0("high",param,"2"),paste0("high",param,"norm"))

  # KM
  surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))
  for(v in vars){
    var_name<- gsub("high", "",v)
    fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", v)), data = ICB_data)
    legend_text<- paste0(c("Low", "High"),paste0(" (n=",fit$n,")"))
    p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9), title=var_name)$plot
    assign(paste0("p_", var_name), p)
  }
  
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
  
  if(param=="neoAgB") HR_df$cond<- c("high TMB", "high MHC-I NeoAgB", "high MHC-II NeoAgB", "high diff. NeoAgB")
  if(param=="pMHC") HR_df$cond<- c("high TMB", "high norm. MHC-I NeoAgB", "high norm. MHC-II NeoAgB", "high norm. diff. NeoAgB")
  HR_df$cond<- factor(HR_df$cond, levels=rev(HR_df$cond))
  fp <- ggplot(data=HR_df, aes(x=cond, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
    geom_pointrange(size=0.5, fatten = 5, shape = 18) + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
    geom_text(aes(y=3, label=paste0("P=",signif(p,2))), fontface="italic", vjust=2, hjust=1, size=6/3) +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("") +
    ylab("Hazard Ratio (Mean +/- 95% CI)") +
    scale_y_continuous(limits = c(0.1,3), trans='log10') +
    theme_bw() +  # use a white background
    theme(
      axis.text = element_text(size=7),
      axis.title = element_text(size=7)
    )
  assign(paste0("fp_", param), fp)
  
}


# Save
#########
save.image(file="results/data/manuscript_surv_neoAgB.RData")
saveRDS(list(p_TMB=p_TMB, p_pMHC1=p_pMHC1, p_pMHC2=p_pMHC2, p_pMHCnorm=p_pMHCnorm, p_neoAgB1=p_neoAgB1, p_neoAgB2=p_neoAgB2, p_neoAgBnorm=p_neoAgBnorm, fp_pMHC=fp_pMHC, fp_neoAgB=fp_neoAgB), file="results/data/manuscript_surv_neoAgB.rds")
