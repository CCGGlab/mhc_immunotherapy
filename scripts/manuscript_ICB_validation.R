###########################
# manuscript_ICB_mhc2_surv
###########################

library(survminer)
library(survival)
library(reshape2)

# Load data
############
# source("scripts/functions/plot_CPH_forest.R")
load("data/MHC_immunotherapy.RData")

# Hugo et al; Melanoma study
# Rizvi et al; NSCLC study
##############################
p_val_ls<- list()
HR_df_all<- NULL

studies<- c("Liu_2019", "Hugo_2016", "Gide_2019", "Riaz_2017", "Rizvi_2015", "Snyder_2014","VanAllen_2015") 
for(s in studies){
  
  for(hasPreIpi in c("all", T, F)){
    
    ICB_data<- ICB_study[ICB_study$study==s,]
    
    # MHC binders
    ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
    ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
    ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
    ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)
    
    if(hasPreIpi!="all"&!s%in%c("Liu_2019","Riaz_2017")) next
    if(hasPreIpi!="all") ICB_data<- ICB_data[ICB_data$preIpi==hasPreIpi,]
    
    # Risk groups
    ICB_data$riskGroup<- "Intermediary"
    ICB_data$riskGroup[ICB_data$highTMB==F&ICB_data$strongMHCnorm==F]<- "High"
    ICB_data$riskGroup[ICB_data$highTMB==T&ICB_data$strongMHCnorm==T]<- "Low"
    ICB_data$riskGroup<- factor(ICB_data$riskGroup, levels=c("Low","Intermediary","High"))
    
    # Sruvival
    vars<- c("highTMB","strongMHC1", "strongMHC2","strongMHCnorm")
    surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))
    
    # KM: observed
    if(!is.na(ICB_data$dead[1])){
      for(v in vars){
        if(is.na(ICB_data[[v]][1])) next
        var_name<- gsub("high", "",v)
        fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", v)), data = ICB_data)
        legend_text<- paste0(c("Low", "High"),paste0(" (n=",fit$n,")"))
        p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9), title=var_name)$plot
        p_val_ls[[s]][["OS"]][[v]][[hasPreIpi]]<- p
        
      if(v=="highTMB"){
        fit <- survfit(Surv(time = overall_survival, event = dead) ~ riskGroup, data = ICB_data) 
        legend_text<- paste0(levels(ICB_data$riskGroup),paste0(" (n=",fit$n,")"))
        p_RG<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="Risk group", font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8))$plot
        p_val_ls[[s]][["OS"]][["RG"]][[hasPreIpi]]<- p_RG
      }
      
      }
    }
    # Median surv
    # highTMB 2.291581 2.683094 
    # strongMHC1 2.291581 2.666667 
    # strongMHC2 2.666667 1.891855 
    # strongMHCnorm 1.891855 2.683094 
  
    # KM: progression-free
    if(!is.na(ICB_data$progression[1])){
      for(v in vars){
        if(is.na(ICB_data[[v]][1])) next
        var_name<- gsub("high", "",v)
        fit <- survfit(as.formula(paste0("Surv(time = progression_free, event = progression) ~ ", v)), data = ICB_data)
        legend_text<- paste0(c("Low", "High"),paste0(" (n=",fit$n,")"))
        p<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9), title=var_name)$plot
        p_val_ls[[s]][["PFS"]][[v]]<- p
        
        # RG
        if(v=="highTMB"){
          fit <- survfit(Surv(time = progression_free, event = progression) ~ riskGroup, data = ICB_data) 
          legend_text<- paste0(levels(ICB_data$riskGroup),paste0(" (n=",fit$n,")"))
          p_RG<- ggsurvplot(fit, data = ICB_data, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="Risk group", font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.8))$plot
          p_val_ls[[s]][["PFS"]][["RG"]][[hasPreIpi]]<- p_RG
        }
      }
    }
    
    # Hazard ratios
    HR_df<- data.frame(
      cond=vars,
      HR=NA,
      ci_l=NA,
      ci_h=NA,
      p=NA,
      n1=NA,
      n=NA,
      preIpi=hasPreIpi,
      study=s
    )
    
    for(i in 1:nrow(HR_df)){
      if(s=="Gide_2019"&HR_df$cond[i]=="highTMB") next # No MR
      ICB_data$cond<- ICB_data[,HR_df$cond[i]]
      if(s!="Rizvi_2015") fit.coxph <- coxph(Surv(time = overall_survival, event = dead) ~ cond, data = ICB_data)
      else fit.coxph <- coxph(Surv(time = progression_free, event = progression) ~ cond, data = ICB_data) # Only PFS data
      HR_df[i,c("HR","ci_l","ci_h")]<- summary(fit.coxph)$conf.int["condTRUE", c(1,3,4)]
      HR_df[i,"p"]<- summary(fit.coxph)$coefficients["condTRUE",5]
      HR_df[i,c("n1","n")]<- c(summary(fit.coxph)$nevent, summary(fit.coxph)$n)
    }
    HR_df$cond_name<- c("high TMB", "strong MHC1", "strong MHC2", "strong MHCnorm")
    HR_df_all<- rbind(HR_df_all, HR_df)
  
    # Create data frame
    for(v in vars) ICB_data[,v]<- as.numeric(as.logical(as.character(ICB_data[,v])))
    resp_t_df<- melt(ICB_data[,c("patient","RECIST2",vars)])
    
    # Make high/low discrete
    resp_t_df$value[resp_t_df$value==1]<- "Hi"
    resp_t_df$value[resp_t_df$value==0]<- "Lo"
    resp_t_df$variable<- gsub("high ", "",resp_t_df$variable)
    resp_t_df<- resp_t_df[!is.na(resp_t_df$value),]
    
    # Chisq test
    if(is.na(ICB_data$RECIST2[1])) next
    resp_t_all<- table(resp_t_df$variable,resp_t_df$RECIST2,resp_t_df$value)
    # resp_t_all<- resp_t_all[,c("PD","MR/SD","PR"),] # Rizvi, no CR, Chisq returns NA
    # resp_t_all<- resp_t_all[,c("PD","PR","CR"),] # Hugo, no MR/SD, Chisq returns NA
    # chisq.test(resp_t_all["highTMB",,]) # p-value = NA (Gide); 0.013 (Rizvi); 0.44 (Hugo)
    chisq.test(resp_t_all["strongMHC1",,]) # p-value = 0.77; 0.72; 0.40
    chisq.test(resp_t_all["strongMHC2",,]) # p-value = 0.0076; 0.23; 0.90
    chisq.test(resp_t_all["strongMHCnorm",,]) # p-value = 0.027; 0.14; 0.25
    
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
      # geom_text(aes(x=value, y=isR, label=isR, vjust=-0.2), size=6/3) +
      facet_grid(.~variable_label) +
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
    p_val_ls[[s]][["RECIST"]][[hasPreIpi]]<- p_RECIST
  }
}

# Create overall FP
HR_df_all$cond_name<- factor(HR_df_all$cond_name, levels=c("high TMB", "strong MHC1", "strong MHC2", "strong MHCnorm"))
HR_df_all$study<- factor(HR_df_all$study, levels=rev(studies))
HR_df_all$treatment<- "anti-PD1"
HR_df_all$treatment[HR_df_all$study%in%c("Snyder_2014","VanAllen_2015")]<- "anti-CTLA4"
HR_df_all$treatment<- factor(HR_df_all$treatment, levels=c("anti-PD1","anti-CTLA4"))
HR_df_all$ci_l[!is.na(HR_df_all$ci_l)&HR_df_all$ci_l<0.1]<- 0.1 # FOr visualization, otherwise CI not drawn
fp <- ggplot(data=subset(HR_df_all,preIpi=="all"), aes(x=study, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.5, fatten = 5, shape = 18) + 
  # geom_pointrange(size=.3, fatten = 5) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
  geom_text(aes(y=10, label=paste0("P=",signif(p,2))), fontface="italic", hjust=1, vjust=2, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_grid(treatment~cond_name, scales = "free")+
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.1,10), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
p_val_ls$fp<- fp

# Create preIpi FP
fp_preIpi <- ggplot(data=subset(HR_df_all,preIpi!="all"), aes(x=study, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.5, fatten = 5, shape = 18) + 
  # geom_pointrange(size=.3, fatten = 5) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
  geom_text(aes(y=10, label=paste0("P=",signif(p,2))), fontface="italic", hjust=1, vjust=2, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_grid(preIpi~cond_name, scales = "free")+
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.1,10), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
p_val_ls$fp_preIpi<- fp_preIpi

# Save
#########
save.image(file="results/data/manuscript_validation.RData")
saveRDS(list(p_val_ls=p_val_ls, HR_df_all=HR_df_all), file="results/data/manuscript_validation.rds")



