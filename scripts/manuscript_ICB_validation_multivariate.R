# Use the Liu logistic regression model (ipi pretreated) to independetly predict response in other datasets
###########################################################################################################

# Load data
load("data/MHC_immunotherapy.RData")
p_ls<- list()

# Get model from Liu - preIpi
ICB_data<- ICB_study[!is.na(ICB_study$study)&ICB_study$study=="Liu_2019",]
ICB_data<- ICB_data[ICB_data$preIpi==T,]

ICB_data$MGBS2_z<- (mean(ICB_data$MGBS2, na.rm=T)-ICB_data$MGBS2)/sd(ICB_data$MGBS2, na.rm=T) # Z-score normalization
ICB_data$MHC2_z<- (mean(ICB_data$MHC2, na.rm=T)-ICB_data$MHC2)/sd(ICB_data$MHC2, na.rm=T) # Z-score normalization
model<- glm(R ~ MHC2_z + MGBS2_z, data = ICB_data, family = binomial)
p_ls$model<- model
# summary(model)


# Roc from model
pred_prob <- predict(model, ICB_data, type = "response")
rocobj<- pROC::roc(response = ICB_data$R, predictor = pred_prob)
auc_value <- signif(pROC::auc(rocobj), 3)

p_ROC<- pROC::ggroc(rocobj, colour = 'steelblue', size = 2) +
  ggtitle("R ~ MGBS-II + MHC-II") + 
  geom_abline(intercept = 1, slope=1, linetype="dashed", size=1.5) +
  theme_classic() +
  annotate("text",x=-Inf, y=-Inf, label= paste0('AUC = ', auc_value), size= 8/.pt, hjust = 1, vjust = -1) +
  # annotate("text",x=-Inf, y=-Inf, label= paste0('cross-val AUC (mean) = ', signif(cv_auc, 2)), size= 8/.pt, vjust = -3, hjust = 1) +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    plot.title = element_text(size=8)
  )
p_ls$ROC<- p_ROC

# Predict response & plot survival
for(s in c("Liu_2019", "Riaz_2017", "Hugo_2016", "Gide_2019")){
  for(preIpi in c("all", T, F)){
    # Get data
    ICB_data<- ICB_study[!is.na(ICB_study$study)&ICB_study$study==s,]
    if(!s%in%c("Liu_2019", "Riaz_2017")&preIpi!="all") next
    if(preIpi!="all") ICB_data<- ICB_data[ICB_data$preIpi==preIpi,]
    # z score normalization
    ICB_data$MGBS2_z<- (mean(ICB_data$MGBS2, na.rm=T)-ICB_data$MGBS2)/sd(ICB_data$MGBS2, na.rm=T) # Z-score normalization
    ICB_data$MHC2_z<- (mean(ICB_data$MHC2, na.rm=T)-ICB_data$MHC2)/sd(ICB_data$MHC2, na.rm=T) # Z-score normalization
    # predict
    pred_prob <- predict(model, ICB_data, type = "response")
    mydat<- cbind(ICB_data, R_pred=pred_prob)
    mydat$R_pred_class<- mydat$R_pred > 0.5
    # Plot survival
    fit <- survfit(Surv(time = overall_survival, event = dead) ~ R_pred_class, data = mydat)
    legend_text<- paste0(c("Non-responder", "Responder"),paste0(" (n=",fit$n,")"))
    surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))
    p_surv<- ggsurvplot(fit, pval = TRUE, conf.int = F, xlab = "Time (years)", legend.labs=legend_text, legend.title="Prediction", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9))$plot
    # Add to list
    p_ls$predict_surv[[s]][[preIpi]]<- p_surv
    cat(s, preIpi, legend_text, "\n")
  }
}

p_ls$predict_surv$Liu_2019$all
p_ls$predict_surv$Liu_2019$'TRUE'
p_ls$predict_surv$Liu_2019$'FALSE'

p_ls$predict_surv$Riaz_2017$all
p_ls$predict_surv$Riaz_2017$'TRUE'
p_ls$predict_surv$Riaz_2017$'FALSE'

p_ls$predict_surv$Hugo_2016$all
p_ls$predict_surv$Gide_2019$all

# Liu_2019 all Non-responder (n=58) Responder (n=62) 
# Liu_2019 TRUE Non-responder (n=22) Responder (n=24) 
# Liu_2019 FALSE Non-responder (n=33) Responder (n=41) 
# Riaz_2017 all Non-responder (n=18) Responder (n=27) 
# Riaz_2017 TRUE Non-responder (n=10) Responder (n=11) 
# Riaz_2017 FALSE Non-responder (n=9) Responder (n=15) 
# Hugo_2016 all Non-responder (n=12) Responder (n=14) 
# Gide_2019 all Non-responder (n=13) Responder (n=22) 

# Boxplot heterogeneity
##########################

resp_t_df<- melt(ICB_study[,c("patient","preIpi","study","heterogeneity")])
resp_t_df<- resp_t_df[!is.na(resp_t_df$value),]

p <- ggboxplot(data = resp_t_df, x = "study", y = "value", col="preIpi",
               add = "jitter", xlab = F, ylab="Heterogeneity", add.params = list(size=0.5))+ 
  stat_compare_means(label = "p.format", label.x=1.5, hjust=0.5, size=7*0.35) +
  # facet_grid(.~variable) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7),
  )
p_ls$bp_heterogeneity<- p

# Save
#########
save.image(file="results/data/manuscript_validation_multivariate.RData")
saveRDS(p_ls, file="results/data/manuscript_validation_multivariate.rds")


