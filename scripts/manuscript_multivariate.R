# Manuscript multivariate
##########################

# Check added value of MGBS-II & related as compared to other previously described biomarkers

library(survminer)
library(survival)
library(reshape2)
library(cowplot)
library(caret)
library(pROC)

# Load data
############
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]

# Remove NA
############
ICB_data<- na.omit(ICB_data)

# Binarize
#############
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)
ICB_data$highPD1<- ICB_data$PD1_noBatch > quantile(ICB_data$PD1_noBatch, 0.5, na.rm=T)
ICB_data$highCYT<- ICB_data$CYT_noBatch > quantile(ICB_data$CYT_noBatch, 0.5, na.rm=T)
ICB_data$highMHC1<- ICB_data$MHC1_noBatch > quantile(ICB_data$MHC1_noBatch, 0.5, na.rm=T)
ICB_data$highMHC2<- ICB_data$MHC2_noBatch > quantile(ICB_data$MHC2_noBatch, 0.5, na.rm=T)
ICB_data$highPurity<- ICB_data$purity > quantile(ICB_data$purity, 0.5, na.rm=T)
ICB_data$highPloidy<- ICB_data$ploidy > quantile(ICB_data$ploidy, 0.5, na.rm=T)
ICB_data$highHeterogeneity<- ICB_data$heterogeneity > quantile(ICB_data$heterogeneity, 0.5, na.rm=T)
ICB_data$highLDH<- ICB_data$LDH==1 

# Cox
######

Cox_model_ls<- list()

ICB_data$hasLOHB2M<- as.numeric(ICB_data$hasLOHB2M)

# Multivariate= PD1 & StrongMHC2
vars<- c("highTMB","strongMHC1","strongMHC2","strongMHCnorm","highHeterogeneity","highPurity","highPloidy", "highMHC1", "highMHC2", "highLDH", "LN_Met","hasLOHB2M","hasLOHHLA","highMHC1_zyg","highMHC2_zyg","highCYT","highPD1")
vars_dict<- c(
  highTMB='High TMB', 
  strongMHC1='High MGBS-I', 
  strongMHC2='High MGBS-II', 
  strongMHCnorm='High MGBS-d', 
  highHeterogeneity='High Heterogeneity', 
  highPurity='High Purity', 
  highPloidy='High Ploidy', 
  highMHC1='High MHC-I expr.', 
  highMHC2='High MHC-II expr.', 
  highLDH='High LDH', 
  LN_Met='LN Met', 
  hasLOHB2M='LOH B2M', 
  hasLOHHLA='LOH HLA', 
  highMHC1_zyg='High MHC-I Zygosity', 
  highMHC2_zyg='High MHC-II Zygosity', 
  highCYT='High cyt. act.', 
  highPD1='High PD-L1 expr.'
  )

# Univariate
HR_uni<- NULL
for(v in vars){
  for(isPreIpi in c(F,T)){
    fit_coxph<- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", v)), data = subset(ICB_data,preIpi==isPreIpi))
    HR_uni<- rbind(HR_uni, c(v, isPreIpi,summary(fit_coxph)$conf.int[,c(1,3,4)], summary(fit_coxph)$coefficients[,"Pr(>|z|)"]))
  }
}
colnames(HR_uni)<- c("variable","preIpi","HR","ci_l","ci_h","p")
HR_uni<- as.data.frame(HR_uni)
HR_uni$HR<- round(as.numeric(HR_uni$HR),2)
HR_uni$ci_l<- round(as.numeric(HR_uni$ci_l),2)
HR_uni$ci_h<- round(as.numeric(HR_uni$ci_h),2)
HR_uni$p<- signif(as.numeric(HR_uni$p),2)
HR_uni$HR_ci<- paste0(HR_uni$HR, " (",HR_uni$ci_l, "-",HR_uni$ci_h,")")

HR_uni_preIpi<- HR_uni[HR_uni$preIpi==T,c("variable","HR_ci","p")]
HR_uni_naive<- HR_uni[HR_uni$preIpi==F,c("variable","HR_ci","p")]

HR_uni_preIpi<- HR_uni_preIpi[order(HR_uni_preIpi$p),]
HR_uni_naive<- HR_uni_naive[order(HR_uni_naive$p),]

t1<- ggtexttable(HR_uni_naive[,-1], rows = vars_dict[HR_uni_naive$variable], cols = c("HR (95% CI)","P"), theme = ttheme("light",base_size = 8,padding = unit(c(2, 2), "mm")))
t1<- tab_add_title(t1, text = "Treatment-naive", face = "bold", size=10)
t2<- ggtexttable(HR_uni_preIpi[,-1], rows = vars_dict[HR_uni_preIpi$variable], cols = c("HR (95% CI)","P"), theme = ttheme("light",base_size = 8,padding = unit(c(2, 2), "mm")))
t2<- tab_add_title(t2, text = "Ipilimumab pretreated", face = "bold", size=10)

Cox_model_ls[["FALSE"]][["t_uni"]]<- t1
Cox_model_ls[["TRUE"]][["t_uni"]]<- t2

# Multivariate FW model
HR_multi_df<- NULL
for(isPreIpi in c(T,F)){
  HR_tmp<- HR_uni[HR_uni$preIpi==isPreIpi,]
  HR_tmp<- HR_tmp[order(HR_tmp$p),]
  vars_incl<- HR_tmp$variable[1]
  p_var_tmp<- NULL
  p_var_to_consider<- 0
  while(p_var_to_consider<0.05){
    for(v in vars){
      if(v %in% vars_incl) next # Leave norm out of analysis, dependent on MHCI & MHCII?
      vars_tmp<- c(vars_incl,v)
      fit_coxph<- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", paste0(vars_tmp, collapse = " + "))), data = subset(ICB_data,preIpi==isPreIpi))
      p_tmp<- rev(summary(fit_coxph)$coefficients[,"Pr(>|z|)"])[1]
      names(p_tmp)<- v
      p_var_tmp<- c(p_var_tmp, p_tmp)
    }
    p_var_to_consider<- sort(p_var_tmp)[1]
    vars_incl<- c(vars_incl, names(p_var_to_consider))
    p_var_tmp<- NULL
  }
  vars_incl<- rev(vars_incl)[-1]
  fit_coxph<- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", paste0(vars_incl, collapse = " + "))), data = subset(ICB_data,preIpi==isPreIpi))
  # summary(fit_coxph)
  
  # Add to df
  HR_multi_df_tmp<- data.frame(
    cond=vars_incl,
    HR=summary(fit_coxph)$conf.int[,1],
    ci_l=summary(fit_coxph)$conf.int[,3],
    ci_h=summary(fit_coxph)$conf.int[,4],
    p=summary(fit_coxph)$coefficients[,5],
    preIpi=isPreIpi
  )
  if(is.null(HR_multi_df)) HR_multi_df<- HR_multi_df_tmp
  else HR_multi_df<- rbind(HR_multi_df, HR_multi_df_tmp)
}

HR_multi_df$ipiLMGBS<- "Ipilimumab pretreated"
HR_multi_df$ipiLMGBS[HR_multi_df$preIpi==F]<- "Treatment naive"
HR_multi_df$ipiLMGBS<- factor(HR_multi_df$ipiLMGBS, rev(unique(HR_multi_df$ipiLMGBS)))
HR_multi_df$condLMGBS<- vars_dict[HR_multi_df$cond]

fp_multi <- ggplot(data=HR_multi_df, aes(x=condLMGBS, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.5, fatten = 10, shape = 18) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-2, size=7/3) +
  geom_text(aes(y=6.5, label=paste0("P=",signif(p,2))), fontface="italic", vjust=2, hjust=0, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_wrap(.~ipiLMGBS, scales = "free") +  
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.15,10), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
fp_multi

Cox_model_ls[["fp_multi"]]<- fp_multi

# Logistic regression
####################

LR_model_ls<- list()

# Univariate
vars<- c("TMB","MGBS1","MGBS2","MGBSnorm","heterogeneity","purity","ploidy", "MHC1_noBatch", "MHC2_noBatch", "LDH", "LN_Met", "hasLOHHLA", "hasLOHB2M", "MHC1_zyg", "MHC2_zyg", "CYT_noBatch","PD1_noBatch")
vars_dict<- c(
  TMB='TMB', 
  MGBS1='MGBS-I', 
  MGBS2='MGBS-II', 
  MGBSnorm='MGBS-d', 
  heterogeneity='Heterogeneity', 
  purity='Purity', 
  ploidy='Ploidy', 
  MHC1_noBatch='MHC-I expr.', 
  MHC2_noBatch='MHC-II expr.', 
  LDH='LDH', 
  LN_Met='LN Met', 
  hasLOHHLA='LOH HLA', 
  hasLOHB2M='LOH B2M', 
  MHC1_zyg='MHC-I Zygosity', 
  MHC2_zyg='MHC-II Zygosity', 
  CYT_noBatch='Cyt. act.', 
  PD1_noBatch='PD-L1 expr.'
)

LR_uni_df<- NULL
for(v in vars){
  ICB_data_tmp<- ICB_data[!is.na(ICB_data[[v]]),]
  ICB_data_tmp$z<- (mean(ICB_data_tmp[[v]])-ICB_data_tmp[[v]])/sd(ICB_data_tmp[[v]]) # Z-score normalization
  for(isPreIpi in c(F,T)){
    # z-score normalize to standardize the scale
    model <- glm(as.formula(paste0("R ~ ", "z")), data = ICB_data_tmp[ICB_data_tmp$preIpi==isPreIpi,], family = binomial)
    OR_tmp<- exp(summary(model)$coefficients["z","Estimate"])
    p_tmp<- summary(model)$coefficients["z","Pr(>|z|)"]
    ci<- confint(model, level = 0.95, method = "Wald")
    ci<- exp(ci[grep("z",rownames(ci)),])
    df_tmp<- data.frame(
      variable = v,
      preIpi = isPreIpi,
      OR = OR_tmp,
      p = p_tmp,
      ci_l = ci[1],
      ci_h = ci[2]
    )
    LR_uni_df<- rbind(LR_uni_df,df_tmp)
  }
}

LR_uni_df$OR<- round(as.numeric(LR_uni_df$OR),2)
LR_uni_df$ci_l<- round(as.numeric(LR_uni_df$ci_l),2)
LR_uni_df$ci_h<- round(as.numeric(LR_uni_df$ci_h),2)
LR_uni_df$p<- signif(as.numeric(LR_uni_df$p),2)
LR_uni_df$OR_ci<- paste0(LR_uni_df$OR, " (",LR_uni_df$ci_l, "-",LR_uni_df$ci_h,")")

LR_uni_preIpi<- LR_uni_df[LR_uni_df$preIpi==T,c("variable","OR_ci","p")]
LR_uni_naive<- LR_uni_df[LR_uni_df$preIpi==F,c("variable","OR_ci","p")]

LR_uni_preIpi<- LR_uni_preIpi[order(LR_uni_preIpi$p),]
LR_uni_naive<- LR_uni_naive[order(LR_uni_naive$p),]

t1<- ggtexttable(LR_uni_naive[,-1], rows = vars_dict[LR_uni_naive$variable], cols = c("OR (95% CI)","P"), theme = ttheme("light",base_size = 8,padding = unit(c(2, 2), "mm")))
t1<- tab_add_title(t1, text = "Treatment-naive", face = "bold", size=10)
t2<- ggtexttable(LR_uni_preIpi[,-1], rows = vars_dict[LR_uni_preIpi$variable], cols = c("OR (95% CI)","P"), theme = ttheme("light",base_size = 8,padding = unit(c(2, 2), "mm")))
t2<- tab_add_title(t2, text = "Ipilimumab pretreated", face = "bold", size=10)

LR_model_ls[["FALSE"]][["t_uni"]]<- t1
LR_model_ls[["TRUE"]][["t_uni"]]<- t2

# Multivariate Forward model
for(isPreIpi in c(T,F)){
  
  # First step
  full_model <- glm(as.formula(paste0("R ~ ", paste0(vars, collapse = " + "))), data = ICB_data[ICB_data$preIpi==isPreIpi,] , family = binomial) #SIGN
  nothing <- glm(R ~ 1, data = ICB_data[ICB_data$preIpi==isPreIpi,], family=binomial)
  fw<- step(nothing,scope=list(lower=formula(nothing),upper=formula(full_model)), direction="forward")
  # summary(fw)
  
  # Iterate untill all vars significant
  var_coefs<- summary(fw)$coefficients[-1,] 
  var_to_excl<- rownames(var_coefs)[var_coefs[,"Pr(>|z|)"]==max(var_coefs[,"Pr(>|z|)"])]
  while(length(var_to_excl)>0){
    full_model <- glm(as.formula(paste0("R ~ ", paste0(setdiff(gsub("TRUE","",rownames(var_coefs)),var_to_excl), collapse = " + "))), data = ICB_data[ICB_data$preIpi==isPreIpi,] , family = binomial) #SIGN
    fw<- step(nothing,scope=list(lower=formula(nothing),upper=formula(full_model)), direction="forward")
    var_coefs<- summary(fw)$coefficients[-1,] 
    var_to_excl<- rownames(var_coefs)[var_coefs[,"Pr(>|z|)"]==max(var_coefs[,"Pr(>|z|)"])&var_coefs[,"Pr(>|z|)"]>0.05]
    # summary(fw)
  }
  LR_model_ls[[as.character(isPreIpi)]][["fw"]]<- fw
}

# Compare AUCs
################

# Helper functions
kfold_glm <- function(sel_data, k, model_formula) {
  folds <- caret::createFolds(sel_data$R, k = k)
  
  # Initialize a vector to store the predictions and ground truth for each fold
  results <- list()
  
  # Loop over the folds
  for (i in 1:k) {
    # Define the indices for the test set
    test_indices <- folds[[i]]
    
    # Split the data into train and test sets
    train_data <- sel_data[-test_indices, ]
    test_data <- sel_data[test_indices, ]
    
    # Train the glm model with binomial family
    model <- glm(model_formula, data = train_data, family = "binomial")
    
    # Predict the probabilities for the test set
    pred_prob <- predict(model, test_data, type = "response")
    
    # Store the predictions and ground truth
    results[[i]] <- list(predictions = pred_prob, truth = test_data$R)
  }
  
  results
}

get_roc_auc <- function(l) {
  roc_obj <- pROC::roc(response = l$truth, predictor = l$predictions)
  pROC::auc(roc_obj)
}

construct_formula <- function(var_names) {
  lhs <- "R"
  rhs <- paste0(var_names, collapse=" + ")
  paste0(lhs, " ~ ", rhs)
}

get_cv_results <- function(model_formula, sel_data) {
  # Repeat the 5-fold cross validation 100 times 
  res <- replicate(n = 100, kfold_glm(sel_data, k = 5, model_formula = model_formula))
  
  # Obtain average cross-validation ROC AUC for each of the 100 repetitions
  res_auc <- apply(res, MARGIN = 2, function(x) {
    # mean over 5 folds
    mean(sapply(x, get_roc_auc))
  })
  
  list(auc = res_auc, all_data = res)
}

calc_roc_auc_ci <- function(cv_results) {
  # Calculate the mean and standard error (of sample mean)
  mean <- mean(cv_results)
  se <- sd(cv_results) / sqrt(length(cv_results))
  
  # Calculate the lower and upper limits of the 95% confidence interval
  lower <- mean - 1.96 * se
  upper <- mean + 1.96 * se
  
  # Return a named vector with all calculated statistics
  c(mean = mean, se = se, lower = lower, upper = upper)
}

format_model_names <- function(model_names) {
  # split the ind column by " + " and store it as a list
  ind_list <- strsplit(as.character(model_names), " \\+ ")
  
  # loop over model names
  ind_list <- lapply(ind_list, function(x) ifelse(x %in% names(vars_dict), vars_dict[x], x))
  
  # concatenate the list elements by " + " and store it as a vector
  ind_clean <- sapply(ind_list, paste, collapse = " + ")
  
  ind_clean
}

# Models to compare
ffw_vars <-
  list(
    `FALSE` = c("heterogeneity", "ploidy", "MGBS1"),
    `TRUE` = c("MGBS2", "hasLOHB2M", "MHC2_noBatch", "heterogeneity")
  )

ffw_vars_rev <-
  list(
    `TRUE` = c("heterogeneity", "ploidy", "MGBS1"),
    `FALSE` = c("MGBS2", "hasLOHB2M", "MHC2_noBatch", "heterogeneity")
  )

MHC_vars <-
  list(
    `FALSE` = c("MGBS2", "MHC2_noBatch"),
    `TRUE` = c("MGBS2", "MHC2_noBatch")
  )

ffw_nomgbs_vars <-
  list(
    `FALSE` = c("heterogeneity", "ploidy"),
    `TRUE` = c("hasLOHB2M", "MHC2_noBatch", "heterogeneity")
  )

liu_vars <-
  list(
    `FALSE` = c("heterogeneity", "ploidy", "purity"),
    `TRUE` = c("MHC2_noBatch", "LN_Met", "LDH")
  )

# Create an empty list to store the results
df_roc_auc_lst <- list()

# Loop over the preIpi conditions
for (isPreIpi in c(T, F)) {
  # Data for the considered preIpi condition
  sel_data <- ICB_data[ICB_data$preIpi == isPreIpi,]
  
  # Multivariate model of Liu
  liu_model <- construct_formula(liu_vars[[as.character(isPreIpi)]])
  
  # Best model after feed forward feature regression
  # ffw_model <- construct_formula(ffw_vars[[as.character(isPreIpi)]])
  
  # Best model according to fig S4 (FFW + removal of non-influential variables)
  ffw_model <- construct_formula(ffw_vars[[as.character(isPreIpi)]])
  
  ffw_nomgbs_model <- construct_formula(ffw_nomgbs_vars[[as.character(isPreIpi)]])

  # Other
  MHC2_model <- construct_formula(MHC_vars[[as.character(isPreIpi)]])
  
  ffw_rev_model <- construct_formula(ffw_vars_rev[[as.character(isPreIpi)]])
  
  # Univariate models
  univariate_models <- unname(sapply(vars, construct_formula))
  
  # Create list with all models
  models <- c(univariate_models, liu_model, ffw_model, ffw_nomgbs_model,MHC2_model,ffw_rev_model)
  model_names <- model_names <- sub("^R ~ ", "", models)
  
  # Cross-validation and calculation of ROC AUC
  cv_results <- lapply(models, get_cv_results, sel_data)
  names(cv_results) <- model_names
  roc_auc_cv <- lapply(cv_results, function(x) x$auc)
  # all_data[[as.character(isPreIpi)]] <- lapply(cv_results, function(x) x$all_data)
  
  # Convert the list of values to a data frame (stacked by model)
  df_cv <- stack(roc_auc_cv)
  
  # Obtain the confidence interval of the ROC AUC for each model
  df_roc_auc <- aggregate(values ~ ind, df_cv, calc_roc_auc_ci)
  
  # Add the preIpi condition as a column
  if (isPreIpi) {
    df_roc_auc$isPreIpi <- "Ipilimumab pretreated"
  } else {
    df_roc_auc$isPreIpi <- "Treatment-naive"
  }
  
  # Rename the single variable models on the bar plot to clean names
  levels(df_roc_auc$ind) <- format_model_names(levels(df_roc_auc$ind))
  
  # Store the data frame in the list
  df_roc_auc_lst[[as.character(isPreIpi)]] <- df_roc_auc
}

# Combine the data frames for preIpi True and False
# (to create one plot with facets)
df_combined <- do.call(rbind, df_roc_auc_lst)

# Make the data frame structure more convenient to work with
df_clean <- data.frame(model_name = df_combined$ind, 
                       isPreIpi = df_combined$isPreIpi, 
                       as.data.frame(df_combined$values))

# Plot
df_clean$model_name <- reorder(interaction(df_clean$isPreIpi, df_clean$model_name, sep=":", drop=T), df_clean$mean)
df_clean$isPreIpi <- factor(df_clean$isPreIpi, levels = c("Treatment-naive", "Ipilimumab pretreated"))

AUC_bp<- ggplot(df_clean, aes(x = model_name, y = mean)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.9)) +
  geom_text(aes(label = round(mean, 2), y = upper), position = position_dodge(0.9), vjust = -0.5, size=6/3) +
  labs(x = "", y = "ROC AUC") +
  facet_wrap(.~isPreIpi, scales = "free_x") +
  scale_x_discrete(label = function(x) {gsub('^[^:]*:', '', x)}) +
  coord_cartesian(ylim = c(0.4, 1)) +  
  theme_pubclean() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7)
    )
LR_model_ls[["AUC_bp"]]<- AUC_bp

AUC_bp_preIpi<- ggplot(subset(df_clean, isPreIpi=="Ipilimumab pretreated"), aes(x = model_name, y = mean)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.9)) +
  geom_text(aes(label = round(mean, 2), y = upper), position = position_dodge(0.9), vjust = -0.5, size=6/3) +
  labs(x = "", y = "ROC AUC") +
  scale_x_discrete(label = function(x) {gsub('^[^:]*:', '', x)}) +
  coord_cartesian(ylim = c(0.5, 0.9)) +  
  theme_pubclean() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7)
  )
LR_model_ls[["TRUE"]][["AUC_bp"]]<- AUC_bp_preIpi

AUC_bp_noIpi<- ggplot(subset(df_clean, isPreIpi=="Treatment-naive"), aes(x = model_name, y = mean)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.9)) +
  geom_text(aes(label = round(mean, 2), y = upper), position = position_dodge(0.9), vjust = -0.5, size=6/3) +
  labs(x = "", y = "ROC AUC") +
  scale_x_discrete(label = function(x) {gsub('^[^:]*:', '', x)}) +
  coord_cartesian(ylim = c(0.4, 0.9)) +  
  theme_pubclean() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 7)
  )
LR_model_ls[["FALSE"]][["AUC_bp"]]<- AUC_bp_noIpi


# ROC + KM curve for entire dataset
######################################
full_roc <- list()
for (isPreIpi in c(F, T)) {
  # Data for the considered preIpi condition
  sel_data <- ICB_data[ICB_data$preIpi == isPreIpi,]
  
  model_formula <- construct_formula(ffw_vars[[as.character(isPreIpi)]])
  model_name <- format_model_names(list(sub("^R ~ ", "", model_formula)))[[1]]
  
  # Train the glm model with binomial family
  model <- glm(model_formula, data = sel_data, family = "binomial")
  
  # Predict the probabilities for the test set
  pred_prob <- predict(model, sel_data, type = "response")
  
  rocobj <- pROC::roc(response = sel_data[,"R"], predictor = pred_prob)
  auc_value <- signif(pROC::auc(rocobj), 2)
  cv_aucs <- df_roc_auc_lst[[as.character(isPreIpi)]]
  cv_auc <- cv_aucs[cv_aucs$ind == model_name,"values"][[1,"mean"]]
  
  p_ROC<- pROC::ggroc(rocobj, colour = 'steelblue', size = 2) +
    ggtitle(model_name) + 
    geom_abline(intercept = 1, slope=1, linetype="dashed", size=1.5) +
    theme_classic() +
    annotate("text",x=-Inf, y=-Inf, label= paste0('AUC = ', auc_value), size= 8/.pt, hjust = 1, vjust = -1) +
    annotate("text",x=-Inf, y=-Inf, label= paste0('cross-val AUC (mean) = ', signif(cv_auc, 2)), size= 8/.pt, vjust = -3, hjust = 1) +
    theme(
      axis.text = element_text(size=6),
      axis.title = element_text(size=7),
      plot.title = element_text(size=8)
    )
  LR_model_ls[[as.character(isPreIpi)]][["p_ROC"]]<- p_ROC

  # Survival
  mydat<- cbind(sel_data,R_pred=pred_prob)
  mydat$R_pred_class<- mydat$R_pred > 0.5
  fit <- survfit(Surv(time = overall_survival, event = dead) ~ R_pred_class, data = mydat)
  legend_text<- paste0(c("Non-responder", "Responder"),paste0(" (n=",fit$n,")"))
  surv_cols<- c(rgb(231/256, 184/256, 0/256),rgb(46/256, 159/256, 223/256))
  p_surv<- ggsurvplot(fit, pval = TRUE, conf.int = T, xlab = "Time (years)", legend.labs=legend_text, legend.title="Prediction", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=7, font.title=8, legend=c(0.8,0.9))$plot
  LR_model_ls[[as.character(isPreIpi)]][["p_surv"]]<- p_surv
}

# Save
######
save.image(file="results/data/manuscript_multivariate.RData")
saveRDS(list(Cox_model_ls=Cox_model_ls, LR_model_ls=LR_model_ls), file="results/data/manuscript_multivariate.rds")
# load(file="results/data/manuscript_multivariate.RData")
