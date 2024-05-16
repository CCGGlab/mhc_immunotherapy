# Load the required packages
library(tidyverse)
library(scales)
library(survival)

# Load expression data and metadata
## Load ICB metadata in an encapsulated object "all_data"
all_data <- new.env()
load("data/MHC_immunotherapy.RData", all_data)
metadata <- all_data$ICB_study

# Global variables
## Lab plot theme
theme_ccgg <- theme_classic() +
  theme(
    # Remove x axis ticks and line
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    # Set y axis text size
    axis.text.y = element_text(size = 6),
    
    # Set legend position and properties
    legend.position = "bottom",
    legend.key.size =  unit(4, "mm"),
    legend.title = element_text(size = 7),
    
    # Set global text size
    text = element_text(size = 6),
    
    # Remove facet strip background and set text properties
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    
    # Remove spacing between facet panels
    panel.spacing.x = unit(0, "lines")
  )

## Mapping dictionaries for studies and predictors
studies_dict <- list(
  "Liu_2019" = "Liu et al., 2019",
  "Hugo_2016" = "Hugo et al., 2016",
  "Gide_2019" = "Gide et al., 2019",
  "Riaz_2017" = "Riaz et al., 2017"
)

predictors_binary_dict<- c(
  highTMB='High TMB',
  strongMHC1='High MGBS-I',
  strongMHC2='High MGBS-II',
  strongMHCnorm='High MGBS-d',
  highHetero='High Heterogeneity',
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
  highPDL1='High PD-L1 expr.'
)

predictors_continuous_dict <- c(
  TMB='TMB', 
  MGBS1='MGBS-I', 
  MGBS2='MGBS-II', 
  MGBSnorm='MGBS-d', 
  heterogeneity='Heterogeneity', 
  purity='Purity', 
  ploidy='Ploidy', 
  MHC1='MHC-I expr.', 
  MHC2 ='MHC-II expr.', 
  LDH='LDH', 
  LN_Met='LN Met', 
  hasLOHHLA='LOH HLA', 
  hasLOHB2M='LOH B2M', 
  MHC1_zyg='MHC-I Zygosity', 
  MHC2_zyg='MHC-II Zygosity', 
  CYT='Cyt. act.', 
  PD1='PD-L1 expr.'
)

# Helper functions
binarize_vars <- function(df, cutoff = 0.5) {
  # Mapping of original variable names to their new names
  var_name_mapping <- c(
    # With high prefix
    TMB = "highTMB",
    CYT = "highCYT",
    MHC1 = "highMHC1",
    MHC2 = "highMHC2",
    # Different capitalization
    purity = "highPurity",
    ploidy = "highPloidy",
    # Special cases with high:
    PD1 = "highPDL1",
    heterogeneity = "highHetero",
    # With "strong"
    MGBS1 = "strongMHC1",
    MGBS2 = "strongMHC2",
    MGBSnorm = "strongMHCnorm"
  )
  
  # Iterate over each variable and binarize based on the median
  for (var in names(var_name_mapping)) {
    y <- (df[[var]] > quantile(df[[var]], cutoff, na.rm = TRUE))
    binary_var_name <- var_name_mapping[[var]]
    df[[binary_var_name]] <- y
  }
  
  # Special case for LDH
  df$highLDH <- df$LDH == 1
  
  return(df)
}

univariate_coxph_analysis <- function(df, outcome_var, predictor_vars) {
  results <- map_dfr(predictor_vars, function(p) {
    # Remove rows where the predictor is missing
    v <- drop_na(df, !!p)
    
    if (nrow(v) == 0) {
      return(tibble(predictor = p))
    }
    
    # Fit Cox model
    fit.coxph <- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", p)), data = v)
    # Obtain results and convert to data frame
    res_stats <- broom::tidy(fit.coxph)
    
    HR <- exp(res_stats$estimate[[1]])
    pval <- res_stats$p.value[[1]]
    
    # Return the result
    tibble(predictor = p, HR, p.value = pval)
  })
  
  return(results)
}

univariate_logreg_analysis <- function(df, outcome_var, predictor_vars) {
  results <- map_dfr(predictor_vars, function(p) {
    # Remove rows where the predictor is missing
    v <- drop_na(df, !!p)
    
    if (nrow(v) == 0) {
      return(tibble(predictor = p))
    }
    
    # z score normalization of predictor variable
    v$z <- (mean(v[[p]])-v[[p]])/sd(v[[p]])
    
    # Fit logistic regression model
    res_glm <- glm(as.formula(paste0(outcome_var, " ~ z")), data = v, family = binomial())
    res_stats <- broom::tidy(res_glm)
    
    # Calculate OR and p value
    OR <- exp(res_stats$estimate[[2]])
    pval <- res_stats$p.value[[2]]
    
    # Calculate confidence intervals
    ci <- confint(res_glm, level = 0.95, method = "Wald")
    ci <- exp(ci[grep("z",rownames(ci)),])
    
    # Return the result
    tibble(predictor = p, OR, p.value = pval, !!!ci)
  })
  
  return(results)
}

# Process metadata
metadata_de <- metadata %>%
  filter(study %in% names(studies_dict)) %>%
  # Make sure the rows of the metadata correspond to the columns of combined_expr
  mutate(join_key = paste0(study, ".", patient)) %>%
  mutate(R = as.factor(R))

res_predictors_survival <- metadata_de %>%
  group_by(study, preIpi) %>%
  group_modify(function(df_study, study) {
    df_study <- binarize_vars(df_study)
    univariate_coxph_analysis(df_study, "overall_survival", names(predictors_binary_dict))
  })

res_predictors_response <- metadata_de %>%
  group_by(study, preIpi) %>%
  group_modify(function(df_study, study) {
    univariate_logreg_analysis(df_study, "R", names(predictors_continuous_dict))
  })

# Visualization
## Create the heatmap
plot_univariate_analysis <- function(df, outcome_var, predictor_dict, title) {
  df <- df %>%
    mutate(preIpi = factor(
      # Set labels for pretreatment status as wanted
      if_else(preIpi, "Pretreated", "Treatment naive") %>% replace_na("Unknown"),
      # Define order (for facets)
      levels = c("Treatment naive", "Pretreated", "Unknown")
    ),
    study = factor(
      unlist(studies_dict[study]),
      levels = c(
        "Liu et al., 2019",
        "Riaz et al., 2017",
        "Hugo et al., 2016",
        "Gide et al., 2019"
      )
    )) %>%
    ungroup %>%
    arrange(study, factor(preIpi, levels = c("Pretreated", "Treatment naive", "Unknown")), desc(!!sym(outcome_var)), p.value) %>%
    # Rename the predictors and choose the current order as the order of levels
    mutate(predictor = unlist(predictor_dict[predictor])) %>%
    mutate(predictor = factor(predictor, levels = unique(predictor), ordered = T))
  
  plt <- ggplot(df, aes(
    x = study,
    y = predictor,
    fill = !!sym(outcome_var),
    colour = if_else(is.na(p.value), "FALSE", as.character(p.value < 0.05)),
    size = abs(log10(!!sym(outcome_var)))
  )) +
    scale_x_discrete() +
    scale_y_discrete(limits = rev) +
    theme_ccgg +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust=1, vjust=1)
    ) +
    scale_fill_gradient2(
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      trans = scales::log_trans(base = 10),
      breaks = scales::log_breaks(n = 5, base = 10),
      midpoint = 0,
      na.value = NA,
      limits = c(0.5, 2),
      oob = scales::oob_squish
    ) +
    scale_colour_manual(
      breaks = c("TRUE", "FALSE"),
      # green or transparent
      values = c(scales::muted("green", l = 100, c = 127), "#FFFFFF00"),
      name = "p-value < 0.05"
    ) +
    scale_size(guide = guide_none()) +
    geom_point(pch = 21, stroke = 2) +
    facet_wrap(~ preIpi, scales = "free_x") +
    ggtitle(title)
  
  return(plt)
}

# Generate and save plot for survival
plt_predictors_survival <- plot_univariate_analysis(
  res_predictors_survival,
  "HR",
  predictors_binary_dict,
  "Survival"
)

# Generate and save plot for response
plt_predictors_response <- plot_univariate_analysis(
  res_predictors_response,
  "OR",
  predictors_continuous_dict,
  "Response"
)

p_ls<- list(pred_surv=plt_predictors_survival, pred_R=plt_predictors_response)

# Save
#########
save.image(file="results/data/manuscript_validation_compare_variables.RData")
saveRDS(p_ls, file="results/data/manuscript_validation_compare_variables.rds")

