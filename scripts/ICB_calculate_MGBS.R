# Calculate MGBS scores, defined as the proportion of peptides with Kd < 500 nM
#################################################################################

# Load MHC genotypes
MHC_gt <- as.data.frame(readRDS("temp/hla_calls.rds")$Liu_2019)

# Load affinities per allele
MHC1<- readRDS("temp/mhc1_rand_matrix.rds")
MHC2<- readRDS("temp/mhc2_rand_matrix.rds")

# Define scores (proportion MHC binders) per allele
MGBS1_alleles<- colMeans(MHC1<500)
MGBS2_alleles<- colMeans(MHC2<500)

# MHC-I genotypes in same format as affinities
for(a in c("A","B","C")){
  for(i in 1:2) MHC_gt[[paste0(a,".", i)]]<- paste0("HLA-",a, MHC_gt[[paste0(a,".", i)]])
}
rownames(MHC_gt)<- MHC_gt$patient
MHC1_gt<- MHC_gt[,2:7]

# Get MGBS-I
MGBS1<- apply(MHC1_gt, 1, function(x) mean(MGBS1_alleles[as.character(x)]))

# MHC-II genotypes in same format as affinities
# First get heterodimers
for(i in grep("D",colnames(MHC_gt))) MHC_gt[[i]]<- gsub("\\:","",MHC_gt[[i]])
MHC2_gt<- data.frame(
  row.names=MHC_gt$patient,
  HLA_DP_11=paste0("HLA-","DPA1",MHC_gt[["DPA1.1"]],"-","DPB1",MHC_gt[["DPB1.1"]]), 
  HLA_DP_12=paste0("HLA-","DPA1",MHC_gt[["DPA1.1"]],"-","DPB1",MHC_gt[["DPB1.2"]]), 
  HLA_DP_21=paste0("HLA-","DPA1",MHC_gt[["DPA1.2"]],"-","DPB1",MHC_gt[["DPB1.1"]]), 
  HLA_DP_22=paste0("HLA-","DPA1",MHC_gt[["DPA1.2"]],"-","DPB1",MHC_gt[["DPB1.2"]]), 
  HLA_DQ_11=paste0("HLA-","DQA1",MHC_gt[["DQA1.1"]],"-","DQB1",MHC_gt[["DQB1.1"]]), 
  HLA_DQ_12=paste0("HLA-","DQA1",MHC_gt[["DQA1.1"]],"-","DQB1",MHC_gt[["DQB1.2"]]), 
  HLA_DQ_21=paste0("HLA-","DQA1",MHC_gt[["DQA1.2"]],"-","DQB1",MHC_gt[["DQB1.1"]]), 
  HLA_DQ_22=paste0("HLA-","DQA1",MHC_gt[["DQA1.2"]],"-","DQB1",MHC_gt[["DQB1.2"]]), 
  HLA_DR_1=paste0("DRB1_",MHC_gt[["DRB1.1"]]),
  HLA_DR_2=paste0("DRB1_",MHC_gt[["DRB1.2"]])
)

# Get MGBS-II
MGBS2<- apply(MHC2_gt, 1, function(x) mean(MGBS2_alleles[as.character(x)]))

# Merge results
MGBS_df<- data.frame(
  patient = names(MGBS1),
  MGBS1 = MGBS1,
  MGBS2 = MGBS2[names(MGBS1)]
)

# Add MGBS-d
calculate_MGBS_norm<- function(df, var_MHC1, var_MHC2){
  
  library(tidyverse)
  
  # Get aff vectors
  m1<- df[,var_MHC1]
  m2<- df[,var_MHC2]
  
  # ECDF function
  m1f <- ecdf(discard(m1, is.na))
  m2f <- ecdf(discard(m2, is.na))
  
  # Calculate
  mdnorm <- set_names(m1f(m1) - m2f(m2), rownames(df))
  
  # Return
  return(mdnorm)
  
}
MGBS_df$MGBSd<- calculate_MGBS_norm(MGBS_df, "MGBS1","MGBS2")

# Save
saveRDS(MGBS_df,  "temp/MGBS_scores.rds")
