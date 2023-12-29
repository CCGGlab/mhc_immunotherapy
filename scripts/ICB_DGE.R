library(tidyverse)

# Get MGBS scores
load("data/MHC_immunotherapy.RData")
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]

# Get downloaded expression data
ICB_expr<- as.data.frame(read_tsv("downloads/pub/campbell_2023/MORRISON-1-public/RNASeq/data/RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv.zip"))
gene_names<- ICB_expr$gene.hgnc.symbol

# Standardize sample IDs
ICB_meta<- read.table("downloads/pub/campbell_2023/MORRISON-1-public/RNASeq/RNA-CancerCell-MORRISON1-metadata.tsv", header=T, row.names = 1)
colnames(ICB_expr)<- ICB_meta[colnames(ICB_expr),"subject.id"]
ICB_expr<- ICB_expr[,grep("_liu", colnames(ICB_expr))]         
colnames(ICB_expr)<- gsub("_liu","",colnames(ICB_expr))
rownames(ICB_expr)<- gene_names

# Stratify in 2 groups
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)

# Remove missing data
ICB_data<- ICB_data[!is.na(ICB_data$strongMHC2),]

# Only patients that contain expression information
common_patients<- intersect(ICB_data$patient, colnames(ICB_expr))
ICB_data<- ICB_data[ICB_data$patient%in%common_patients,]
ICB_expr<- ICB_expr[,ICB_data$patient]

# Limma
library(limma)
DGE_ls<- list()
vars<- c("highTMB","strongMHC1", "strongMHC2","strongMHCnorm")
for(v in vars){
  cat(v," ")
  for(preIpi in c(NA,T,F)){
    if(is.na(preIpi)) ICB_data_tmp<- ICB_data
    else ICB_data_tmp<- ICB_data[ICB_data$preIpi==preIpi,]
    design <- model.matrix(~ICB_data_tmp[[v]])
    ICB_expr_tmp<- ICB_expr[,ICB_data_tmp$patient]
    fit<- lmFit(ICB_expr_tmp, design)
    fit<- eBayes(fit, trend=TRUE)
    res<- topTable(fit, coef=ncol(design), number=Inf)
    
    # GSEA
    stat<- res$t
    names(stat)<- rownames(res)
    stat<- sort(stat,decreasing = T)
    
    GSEA_res<- fgsea::fgseaMultilevel(pathways = Ha_ls,
                                      stats = stat,
                                      minSize=15,
                                      maxSize=500,
                                      eps=0,
                                      nPermSimple = 10000,
                                      nproc=1)
    GSEA_res<- GSEA_res[order(GSEA_res$pval),]
    
    # Save
    if(is.na(preIpi)) preIpi<- "all"
    DGE_ls[[v]][[as.character(preIpi)]][["DGE"]]<- res
    DGE_ls[[v]][[as.character(preIpi)]][["GSEA"]]<- GSEA_res
  }
}

# Save
saveRDS(DGE_ls,"temp/DGE_ls.rds")
