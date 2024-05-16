load("data/MHC_immunotherapy.RData")
library(tidyverse)

# Added: mapping table as ENSG -> HGNC must occur here (created by create_ENSG_HGNC_map_table.R)
ENSG_HGNC_map_table <- ENSG_HGNC_map_table %>%
  select(gene_id = ensembl_gene_id, gene = external_gene_name, gene_biotype)

# Get MGBS scores for Liu
ICB_data<- ICB_study[ICB_study$study=="Liu_2019",]

# Convert to txi object to DGEList for Limma
dgelist <- edgeR::DGEList(txi_object)

# Stratify in 2 groups
ICB_data$highTMB<- ICB_data$TMB > quantile(ICB_data$TMB, 0.5, na.rm=T)
ICB_data$strongMHC1<- ICB_data$MGBS1 > quantile(ICB_data$MGBS1, 0.5, na.rm=T)
ICB_data$strongMHC2<- ICB_data$MGBS2 > quantile(ICB_data$MGBS2, 0.5, na.rm=T)
ICB_data$strongMHCnorm<- ICB_data$MGBSnorm > quantile(ICB_data$MGBSnorm, 0.5, na.rm=T)

# Remove missing data
ICB_data<- ICB_data[!is.na(ICB_data$strongMHC2),]

# Only patients that contain expression information
common_patients<- intersect(ICB_data$patient, colnames(txi_object))
ICB_data<- ICB_data[ICB_data$patient%in%common_patients,]

# Limma
library(edgeR) # needed for handling DGEList objects and preprocessing functions
library(limma)
DGE_ls<- list()
vars<- c("highTMB","strongMHC1", "strongMHC2","strongMHCnorm")
for(v in vars){
  cat(v," ")
  for(preIpi in c(NA,T,F)){
    if(is.na(preIpi)) {
      ICB_data_tmp<- ICB_data
    } else {
      ICB_data_tmp<- ICB_data[ICB_data$preIpi==preIpi,]
    }
    
    design <- model.matrix(~ICB_data_tmp[[v]])
    
    # Load the correct DGEList and take the requested subset
    dge <- dgelist[,ICB_data_tmp$patient]
    # Remove rows that consistently have zero or very low counts
    # (As recommended in Limma manual: 15.3)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, ]

    # Normalize and run voom transformation (15.3 and 15.5)
    dge <- calcNormFactors(dge)
    dge <- voom(dge, design)
    # dge is now ready for lmFit()
    
    # Below, just changed ICB_expr to "dge"
    fit<- lmFit(dge, design)
    # We now applied the voom transform do not add "trend"
    fit<- eBayes(fit)
    # Obtain results and map gene IDs (from ENSG to HGNC)
    res<- topTable(fit, coef=ncol(design), number=Inf) %>%
      # We still have Ensembl IDs, make the mapping here
      rownames_to_column("gene_id") %>%
      left_join(ENSG_HGNC_map_table, by = "gene_id") %>%
      select(-gene_id) %>%
      # Only keep protein coding genes
      filter(gene_biotype == "protein_coding") %>%
      select(-gene_biotype) %>%
      select(gene, everything()) %>%
      # Mapping ENSG to HGNC might not be unique
      # Keep the entry with lowest p value
      group_by(gene) %>%
      slice_min(P.Value) %>%
      ungroup %>%
      column_to_rownames("gene")

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
