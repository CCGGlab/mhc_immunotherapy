Analysis underlying the results reported in Claeys A and Van den Eynden J. MHC class II genotypes are independent predictors of anti-PD1 immunotherapy response in melanoma, Communications Medicine, 2024.

# Background & Aim

This study aims to find a correlation between MHC class I and class II genotypes, quantified based on their HLA binding properties.

The study is primarily based on data obtained from 144 anti-PD1 treated melanoma patients, as reported by [Liu et al., 2019](https://www.nature.com/articles/s41591-019-0654-5)

Other data sets that were used

- [Hugo et al., 2016](https://doi.org/10.1016/j.cell.2016.02.065): 37 anti-PD1 treated melanoma patients
- [Gide et al., 2019](https://doi.org/10.1016/j.ccell.2019.01.003): 36 anti-PD1 treated melanoma patients
- [Rizvi et al., 2015](https://www.science.org/doi/10.1126/science.aaa1348): 34 anti-PD1 treated lung cancer patients
- [Snyder et al., 2014](https://doi.org/10.1056/NEJMoa1406498): 35 anti-CTLA4 treated melanoma patients
- [Van Allen et al., 2015](https://www.science.org/doi/10.1126/science.aad0095): 110 anti-CTLA4 treated melanoma patients
- [Riaz et al., 2017](https://doi.org/10.1016/j.cell.2017.09.028): 68 anti-PD1 treated melanoma patients

# Environment
  
Analysis was performed in a Conda environment. See **mhc_immunotherapy.yml** for details. **scripts/Rpacks** describes R packages that were installed independently.

# Data processing

## MHC Genotyping

MHC Genotyping was performed using Optitype v1.3.5 (MHC-I) and HLA-HD v1.3 (MHC-II) using non-tumoral, normal blood cell control DNA-Seq data

MHC Genotyping on the Gide data was performed using ArcasHLA v0.2.0 using tumoral RNA-Seq data.

DNA/RNA sequencing data was derived from dbGaP:

- [phs000452.v3.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000452.v3.p1): 
- [phs000980.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000980.v1.p1)
- [phs001041.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phsphs001041.v1.p1)
- [PRJNA307199](https://www.ebi.ac.uk/ena/browser/view/PRJNA307199)
- [PRJNA343789](https://www.ebi.ac.uk/ena/browser/view/PRJNA343789)
- [PRJEB23709](https://www.ebi.ac.uk/ena/browser/view/PRJEB23709)
- [PRJNA359359](https://www.ebi.ac.uk/ena/browser/view/PRJNA359359)
- [PRJNA356761](https://www.ebi.ac.uk/ena/browser/view/PRJNA356761)

Data are stored in a patient x HLA gene data matrix with each cell containing the HLA allele.

## MGBS calculation

For each patient in the dataset the MGBS score was calculated for MHC-I and MHC-II genotypes. This score is the proportion of HLA binders (defined as Kd<500nM) in 1 million random peptides (9-mers for MHC-I and 15-mers for MHC-II)

- MGBS-I (MHC-I): 2 alleles from HLA-A, HLA-B and HLA-C
- MGBS-II (MHC-II): 2 alleles from HLA-DPA1 HLA-DPB1, HLA-DQA1, HLA-DQB1 and HLA-DRB1 alleles. For HLA-DP & HLA HLA-DQ the 4 different heterodimers between A and B chain alleles were considered. 

### Generate fasta files for 1 mill random peptides

```{r}
source("scripts/ICB_get_random_peptides.R")
```

### Run NetMHCPan to determine allele-specific affinities for these peptides 

```{sh}
# MHC-I
scripts/other/run_netmhcpan_v40_random.sh
# MHC-II
scripts/other/run_netmhciipan_v32_random.sh
```

### Merge results in a peptide x HLA allele affinity (Kd) matrix

```{r}
# MHC-I
source("scripts/other/process_netmhcpan_v40_random.R")
# MHC-II
source("scripts/other/process_netmhciipan_v32_random.R")
```

### Calculate MGBS-I and MGBS-II

```{r}
source("scripts/ICB_calculate_MGBS.R")
```

## Clinical Metadata

The following clinical metadata were derived from supplements from the respective studies.
- Survival 
- RECIST 
- TMB: defined as the total number of nonsynonymous mutations
- LOH HLA and B2M
- ploidy, purity, heterogeneity
- LN metastatic state

MHC-I and MHC-II zygosity, defined as the total number of different alleles (maximal 6 and 10 for MHC-I and MHC-II respectively) was calculated from the MHC genotypes.   

```{r}
source("scripts/ICB_data_process.R")
```


## Gene expression quantification

RSEM calculations
```{sh}
"scripts/other/run_rsem.sh"
```

Process RSEM output 
```{r}
source("scripts/other/create_txi_objects.R")
```

Batch correction
```{r}
source("scripts/other/create_icb_expr_rsem_converted-to-tpm_batch.R")
```

# Manuscript

Data used in the manuscript are directly available in *data/MHC_immunotherapy.RData*

## Power analysis

Power analysis to calculate required sample numbers for classical single allele - ICB response association studies:
- Fig. S1

```{r}
source("scripts/manuscript_power_analysis.R")
```

## MGBS survival analysis

Correlation between MGBS scores & survival:
- Fig. 1b,c,e,g: baseline 
- Fig. 1d, FigS2: influence of different MGBS thresholds 
- Fig. 2: stratified for ipilimumab treatment

```{r}
source("scripts/manuscript_ICB_surv.R")
```

## MGBS RECIST analysis 

Correlation between MGBS scores & clinical responses (RECIST criteria)
- Fig. 1h-i: baseline

```{r}
source("scripts/manuscript_RECIST.R")
```

## Multivariate analysis

Multivariate analysis to correlate all available & previously described biomarkers with MGBS
- Fig. 3 & Fig. S6: Logistic regression (response variable = responder)
- Fig. S7: Cox regression (response variable = survival)

```{r}
source("scripts/manuscript_multivariate.R")
```

## MGBS gene expression analysis

GSEA of genes diff. expressed between high & low MGBS scores:
- Fig. 4: comparison of GSEA in TMB with MGBS2
- Fig. S9: complete analysis, stratified for ipilimumab

```{r}
# Perform DGE & GSEA
source("scripts/ICB_DGE.R")
# Analysis as reported in manuscript
source("scripts/manuscript_DGE_analysis.R")
```

## Biomarker correlation analysis

Correlation between available & previously described biomarkers with MGBS
- Fig. S3

```{r}
source("scripts/manuscript_correlation_analysis.R")
```

## Neoantigen burden survival analysis

Alternative approaches to MGBS: absolute & relative neoantigen burden
- Fig. S5

```{r}
source("scripts/manuscript_surv_neoAgB.R")
```

## Survival analysis on validation data

MGBS analysis on other ICB datasets
- Fig. S4

```{r}
source("scripts/manuscript_ICB_validation.R")
```

Multivariate analysis on independent ICB datasets
- Fig. S8

```{r}
source("scripts/manuscript_ICB_validation_multivariate.R")
```

```{r}
source("scripts/manuscript_ICB_validation_compare_variables.R")
```

## Summary table

MGBS scores of all studies used in study:
- Table S1

```{r}
source("scripts/summary table")
```

## Creation of manuscript figures

```{r}
source("scripts/manuscript_create_figs.R")
```

