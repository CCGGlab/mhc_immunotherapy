# Summary table
load("data/MHC_immunotherapy.RData")

# Get scores
ICB_data<- ICB_study[,c("study", "patient", "MGBS1", "MGBS2", "MGBSnorm")]
colnames(ICB_data)<- c("study", "patient", "MGBS-I", "MGBS-II", "MGBS-d")

# STudies included
Liu_2019<- ICB_data[ICB_data$study=="Liu_2019",-1]
Gide_2019<- ICB_data[ICB_data$study=="Gide_2019",-1]
Hugo_2016<- ICB_data[ICB_data$study=="Hugo_2016",-1]
Riaz_2017<- ICB_data[ICB_data$study=="Riaz_2017",-1]
Rizvi_2015 <- ICB_data[ICB_data$study=="Rizvi_2015",-1]
VanAllen_2015<- ICB_data[ICB_data$study=="VanAllen_2015",-1]
Snyder_2014<- ICB_data[ICB_data$study=="Snyder_2014",-1]

# Save
WriteXLS::WriteXLS(c("Liu_2019", "Gide_2019", "Hugo_2016", "Riaz_2017", "Rizvi_2015", "VanAllen_2015", "Snyder_2014"),"results/tables/manuscript_summary_table.xlsx",row.names = F)
