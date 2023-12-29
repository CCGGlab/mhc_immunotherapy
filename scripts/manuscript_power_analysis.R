# Power analysis to detect survival differences for different MHC alleles

# Load allele frequency data
AF<- readRDS("downloads/AFN/consensus_ggroup.rds")
# install.packages("powerSurvEpi")
library(powerSurvEpi)

# Create plot with allele frequency distribution for MHC-I and MHC-II alleles
#############################################################################

# Caucasian american population

# Data frame
f<- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
MHC<- c("MHC-I","MHC-II")
MHC_AF<- matrix(NA,length(f),length(MHC), dimnames = list(f,MHC))
for(MHC_tmp in MHC){
  if(MHC_tmp=="MHC-I") AF_tmp<- AF[AF$gene%in%c("A","B","C"),]
  if(MHC_tmp=="MHC-II") AF_tmp<- AF[AF$gene%in%c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1"),]
  for(f_tmp in f){
    MHC_AF[as.character(f_tmp),MHC_tmp]<- mean(AF_tmp$`Caucasian American`<f_tmp,na.rm=T)
  }
}

MHC_AF
#        MHC-I    MHC-II
# 0.001 0.8797468 0.7655860
# 0.003 0.9103376 0.8104738
# 0.01  0.9472574 0.8703242
# 0.03  0.9683544 0.9201995
# 0.1   0.9915612 0.9551122
# 0.3   1.0000000 0.9950125

for(i in nrow(MHC_AF):2){
  MHC_AF[i,]<- MHC_AF[i,]-MHC_AF[i-1,]
}

# Plot
MHC_AF_df<- melt(MHC_AF)
colnames(MHC_AF_df)<- c("AF","MHC","prop")
MHC_AF_df$AF_name<- rep(c("<0.1%","0.1-0.3%","0.3-1%","1%-3%","3-10%",">10%"),2)
MHC_AF_df$AF_name<- factor(MHC_AF_df$AF_name, levels=c("<0.1%","0.1-0.3%","0.3-1%","1%-3%","3-10%",">10%"))
p_AF<- ggplot(MHC_AF_df,aes(x=AF_name,y=prop,fill=MHC)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Allele Frequency") +
  scale_y_continuous("Proportion", limits = c(0,1), expand = c(0, 0)) +
  theme_classic2() +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.title = element_blank(),
    legend.text = element_text(size=8)
  )

# Power analysis
#################

# Create dataframe
f<- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
rr<- c(0.4, 0.7,1.2, 1.5, 2, 3)
f_ss_df<- matrix(NA,length(f),length(rr), dimnames = list(f,rr))

for(f in rownames(f_ss_df)){
  for(rr in colnames(f_ss_df)){
    res<- ssizeCT.default(power = 0.8, 
                          k = as.numeric(f)/(1-as.numeric(f)), #E/C, 
                          pE = 0.5, 
                          pC = 0.5, 
                          RR = as.numeric(rr), 
                          alpha = 0.05)
    f_ss_df[f,rr]<- sum(res)
  }
}
f_ss_df<- reshape2::melt(f_ss_df, )
colnames(f_ss_df)<- c("AF","HR","n")
f_ss_df$HR<- as.character(f_ss_df$HR)
f_ss_df$AF_name<- rep(c("0.1%","0.3%","1%","3%","10%","30%"),length(rr))
f_ss_df$AF_name<- factor(f_ss_df$AF_name, levels=c("0.1%","0.3%","1%","3%","10%","30%"))

# Create plot
p_ss<- ggplot(data=f_ss_df, aes(x=AF_name, y=n, group=HR, colour=HR)) +
  geom_line()+
  geom_point() +
  xlab("Allele Frequency") +
  scale_y_continuous("Sample Size", trans='log10',limits = c(10,1000000), expand = c(0, 0), breaks = c(10,100,1000,10000,100000,1000000)) +
  theme_classic2() +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    axis.title.x = element_blank(),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size=8)
  )
  # geom_hline(yintercept = 144)
p_ss

# Save
#######
p_AF_power<- list(SS=p_ss, AF=p_AF)
saveRDS(p_AF_power, file="results/data/manuscript_AF_power_analysis.rds")
