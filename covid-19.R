#devtools::install_github("MRCIEU/MRInstruments")
ieugwasr::api_status()
install.packages('openxlsx')
setwd("C:\\R\\workplace\\cFDR")

library(data.table)
library(TwoSampleMR)
library(openxlsx)
library(friendly2MR)

covid_19<-fread('finngen_R10_F5_ALZHDEMENT',head=T)
covid_19$phenotype<-'Covid'

exposure_dat <- format_data(
  dat=covid_19,
  type = "exposure",
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  eaf_col = "af_alt"
  #samplesize_col = "all_meta_N"
)
exposure_dat$id.exposure<-'finngen_R10_G6_ALZHEIMER'
#pval<e-5 snp
exposure_dat <-subset(exposure_dat,pval.exposure < 0.00000005)
#remove LD
exposure_dat<-clump_data(exposure_dat)



#neuroticism
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006368")

# #Anxiety
# outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="finn-b-F5_PHOBANX")
# 
# #depression
# outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-102")
# 
# #Post-traumatic stress disorder
# outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="finn-b-F5_PTSD")


dat <- harmonise_data(exposure_dat,outcome_dat)
#res_single <- mr_singlesnp(dat)

#MR analysis
res <- mr(dat)
or<-generate_odds_ratios(res);or

#sensitivity analyses 
#异质性
he<-mr_heterogeneity(dat);he
#多效性
ple<-mr_pleiotropy_test(dat);ple

write.csv(he,"he.csv",row.names = F)
write.csv(ple,"ple.csv",row.names = F)
write.csv(or,"or_covid_19_gut_microbiota.csv",row.names = F)


#reverse
#Get instruments
exposure_dat <- extract_instruments("ebi-a-GCST005348")

covid_19<-fread('finngen_R10_G6_ALZHEIMER',head=T)
covid_19$phenotype<-'Alzheimer disease'

outcome_dat <- format_data(
  dat=covid_19,
  type = "outcome",
  snps=exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  samplesize_col = "n",
  #  eaf_col = "all_meta_AF"
)
outcome_dat$id.outcome<-'finngen_R10_G6_ALZHEIMER'


#Remove invalid SNPs
dat <- harmonise_data(exposure_dat,outcome_dat)
#res_single <- mr_singlesnp(dat)

#MR analysis
res <- mr(dat)
or<-generate_odds_ratios(res);or
#sensitivity analyses 
#异质性
he<-mr_heterogeneity(dat);he
#多效性
ple<-mr_pleiotropy_test(dat);ple
write.csv(he,"he.csv",row.names = F)
write.csv(ple,"ple.csv",row.names = F)
write.csv(or,"or_covid_19_gut_microbiota.csv",row.names = F)




