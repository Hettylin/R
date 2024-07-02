#devtools::install_github("MRCIEU/MRInstruments")
setwd("C:/R/workplace/Alzheimer's disease and Sarcopenia")
ieugwasr::api_status()

#install.packages('remotes')
#remotes::install_github('MRCIEU/TwoSampleMR')
library(TwoSampleMR)
#List available GWASs
#ao<-available_outcomes()
#Get instruments
exposure_dat<-extract_instruments("prot-a-114")
#exposure_dat<-extract_instruments(outcomes = "ebi-a-GCST002245")
#Get effects of instruments on outcome
outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,outcomes = 'ieu-b-5118')
#Harmonise the exposure and outcome dataHarmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat,outcome_dat)
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
write.csv(or,"or.csv",row.names = F)

