library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(data.table)
library(readr)
set.seed(1234)

rm(list = ls())
###################prepare proteins GWAS data#######################
folder_path <- "C:/Users/JiefeiNiu/Desktop/sun/"
#folder_path <- "/home/ame/jiefei.niu/tools/project/pQTL/" #put in server

# extract files' name
file_names <- list.files(folder_path)[1]
#file_names <- c(file_names[5])

result_list <- list()
for (file in file_names) {
  file_path <- paste0("C:/Users/JiefeiNiu/Desktop/sun/", file)
  data <- fread(file_path)
  #p_value < 5e-08 or logpvalue >8
  ex <- subset(data, p_value < 5e-08)
  # add files name as a column
  ex$FileName <- gsub(".tsv", "", file)
  # put results in list
  result_list[[file]] <- ex
}

#combine all res together
colnames(data)
head(data)
result <- rbindlist(result_list)
colnames(result)
#save it
#write.csv(result, file = "C:/project/pro/onlypro/result/mr/result.csv")
#write.csv(result, file = "C:/project/pro/onlypro/result/mr/result_LGALS3BP.csv")


###########################proteins lead to BMI###################################
#Ferkingstad, E. et al. Large-scale integration of the plasma proteome with genetics and disease. 
result <- read.csv("C:/project/pro/onlypro-final/result/mr/result.csv")
#check data
colnames(result)
table(result$FileName)
head(result)
#include 9 PROTEINS: A2M, APCS, ADIPOQ, AFM, APOD, AZGP1, CRP, S100A9, SERPINA6
result <- subset(result, genename != "SERPINA6") #SERPINA6 cann't clumped
exp_dat <- format_data(
  result,
  type = "exposure",
  header = TRUE,
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  pval_col = "Pval",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  ncase_col = "N",
  id_col = "genename")

#APOF "prot-a-133"; APOM "prot-a-136"; CFH "prot-a-520"; A2M 'prot-a-1781';
#SERPINF1 "prot-a-2698"; 
#put the id one by one and get the results
exp_dat <- extract_instruments(outcomes = 'prot-a-2698')
snp_exp_dat <- clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.001,
                          clump_p1 = 1, clump_p2 = 1, pop = "EUR")
# outcomes data
out_dat <- extract_outcome_data(
  snps = snp_exp_dat$SNP,
  outcomes = 'ebi-a-GCST90029007'
)

#harmonise data
dat <- harmonise_data(
  exposure_dat = snp_exp_dat, 
  outcome_dat = out_dat)
#write.csv(dat, file = "C:/project/pro/onlypro-final/result/mr/result_snps.csv")
res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw")) #"mr_egger_regression",
res
#Sensitivity analyses 
#Heterogeneity statistics
mr_heterogeneity(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_heterogeneity(dat, method_list = c("mr_ivw"))
#可以挑出离群值
#run_mr_presso(dat, NbDistribution = 3000)
#Horizontal pleiotropy
mr_pleiotropy_test(dat)

#PLOT
res_single=mr_singlesnp(dat)
mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))

#write.csv(snp_exp_dat, file = "C:/project/pro/onlypro/result/mr/snp_exp_dat.csv")
#write.csv(res, file = "C:/project/pro/onlypro/result/mr/res.csv")

#from LGALS3BP to BMI
result <- read.csv("C:/project/pro/onlypro-final/result/mr/result_LGALS3BP.csv")
exp_dat <- format_data(
  result,
  type = "exposure",
  header = TRUE,
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  ncase_col = "n")

snp_exp_dat <- clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.001,
                          clump_p1 = 1, clump_p2 = 1, pop = "EUR")
#write.csv(snp_exp_dat, file = "C:/project/pro/onlypro/result/mr/snp_exp_dat.csv")
# outcomes data
out_dat <- extract_outcome_data(
  snps = snp_exp_dat$SNP,
  outcomes = 'ebi-a-GCST90029007'
)

#harmonise data
dat <- harmonise_data(
  exposure_dat = snp_exp_dat, 
  outcome_dat = out_dat)
#write.csv(dat, file = "C:/project/pro/onlypro-final/result/mr/result_snps.csv")
res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
res
#Sensitivity analyses 
#Heterogeneity statistics
mr_heterogeneity(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_heterogeneity(dat, method_list = c("mr_ivw"))
#可以挑出离群值
#run_mr_presso(dat, NbDistribution = 3000)
#Horizontal pleiotropy
mr_pleiotropy_test(dat)



###########################BMI lead to proteins###################################
#exposure
exp_dat1 <-extract_instruments(outcomes = 'ebi-a-GCST90029007')
snp_exp_dat1 <- clump_data(exp_dat1, clump_kb = 10000, clump_r2 = 0.001,
                           clump_p1 = 1, clump_p2 = 1, pop = "EUR")
#outcome
#include 9 PROTEINS: A2M, APCS, ADIPOQ, AFM, APOD, AZGP1, CRP, S100A9, SERPINA6
out_dat1 <- read_outcome_data(
  snps = snp_exp_dat1$SNP,
  filename = "C:/project/pro/onlypro-final/result/mr/result.csv",
  sep = ",",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  pval_col = "Pval",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  ncase_col = "N",
  id_col = "genename"
)

#harmonise data
dat <- harmonise_data(
  exposure_dat = snp_exp_dat1, 
  outcome_dat = out_dat1)

res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
res
#Sensitivity analyses 
#Heterogeneity statistics
mr_heterogeneity(dat, method_list = c("mr_ivw"))
#可以挑出离群值
#run_mr_presso(dat, NbDistribution = 3000)
#Horizontal pleiotropy
mr_pleiotropy_test(dat)

#PLOT
res_single=mr_singlesnp(dat)
mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))


#directly find GWAS data from GWAS catalog
#APOF "prot-a-133"; APOM "prot-a-136"; CFH "prot-a-520"; A2M 'prot-a-1781'
out_dat1 <- extract_outcome_data(
  snps = snp_exp_dat1$SNP,
  outcomes = 'prot-a-520'
)

#harmonise data
dat <- harmonise_data(
  exposure_dat = snp_exp_dat1, 
  outcome_dat = out_dat1)

res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
res
#Sensitivity analyses 
#Heterogeneity statistics
mr_heterogeneity(dat, method_list = c("mr_ivw"))
#可以挑出离群值
#run_mr_presso(dat, NbDistribution = 3000)
#Horizontal pleiotropy
mr_pleiotropy_test(dat)

#PLOT
res_single=mr_singlesnp(dat)
mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
