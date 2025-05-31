#Packages
library(TwoSampleMR)
library(openxlsx)
library(dplyr)
library(gridExtra) 
library(data.table)
library(forestplot) 
library(scales) 
library(tidyverse) 
library(RVAideMemoire) 
library(rstatix) 
# Exposures
csf <- read.xlsx("CSF_protein_cis.xlsx")
plasma <- read.xlsx("Plasma_protein_cis.xlsx")
# Outcomes:"pgc-mdd2025-no23andMe-noUKBB-eur"
out <- read_outcome_data(
  snps = unique(c(csf$SNP,plasma$SNP)),
  filename = "PGC_MDD.tsv",
  sep = "\t",
  chr_col = "#CHROM",               
  pos_col = "POS",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "PVAL",
  ncase_col = "NCAS",
  ncontrol_col = "NCON",)
out$id.outcome <-"PGC"
out$outcome <-"MDD"
write.xlsx(out,"MDD_SNPs.xlsx")
# Harmonise the exposure and outcome data
csf_out <- harmonise_data(csf, out)
plasma_out <- harmonise_data(plasma, out)
# Custom analysis function
mr_combined <- function(dat , prop_var_explained = T)
{
  dat$R2 <- 2* dat$eaf.exposure*(1-dat$eaf.exposure)*dat$beta.exposure^2
  dat$F <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)
  dat <-dat[dat$F>20,] #Select instrumental variables greater than 20
  a <-head(dat$Tissue,1)
  b <- "_harmonise_data.xlsx"
  write.xlsx(dat,paste0(a,b)) #Output the data
  dat <- dat[dat$mr_keep == TRUE, ]#Exclude SNPs that did not match
  mr_res <- mr(dat,method_list = c("mr_ivw","mr_wald_ratio"))
  pve <- dat %>% #Calculate the sum of R-squared values for each exposure
    dplyr::select(id.exposure, R2) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum(R2))
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  return(mr_res)
}
# Extract sample information
csf_sample <- unique(csf[,c("id.exposure","samplesize.exposure","Tissue","Author")])
csf_sample <-csf_sample%>%distinct(id.exposure,.keep_all = TRUE)  
plasma_sample <- unique(plasma[, c("id.exposure","samplesize.exposure","Tissue","Author")])
plasma_sample <-plasma_sample%>%distinct(id.exposure,.keep_all = TRUE)  
# Calculate the effects of CSF proteins on MDD by custom analysis function
mr_csf <- mr_combined(csf_out, prop_var_explained = T)
mr_csf <- merge(mr_csf, csf_sample,by="id.exposure",all = TRUE)
mr_csf <- na.omit(mr_csf)
mr_csf$F_statistics <- mr_csf$pve*
  (mr_csf$samplesize.exposure - mr_csf$nsnp-1) / 
  (1 - mr_csf$pve)#Calculate the F-statistic
write.xlsx(mr_csf,"mr_res_csf_MDD.xlsx")
# Calculate the effects of plasma proteins on MDD by custom analysis function
mr_plasma <- mr_combined(plasma_out, prop_var_explained = T)
mr_plasma <- merge(mr_plasma, plasma_sample,by="id.exposure",all = TRUE)
mr_plasma <- na.omit(mr_plasma)
mr_plasma$F_statistics <- mr_plasma$pve*
  (mr_plasma$samplesize.exposure - mr_plasma$nsnp-1) / 
  (1 - mr_plasma$pve)
write.xlsx(mr_plasma,"mr_res_plasma_MDD.xlsx")
# Calculate the significance threshold and filter out significant associations
csf_list <-  mr_csf %>% group_split(id.exposure)
plasma_list <-  mr_plasma %>% group_split(exposure)
m = length(csf_list) + length(plasma_list)
P = 0.05/m
mr_csf_filter <- mr_csf[mr_csf$pval <  P,]
mr_plasma_filter <- mr_plasma[mr_plasma$pval < P,]
mr_res <- rbind(mr_csf_filter, mr_plasma_filter)
OR <-generate_odds_ratios(mr_res)
OR$or <-sprintf("%.2f", OR$or)
OR$or_lci95 <-sprintf("%.2f", OR$or_lci95)
OR$or_uci95 <-sprintf("%.2f", OR$or_uci95)
write.xlsx(OR,"Significant_results.xlsx")
# Create a function to generate volcano plots
volcano_plot <- function(.data, 
                         number_comparasion = m,
                         title = "(A)",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         col_se = "se",
                         legend.position = "none") {
  p_thershold <- 0.05 / number_comparasion
# Determine the color of the points based on the value of F_statistics.
# Use the mutate function of dplyr along with the case_when statement for conditional judgments.
# Classify each observation in the dataframe into different color groups based on the value of F_statistics.
  .data <-.data %>%
    mutate(color_group = case_when(
      F_statistics < 50 ~ "blue",
      F_statistics >= 50 & F_statistics <= 100 ~ "orange",
      F_statistics > 100 ~ "red"
    ))
# Preprocess the data, including renaming columns and creating variables required for plotting.
# Use the rename function of dplyr to rename columns according to specified parameters.
# Use the mutate function to create variables x (corresponding to beta values, i.e., ln(OR)) and y (corresponding to -log10(pval)).
# Set the label variable based on the p-value threshold; retain the original annotation if the p-value is below the threshold, otherwise set it to NA.
  .data <-.data %>%
    rename(beta :=!!col_beta,
           size :=!!col_size,
           label :=!!col_label) %>% 
    mutate(x = beta,
           y = -log10(pval),
           label = ifelse(pval < p_thershold, label, NA))
# Calculate confidence intervals but only for points where pval is less than p_threshold.
  .data <-.data %>%
    mutate(lower_ci = ifelse(pval < p_thershold, beta - 1.96 *!!sym(col_se), NA),
           upper_ci = ifelse(pval < p_thershold, beta + 1.96 *!!sym(col_se), NA))
# Begin creating a volcano plot using the ggplot2 library
  p <-.data %>% 
    ggplot(aes(x = x, y = y)) +
# Draw confidence interval lines, but only for points where pval is less than p_threshold.
    geom_errorbarh(aes(xmin = ifelse(pval < p_thershold, lower_ci, NA), 
                       xmax = ifelse(pval < p_thershold, upper_ci, NA), 
                       y = y), height = 0.07, color = "grey") +
# The size of the markers is determined by the size variable, their color is based on the color_group variable, and the transparency (alpha) is set to 0.5.
    geom_point(aes(size = size, color = color_group), alpha =0.7 ) +
# Add a dashed line perpendicular to the x-axis with an x-intercept set at 0 to delineate different regions.
    geom_vline(xintercept = 0, linetype = 2) +
# Add a dashed line perpendicular to the y-axis, with the y-intercept set to -log10(p_threshold), to delineate different regions.
    geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
# Set the theme style of the plot to classic style
    theme_classic() +
# Further configure parameters related to the theme
# Hide panel grid
    theme(panel_grid = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.position = legend.position) +
    labs(x = "ln(OR)", 
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title) +
    scale_size(name = "PVE",
               breaks = c(0.15 * 1:3),
               range = c(2.5, 12.5)) +
    scale_color_manual(values = c("blue" = "blue", "orange" = "orange", 
                                  "red" = "red"),
                       name = "F-statistics", breaks = c("blue", "orange", "red"),
                       labels = c("20 - 50", 
                                  "50 - 100",
                                  "> 100"),
                       guide = guide_legend(override.aes = list(size = 6))) +
    ggrepel::geom_label_repel(aes(label = label), size = 5)
  plot(p)
}

gridExtra::grid.arrange(
  mr_plasma %>% 
    filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
    volcano_plot(legend.position = c(0.1,0.8),title = "Plasma"),#number_comparasion = m,
  mr_csf %>% 
    filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
    volcano_plot(title = "CSF"),
  ncol = 2,
  nrow = 1)
# Steiger filtering
res <- read.xlsx("Significant_results.xlsx")
dat_plasma <- read.xlsx("Plasma_harmonise_data.xlsx")
dat_CSF <- read.xlsx("CSF_harmonise_data.xlsx")
dat <- rbind(dat_plasma,dat_CSF)
dat_res <- dat[dat$Tissue %in% res$Tissue & dat$exposure %in% res$exposure,]
dat_res1 <- dat_res[dat_res$Tissue=="Plasma",]
dat_res2 <- dat_res[dat_res$Tissue=="CSF",]
filtering1 <- steiger_filtering(dat_res1)
filtering2 <- steiger_filtering(dat_res2)
filtering <- rbind(filtering1,filtering2)
write.xlsx(filtering,"Steiger_filtering_res.xlsx")
# External validation 1
res <- read.xlsx("Significant_results.xlsx")
res_Plasma <- res[res$Tissue=="Plasma",]
res_CSF <- res[res$Tissue=="CSF",]
plasma <- read.xlsx("Plasma_protein_cis.xlsx")
plasma2 <- plasma[plasma$exposure %in% res_Plasma$exposure,]
plasma2$exposure <- paste0(plasma2$Tissue,"_",plasma2$exposure)
plasma2$id.exposure <- paste0(plasma2$Tissue,"_",plasma2$id.exposure)
csf <- read.xlsx("CSF_protein_cis.xlsx")
csf2 <- csf[csf$exposure %in% res_CSF$exposure,]
csf2$exposure <- paste0(csf2$Tissue,"_",csf2$exposure)
csf2$id.exposure <- paste0(csf2$Tissue,"_",csf2$id.exposure)
exp <- rbind(plasma2,csf2)
write.xlsx(exp,"Prominent_protein_cis.xlsx")
#FinnGen Cohort
out2 <- read_outcome_data(
  snps = unique(c(exp$SNP)),
  filename = "finngen_R12_F5_DEPRESSIO.gz",
  sep = "\t",
  chr_col = "#chrom",               
  pos_col = "pos",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")
out2$id.outcome <-"Finngen_R12"
out2$outcome <-"Depression"
write.xlsx(out2,"FinnGen_SNPs_1.xlsx")
FinGen <- harmonise_data(exp,out2)
FinGen$samplesize.outcome=494164
FinGen_res <- mr(FinGen)
FinGen_res <-generate_odds_ratios(FinGen_res)
FinGen_res[,12:14] <- round(FinGen_res[,12:14],2)
FinGen_res <- FinGen_res[,-c(7,8,10,11)]
write.xlsx(FinGen_res,"External_validation_1.xlsx")
FinGen_filtering <- steiger_filtering(FinGen)
write.xlsx(FinGen_filtering,"External_validation_1_steiger_filtering.xlsx")
# External validation 2_PGC_Cohort
rm(list=ls())
exp <- read.xlsx("External_validation_2_proteins_cis.xlsx")
PGC <- read_outcome_data(
  snps = unique(c(exp$SNP)),
  filename = "PGC_MDD.tsv",
  sep = "\t",
  chr_col = "#CHROM",               
  pos_col = "POS",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "PVAL",
  ncase_col = "NCAS",
  ncontrol_col = "NCON",)
PGC$id.outcome <-"PGC"
PGC$outcome <-"MDD"
write.xlsx(PGC,"MDD_SNPs_2.xlsx")
pqtl_MDD <- harmonise_data(exp, PGC)
res1 <- mr(pqtl_MDD,method_list = c("mr_wald_ratio","mr_ivw",
                                    "mr_egger_regression",
                                    "mr_weighted_median"))
res1 <- generate_odds_ratios(res1)
res1[,12:14] <- round(res1[,12:14],2)
res1 <- res1[,-c(7,8,10,11)]
res1$'OR(95%CI)' <- paste(res1$or,"(",res1$or,"-",res1$or_uci95,")")
write.xlsx(res1,"External_validation_2_res_PGC.xlsx")
filtering1 <- steiger_filtering(pqtl_MDD)
write.xlsx(filtering1,"External_validation_2_steiger_filtering1.xlsx")
pleiotropy <-mr_pleiotropy_test(pqtl_MDD)
write.xlsx(pleiotropy,"External_validation_2_pleiotropy1.xlsx")
# External validation 2_FinnGen_Cohort
FinGen <- read_outcome_data(
  snps = unique(c(exp$SNP)),
  filename = "finngen_R12_F5_DEPRESSIO.gz",
  sep = "\t",
  chr_col = "#chrom",               
  pos_col = "pos",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")
FinGen$id.outcome <-"Finngen_R12"
FinGen$outcome <-"Depression"
write.xlsx(FinGen,"FinnGen_SNPs_2.xlsx")
FinGen$samplesize.outcome=494164
pqtl_D <- harmonise_data(exp,FinGen)
res2 <- mr(pqtl_D,method_list = c("mr_wald_ratio","mr_ivw",
                                  "mr_egger_regression",
                                  "mr_weighted_median"))
res2 <-generate_odds_ratios(res2)
res2[,12:14] <- round(res2[,12:14],2)
res2 <- res2[,-c(7,8,10,11)]
res2$'OR(95%CI)' <- paste(res2$or,"(",res2$or,"-",res2$or_uci95,")")
write.xlsx(res2,"External_validation_2_res_FinnGen.xlsx")
filtering2 <- steiger_filtering(pqtl_D)
write.xlsx(filtering2,"External_validation_2_steiger_filtering2.xlsx")
pleiotropy2 <-mr_pleiotropy_test(pqtl_D)
write.xlsx(pleiotropy2,"External_validation_2_pleiotropy2.xlsx")
