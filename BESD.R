library(data.table)
library(dplyr)
exp <- "GCST90426110_ITIH4.tsv.gz"
pqtl_data <- fread(exp)
pqtl_data$P <- 10^(-pqtl_data$neg_log_10_p_value)
pqtl_data=data.frame(SNP=pqtl_data$rs_id,
                     Chr=pqtl_data$chromosome,
                     BP=pqtl_data$base_pair_location,
                     A1=pqtl_data$effect_allele,
                     A2=pqtl_data$other_allele,
                     Freq=pqtl_data$effect_allele_frequency,
                     Probe="4811_33",
                     Probe_Chr=3,
                     Probe_bp=52821825,
                     Gene="ITIH4",
                     Orientation="-",
                     b=pqtl_data$beta,
                     se=pqtl_data$standard_error,
                     p=pqtl_data$P)
gc()
pqtl_data <- pqtl_data[!pqtl_data$SNP=="",]
write.table(pqtl_data, "ITIH4.txt",row.names = F,sep = "\t",quote = F)