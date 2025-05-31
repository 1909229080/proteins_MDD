library(coloc)
library(patchwork)
library(data.table)
library(dplyr)
library(openxlsx)
library(locuscomparer)
# Parameter settings
gene_name <- "ACAA1" 
exp <- "Decode_NEGR1.tsv" #Input file
Tissue <- "Plasma"
c= 1  #chr
rs <- "rs1026566"
kb= 500000 
outcome <- "MDD"
gwas_cases= 357636
gwas_control= 1281936
a=gwas_cases/gwas_control
b=gwas_cases+gwas_control
# Load complete data (pQTLs)
pqtl_data <- fread(exp)
rs_bp <- pqtl_data[pqtl_data$rs_id==rs,]
bp = as.numeric(rs_bp[1,2])
pqtl <- subset (pqtl_data, chromosome == c &
                  base_pair_location >= bp-kb &
                  base_pair_location <= bp+kb)
pqtl$P <- 10^(-pqtl$neg_log_10_p_value)
pqtl=data.frame(SNP=pqtl$rs_id,
                chr=pqtl$chromosome,
                position=pqtl$base_pair_location,
                A1=pqtl$effect_allele,
                A2=pqtl$other_allele,
                beta=pqtl$beta,
                se=pqtl$standard_error,
                p=pqtl$P,
                EAF=pqtl$effect_allele_frequency,
                N=pqtl$n)
# pqtl_data$A1 <- toupper(pqtl_data$A1)
# pqtl_data$A2 <- toupper(pqtl_data$A2)
pqtl$MAF <- pmin(pqtl$EAF, 1 - pqtl$EAF)
pqtl$gene <- gene_name
# Another input format, which can be ignored
pqtl_data <- fread(exp)
head(pqtl_data)
rs_bp <- pqtl_data[pqtl_data$SNP==rs,]
bp = as.numeric(rs_bp[1,3])
pqtl <- subset (pqtl_data, chr == c &
                  pos >= bp-kb &
                  pos <= bp+kb)
pqtl=data.frame(SNP=pqtl$SNP,
                chr=pqtl$chr,
                position=pqtl$pos,
                A1=pqtl$effect_allele.outcome,
                A2=pqtl$other_allele.outcome,
                beta=pqtl$beta.outcome,
                se=pqtl$se.outcome,
                p=pqtl$pval.outcome,
                EAF=pqtl$eaf.outcome,
                N=pqtl$samplesize.outcome)
pqtl$MAF <- pmin(pqtl$EAF, 1 - pqtl$EAF)
pqtl$gene <- gene_name
# Match outcome data
gwas_data <- fread("PCG_depression.tsv")
gwas_data=data.frame(SNP=gwas_data$ID,chr=gwas_data$`#CHROM`,
                     position=gwas_data$POS,A1=gwas_data$EA,
                     A2=gwas_data$NEA,beta=gwas_data$BETA,
                     se=gwas_data$SE,p=gwas_data$PVAL)
gwas <- subset(gwas_data,chr == c )
# Convert GRch37 to GRch38; if unnecessary, it can be ignored.
#BiocManager::install("rtracklayer") # Install rtracklayer
library(rtracklayer)
# Define the path
chain_file_path <- "D:/R/Rlibrary/hg19ToHg38.over.chain"
chain <- import.chain(chain_file_path)
# Create a GRanges object
gr <- GRanges(
  seqnames = paste0("chr", gwas$chr),
  ranges = IRanges(start = gwas$position, width = 1),
  strand = "*",
  originalIndex = seq_len(nrow(gwas)))   
# Perform coordinate transformation
lifted_gr <- liftOver(gr, chain)
lifted_gr <- unlist(lifted_gr)
lifted_df <- as.data.frame(lifted_gr)
lifted_df <- lifted_df[, c("originalIndex", "seqnames", "start")]
gwas$originalIndex <- seq_len(nrow(gwas))
gwas <- merge(gwas, lifted_df, 
                   by = "originalIndex", 
                   all.x = TRUE,
                   suffixes = c("_orig", "_lifted"))
gwas$lift_success <- !is.na(gwas$start)
gwas$new_chr <- gsub("^chr", "", gwas$seqnames)
gwas$new_pos <- gwas$start
# Check conversion success rate
cat("success_rate:", mean(gwas$lift_success, na.rm = TRUE) * 100, "%\n")
gwas <- gwas[,-c(1,3,4,10:12)]
gwas <- gwas %>% rename_at(vars(7:8), ~ c("chr", "position"))
# Intersect exposure and outcome data
merged_data <- merge(gwas, pqtl, 
                     by = c("SNP", "chr","position"),
                     suffixes = c("_gwas", "_pqtl"))
gc()
merged_data = merged_data %>% filter((A1_gwas==A1_pqtl&A2_gwas==A2_pqtl)|
                                       (A1_gwas==A2_pqtl&A2_gwas==A1_pqtl))
# Check and reverse the beta direction to match the GWAS allele
flip <- merged_data$A1_gwas != merged_data$A1_pqtl
merged_data$beta_pqtl[flip] <- -merged_data$beta_pqtl[flip]
region_data <-merged_data%>%distinct(SNP,.keep_all = TRUE)  
# Co-localization analysis (must include beta and varbeta)
dataset_gwas <- list(beta = region_data$beta_gwas,
                     varbeta = region_data$se_gwas^2,
                     type = "cc",     # Case-control study categorical variables
                     s = a,         # Case ratio (adjusted based on actual data)
                     N = b,       # Total sample size
                     snp = region_data$SNP,
                     p_gwas = region_data$p_gwas,
                     position = region_data$position)
dataset_pqtl <- list(beta = region_data$beta_pqtl,
                     varbeta = region_data$se_pqtl^2,
                     type = "quant", 
                     N = region_data$N,         
                     snp = region_data$SNP,
                     p_pqtl = region_data$p_pqtl,
                     position = region_data$position,
                     MAF = region_data$MAF,
                     gene = region_data$gene[1]  
                     )
coloc_result <- coloc.abf(dataset1 = dataset_gwas,dataset2 = dataset_pqtl,
                         p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
print(coloc_result$summary)
# Output results
res <-as.data.frame(coloc_result[["summary"]])
new_res <- as.data.frame(t(res))
dt <- data.frame(Tissue,gene_name)
dt$outcome <- outcome
dt$kb =500
dt$p12 = 1e-05
res_500kb <- cbind(dt,new_res)
write.xlsx(res_500kb,"res_Colocalization.xlsx")
# Visualization
data <- as.data.frame(dataset_pqtl)
values <- data$p_pqtl[data$p_pqtl > 0]  
min_values <- min(values, na.rm = TRUE)  
data$p_pqtl[data$p_pqtl == 0] <- min_values
data2 <-as.data.frame(dataset_gwas)
write.table(data, "pQTL.txt",row.names = F,sep = "\t",quote = F)
write.table(data2, "GWAS.txt",row.names = F,sep = "\t",quote = F)
# Path of the input file
pqtl_fn="D:/R/MDD/pQTL.txt"
gwas_fn="D:/R/MDD/GWAS.txt"
P <-locuscompare(in_fn1=pqtl_fn, in_fn2=gwas_fn, 
                 title1=gene_name,
                 title2=outcome, 
                 marker_col1="snp", pval_col1= "p_pqtl", 
                 marker_col2="snp", pval_col2="p_gwas")
print(P)
file_name <- paste0(Tissue,"-","pQTL-",gene_name,"-",outcome,".pdf")
library(ggplot2)
pdf(file=file_name, width=7.5, height=8)
plot(P)
dev.off()
# Manhattan plot function (displaying significance thresholds)
plot_manhattan <- function(dat, pval_col, title) {
  ggplot(dat, aes(x = position, y = -log10(.dat[[pval_col]]))) +
    geom_point(alpha = 0.6) +geom_hline(yintercept = -log10(5e-8),
                                        color = "red", linetype = "dashed") +
    labs(title = title, x = "Position", y = "-log10(p)") +
    theme_minimal()
}
p_gwas <- plot_manhattan(region_data, "p_gwas", "MDD")#"GWAS Associations"
p_pqtl <- plot_manhattan(region_data, "p_pqtl", 
                         paste("pQTL for", region_data$gene[1]))
combined_plot <- p_gwas / p_pqtl +
  plot_annotation(title = "Co-localization Analysis",tag_levels = "A")
print(combined_plot)
file_name2 <- paste0(Tissue,"-","pQTL-",gene_name,"-",outcome,".tif")
ggsave(file_name2 , combined_plot, width = 10, height = 7)
