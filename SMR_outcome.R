library(data.table)
library(dplyr)
df=fread("finngen_R12_F5_DEPRESSIO.gz") # PGC_MDD.tsv
head(df)
df2 <- df[df$rsids=="",]
df3 <- df[!df$rsids=="",]
# FInnGen
df3=data.frame(SNP=df3$rsids,
               A1=df3$alt,
               A2=df3$ref,
               freq=df3$af_alt,
               b=df3$beta,
               se=df3$sebeta,
               p=df3$pval)
df3$n=494164
df3 <-df3%>%distinct(SNP,.keep_all = TRUE)
df4 <- subset(df3, grepl(",", df3$SNP))
df3 <- subset(df3,! grepl(",", df3$SNP))
gc()
fwrite(df3,"R12_depression.ma",sep = "\t",quote = F,row.names = F)
#
df$n = df$NCAS+df$NCON
df1=data.frame(SNP=df$ID,
               A1=df$EA,
               A2=df$NEA,
               freq=df$EAF,
               b=df$BETA,
               se=df$SE,
               p=df$PVAL,
               n=df$n)
fwrite(df1,"PCG_depression.ma",sep = "\t",quote = F,row.names = F)