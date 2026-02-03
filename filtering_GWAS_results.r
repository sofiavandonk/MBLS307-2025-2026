############################################
## Example of exploring results of GWAS
############################################

############################################
# load the necessary libraries
############################################
library(tidyverse)

############################################
# load results of GWAS, results of GWAS will be in the object "plot_frame"
############################################
# In this object, the data are stored in the tabular format, where every row represents GWAS results for a single SNP:
# pval: -log10 p-value
# chr: chromosome where the SNP is located
# pos: position of the SNP on the chromosome, in base pairs
# z-score: scaled effect of the non-reference allele on the phenotype; high positive value - non-reference allele increases the phenotype value
# SNP_index: index of SNP used to simplify post-GWAS analysis and visualization

load("GWAS_sat/frames/ROS_burst_nlp24_AUC90.out")

############################################
# find most significantly associated SNPs
############################################
plot_frame %>% filter(pval == max(pval, na.rm = TRUE))
#pval chr      pos  zscore snp_index
#1 98.11573   3 77349059 2110177    658784
#2 98.11573   3 77349083 2110177    658785
#3 98.11573   3 77348389 2110177    658794
#4 98.11573   3 77348930 2110177    658797

############################################
# define the QTL borders as +/-5 Mb from the most significant SNP
############################################
QTL_start = 77348930 - 5000000
QTL_end = 77349083 + 5000000

############################################
# load gene annotation data for the lettuce genome
# from https://data.4tu.nl/datasets/af1b751d-23a2-4954-ac01-4eb6c68d895b/1
############################################
# you can ignore the parsing errors
gene_annotation <- read_tsv("gene_annotation/Lactuca_sativa_Salinas_V8.annotation_overview.tsv")

############################################
# select genes in the QTL of interest
############################################
# selection based on the chromosome and position of QTL
goi <- gene_annotation %>%
  filter(`sequence ID` == "NC_056625.1" & `start sequence` >= QTL_start & `start sequence` <= QTL_end) %>%
  filter(type == "mRNA")

# number of protein-coding genes in the QTL region
nrow(goi)

# write the gene descriptions into a tab-separated values file; you can open the file later
write_tsv(goi, "gene_annotation/genes-in-QTL.tsv")
