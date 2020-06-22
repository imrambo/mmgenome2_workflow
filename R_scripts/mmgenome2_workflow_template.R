#!/usr/bin/env Rscript

require(mmgenome2)
require(dplyr)
require(tidyr)

#rm(list=ls())

cols_by_sample <- function(sample, df) {
  ### Return a list of column names for a set of samples
  sample_patt <- paste(sample, collapse = "|")
  sample_vars <- grep(sample_patt, colnames(df), value = TRUE)
  return(sample_vars)
}


setwd("/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/47.concoct.D1105S12metaG")

#Load coverage table
cov_df <- read.table("47.concoct.D1105S12metaG_cov", header = TRUE, sep = "\t")

#Path to genome fasta
genome <- "./47.concoct.D1105S12metaG.fa"

genome_basename <- gsub("\\.fa", "", basename(genome))

#dplyr::arrange(cov_df_long, desc(coverage))

renamed_cols <- gsub("assembly_metaspades_cat_renamed_3kb_|\\.coordsort\\.bam",
                     "", colnames(cov_df))

colnames(cov_df) <- renamed_cols

#Coverage data in long format
cov_df_long <- cov_df %>% tidyr::gather(colnames(cov_df)[2:length(colnames(cov_df))], key = "sample", value = "coverage")

#Boxplot of sample coverage
cov_boxplot <- cov_df_long %>% ggplot(aes(x = sample, y = coverage)) + geom_boxplot() +
  theme(axis.line = element_line(color="black"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(color = "black", angle=65, hjust = 1, size = 7.5),
        axis.text.y=element_text(color="black", size = 10),
        strip.text.x = element_text(size = 11, colour = "black"),
        strip.text.y = element_text(size = 11, colour = "black"))

mm <- mmgenome2::mmload(assembly = genome,
                  coverage = cov_df,
                  verbose = TRUE,
                  kmer_pca = FALSE,
                  kmer_BH_tSNE = FALSE)

#Samples to target
s <- c("D0606U", "D0606S")

svars <- cols_by_sample(sample = s, df = mm)

mmgenome2::mmplot_pairs(mm,
             variables = svars,
             #x = "cov_D0606M02metaG_FD",
             #y = "cov_D0606M12metaG_FD",
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.1,
             size_scale = 0.4,
             textsize = 3)
