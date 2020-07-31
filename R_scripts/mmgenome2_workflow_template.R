#!/usr/bin/env Rscript

library(mmgenome2)
library(tidyverse)

#rm(list=ls())

cols_by_sample <- function(sample, df) {
  ### Return a list of column names for a set of samples
  sample_patt <- paste(sample, collapse = "|")
  sample_vars <- grep(sample_patt, colnames(df), value = TRUE)
  return(sample_vars)
}

gg_boxplot <- function(df, x_vec, y_vec) {
  #Create a box and whisker plot from long-format data
  df %>% ggplot(aes(x = x_vec, y = y_vec)) + geom_boxplot() +
    theme(axis.line = element_line(color="black"),
          axis.title=element_text(size=12),
          axis.text.x=element_text(color = "black", angle=65, hjust = 1, size = 7.5),
          axis.text.y=element_text(color="black", size = 10),
          strip.text.x = element_text(size = 11, colour = "black"),
          strip.text.y = element_text(size = 11, colour = "black"))
}

cov_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/mmgenome2/cov_files"
genome_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/das_tool/DASTool_Run2_concoct-maxbin2-metabat2-vamb_diamond_DASTool_bins"

#Path to genome coverage file
cov_file <- paste(cov_dir,
                  '1.concoct.D0819M02metaG/1.concoct.D0819M02metaG_cov',
                  sep = "/")

cov_df <- read.table(cov_file, header = TRUE, sep = "\t")

cov_cols_renamed <- gsub("assembly_metaspades_cat_renamed_3kb_|\\.coordsort\\.bam",
                         "", colnames(cov_df))

colnames(cov_df) <- cov_cols_renamed

#Coverage data in long format
cov_df_long <- cov_df %>%
  tidyr::gather(colnames(cov_df)[2:length(colnames(cov_df))],
                key = "sample", value = "coverage")

cov_quartiles <- cov_df_long %>% group_by(sample) %>%
  summarise(tibble::enframe(quantile(coverage, c(0.25, 0.5, 0.75)),
                            "quantile", "q_coverage"))

cov_df_long %>% left_join(cov_quartiles) %>% head()

upper_outlier <- upper_quartile + (IQR * 1.5)

cov_lower_quartile <- cov_quartiles[[2]]
cov_upper_quartile <- cov_quartiles[[5]]


#Boxplot of sample coverage
cov_boxplot <- cov_df_long %>% ggplot(aes(x = sample, y = coverage)) + geom_boxplot() +
  theme(axis.line = element_line(color="black"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(color = "black", angle=65, hjust = 1, size = 7.5),
        axis.text.y=element_text(color="black", size = 10),
        strip.text.x = element_text(size = 11, colour = "black"),
        strip.text.y = element_text(size = 11, colour = "black"))

cov_boxplot
#dplyr::arrange(cov_df_long, desc(coverage))
#-----------------------------------------------------------------------------
#Path to genome fasta
genome <- paste(genome_dir, "1.concoct.D0819M02metaG.fa", sep = "/")

genome_basename <- gsub("\\.fa", "", basename(genome))

mm <- mmgenome2::mmload(assembly = genome,
                  coverage = cov_df,
                  verbose = TRUE,
                  kmer_pca = FALSE,
                  kmer_BH_tSNE = FALSE)

#Samples to target
s <- c("D0817M02", "D0817M12")

svars <- cols_by_sample(sample = s, df = mm)


mmgenome2::mmplot(mm, x = svars[1], y = svars[2])

mmgenome2::mmplot_pairs(mm,
             variables = svars,
             #x = "cov_D0606M02metaG_FD",
             #y = "cov_D0606M12metaG_FD",
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.1,
             size_scale = 0.4,
             textsize = 3)
