#!/usr/bin/env Rscript

#Author: Ian Rambo
#Thirteen... that's a mighty unlucky number... for somebody!

#Script to clean 148.concoct.D1103S02metaG

library(mmgenome2)
library(tidyverse)
library(shiny)

rm(list=ls())
#=============================================================================
cols_by_sample_name <- function(sample_name, df) {
  ### Return a list of column names for a set of sample_names
  sample_name_patt = paste(sample_name, collapse = "|")
  sample_name_vars = grep(sample_name_patt, colnames(df), value = TRUE)
  return(sample_name_vars)
}

write_genomes <- function(slist, extension, mmdf, gpn, outdir){
  #Write output genome subsets from a list of shiny selection data frames
  for(i in seq(1, length(slist))){
    magid = paste("mag", as.character(i), sep = "-")
    gpn_full = paste(gpn, paste(magid, extension, sep = "."), sep = "_")
    mmsubset = mmextract(mmdf, selection = mag_selections[[i]])
    outfile = file.path(outdir, gpn_full)
    print(paste("exporting", outfile, sep = " "))
    mmexport(mmsubset, assembly = assembly, file = outfile)
  }
  
}
#=============================================================================
###---Input/Output
#-----------------------------------------------------------------------------
###---EDIT THESE PATHS

#Directory containing depth files
cov_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/mmgenome2/cov_files"
#Directory containing genome files
genome_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/das_tool/DASTool_Run2_concoct-maxbin2-metabat2-vamb_diamond_DASTool_bins"
#Create output directory for cleaned genomes
clean_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/DASTool_Run2_bins_cleaned"
dir.create(clean_dir)

#Pre-clean name of genome file
genome_pre_name <- "148.concoct.D1103S02metaG.fa"
#-----------------------------------------------------------------------------
#Path to genome
genome_pre_path <- file.path(genome_dir, genome_pre_name)
genome_basename <- gsub("\\.fa", "", basename(genome_pre_path))

#Post-clean name of genome file
genome_post_name <- paste(genome_basename, "cleaned", sep = "_")

#Path to genome coverage file
cov_file <- list.files(file.path(cov_dir, genome_basename),
                       pattern = ".*_cov", full.names = TRUE)

#Depth file data frame
cov_df <- read.table(cov_file, header = TRUE, sep = "\t")

colnames(cov_df) <- gsub("assembly_metaspades_cat_renamed_3kb_|\\.coordsort\\.bam",
                         "", colnames(cov_df))

#Coverage data in long format
cov_df_long <- cov_df %>%
  tidyr::gather(colnames(cov_df)[2:length(colnames(cov_df))],
                key = "sample_name", value = "coverage")

#=============================================================================
###---Identify outliers

#Quartile values and Interquartile Ranges
cov_quartiles <- cov_df_long %>% group_by(sample_name) %>%
  summarise(tibble::enframe(quantile(coverage, c(0.25, 0.5, 0.75)),
                            "quantile", "q_coverage")) %>%
  distinct() %>%
  tidyr::spread(quantile, q_coverage) %>%
  rename(quant_25 = `25%`, quant_50 = `50%`, quant_75 = `75%`) %>%
  mutate(iqr = quant_75 - quant_25,
         upper_outlier = quant_75 + (iqr * 1.5),
         lower_outlier = quant_25 - (iqr * 1.5))

# Percent of contig depth values outside of outlier thresholds for
# each sample
cov_pct_outlier <- cov_quartiles %>% left_join(cov_df_long) %>%
  filter(coverage > upper_outlier | coverage < lower_outlier) %>%
  group_by(sample_name) %>%
  count() %>%
  mutate(pct_outlier = n/nrow(cov_df)) %>%
  dplyr::arrange(desc(pct_outlier))

#Boxplot of log10 coverage values by sample
cov_boxplot <- cov_df_long %>% ggplot(aes(x = sample_name,
                                          y = coverage)) +
  geom_boxplot() +
  theme(axis.line = element_line(color="black"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(color = "black", angle=65, hjust = 1, size = 7.5),
        axis.text.y=element_text(color="black", size = 10),
        strip.text.x = element_text(size = 11, colour = "black"),
        strip.text.y = element_text(size = 11, colour = "black")) +
  ggtitle(genome_basename)

#Boxplot of log10 coverage values by sample
cov_boxplot_log10 <- cov_df_long %>% ggplot(aes(x = sample_name,
                                                y = log10(coverage))) +
  geom_boxplot() +
  theme(axis.line = element_line(color="black"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(color = "black", angle=65, hjust = 1, size = 7.5),
        axis.text.y=element_text(color="black", size = 10),
        strip.text.x = element_text(size = 11, colour = "black"),
        strip.text.y = element_text(size = 11, colour = "black")) +
  ggtitle(genome_basename)

#=============================================================================
###---Coverage Graphics 

mm <- mmgenome2::mmload(assembly = genome_pre_path,
                  coverage = cov_df,
                  verbose = TRUE,
                  kmer_pca = FALSE,
                  kmer_BH_tSNE = FALSE) 

###---EDIT PAST THIS POINT---###
#Regex of sample_names to target
sreg <- "S[0-9]{2}"

cov_pct_outlier_smp <- cov_pct_outlier %>% filter(grepl(sreg, sample_name)) %>%
  mutate(sample_name = gsub("^", "cov_", sample_name))

#Sample names for highest and lowest percent outliers
#Number of high and low samples
pr <- 3
pvars <- c(cov_pct_outlier_smp$sample_name[1:pr],
           cov_pct_outlier_smp$sample_name[(nrow(cov_pct_outlier_smp) - pr):nrow(cov_pct_outlier_smp)])

#Sample names based on regex
svars <- cols_by_sample_name(sample_name = sreg, df = mm)

#Plot multiple sample coverages 
mmpairs <- mmgenome2::mmplot_pairs(mm,
             variables = pvars,
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.1,
             size_scale = 0.4,
             textsize = 3)

### Create a scatter plot with 2d density overlay
scatter_density2d_gg <- ggplot(mm, aes(x = log10(cov_D0819S20remetaG_FD),
               y = log10(cov_D0819S12metaG_FD))) +
  geom_point(aes(size = `length`), alpha = 0.6) +
  geom_density_2d()

scatter_density2d_gg <- ggplot(mm, aes(x = log10(cov_D0606S02metaG_FD),
                                       y = log10(cov_D0608S02metaG_FD))) +
  geom_point(aes(size = `length`), alpha = 0.6) +
  geom_density_2d()

mmgenome2::mmplot(mm,
                  x = "cov_D0606S02metaG_FD",
                  y = "cov_D0608S02metaG_FD",
                  x_scale = "log10",
                  y_scale = "log10",
                  locator = TRUE)

#List of shiny selection data frames
mag_selections <- list(data.frame(cov_D0606S02metaG_FD = c(6.563, 7.045, 7.213, 7.069, 6.282, 5.434, 3.985, 3.158, 2.731, 2.519, 2.316, 2.172, 1.829, 1.582, 1.346, 1.129, 1.003, 0.919, 0.868, 0.871, 0.895, 0.922, 0.922, 1.137, 2.485, 2.922, 4.039, 4.622, 6.053),
                                  cov_D0608S02metaG_FD = c(14.611, 12.342, 9.773, 7.008, 5.025, 4.141, 3.167, 2.662, 2.172, 1.947, 1.686, 1.604, 1.596, 1.424, 1.366, 1.57, 1.833, 2.254, 2.701, 3.17, 3.664, 4.278, 5.022, 5.833, 9.384, 12.728, 14.335, 16.566, 16.912))
                   )


write_genomes(slist = mag_selections, extension = "fna", mmdf = mm, gpn = genome_post_name, outdir = clean_dir)
