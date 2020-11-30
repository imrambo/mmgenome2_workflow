#!/usr/bin/env Rscript

#Author: Ian Rambo
#Thirteen... that's a mighty unlucky number... for somebody!

#Script to clean 178.concoct.D0608WATERmetaG_sub

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
genome_pre_name <- "178.concoct.D0608WATERmetaG_sub.fa"
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
sreg <- "WATER"

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
             variables = svars,
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.1,
             size_scale = 0.4,
             textsize = 3)

### Create a scatter plot with 2d density overlay
scatter_density2d_gg <- ggplot(mm, aes(x = log10(cov_D0608WATERmetaG_FD),
               y = log10(cov_D0606WATERmetaG_FD))) +
  geom_point(aes(size = `length`), alpha = 0.6) +
  geom_density_2d()

mmgenome2::mmplot(mm,
                  x = "cov_D0608WATERmetaG_FD",
                  y = "cov_D0606WATERmetaG_FD",
                  x_scale = "log10",
                  y_scale = "log10",
                  locator = TRUE)

#List of shiny selection data frames
mag_selections <- list(data.frame(cov_D0608WATERmetaG_FD = c(456.259, 339.499, 240.854, 190.154, 126.805, 96.788, 67.137, 55.758, 48.441, 48.989, 65.643, 119.193, 209.244, 310.266, 489.441),
                                  cov_D0606WATERmetaG_FD = c(417.228, 284.029, 191.456, 138.277, 87.854, 55.818, 47.204, 47.204, 69.34, 94.132, 139.648, 257.359, 413.135, 618.968, 643.001)),
                       data.frame(cov_D0608WATERmetaG_FD = c(4.21, 3.479, 3.174, 2.942, 4.163, 5.901, 8.315, 11.717, 16.132, 22.21, 25.094, 30.401, 33.364, 37.478, 36.83, 30.401, 24.234, 15.579, 10.676, 5.192),
                                  cov_D0606WATERmetaG_FD = c(2.135, 2.507, 3.872, 5.407, 6.962, 10.434, 16.118, 22.733, 29.871, 41.286, 51.056, 49.032, 41.706, 32.718, 25.152, 16.282, 10.225, 5.63, 3.571, 2.069))
                       )


write_genomes(slist = mag_selections, extension = "fna", mmdf = mm, gpn = genome_post_name, outdir = clean_dir)
