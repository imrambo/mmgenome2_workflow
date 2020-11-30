#!/usr/bin/env Rscript

#Author: Ian Rambo
#Thirteen... that's a mighty unlucky number... for somebody!

#Script to REFINE 178.concoct.D0608WATERmetaG_sub_cleaned_mag-1

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
#cov_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/mmgenome2/cov_files"
#Directory containing genome files
genome_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/DASTool_Run2_bins_cleaned"
#Create output directory for cleaned genomes
clean_dir <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/DASTool_Run2_bins_cleaned"
dir.create(clean_dir)


#Name of cleaned genome to refine
genome_pre_name <- "178.concoct.D0608WATERmetaG_sub_cleaned_mag-1.fna"
#-----------------------------------------------------------------------------
#Path to genome
genome_pre_path <- file.path(genome_dir, genome_pre_name)
genome_basename <- gsub("\\.fna", "", basename(genome_pre_path))

#Post-clean name of genome file
genome_post_name <- paste(genome_basename, "refined", sep = "_")

#Path to genome coverage file
cov_file <- "/Users/ian/Documents/phd_research/MANERR_JGI/analysis/metaG/mmgenome2/mmgenome2/cov_files/178.concoct.D0608WATERmetaG_sub/178.concoct.D0608WATERmetaG_sub_cov"

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
mag_selections <- list(data.frame(cov_D0608WATERmetaG_FD = c(385.052, 404.854, 408.073, 387.09, 295.725, 241.337, 182.92, 131.862, 101.808, 86.211, 78.19, 67.446, 62.476, 56.216, 55.92, 55.92, 56.514, 123.116, 137.187, 143.483),
                                  cov_D0606WATERmetaG_FD = c(489.948, 455.119, 418.719, 374.273, 242.003, 188.458, 134.589, 95.199, 77.043, 71.109, 67.553, 64.381, 67.77, 79.045, 86.192, 96.119, 109.974, 256.38, 306.803, 331.342))
                       )


write_genomes(slist = mag_selections, extension = "fna", mmdf = mm, gpn = genome_post_name, outdir = clean_dir)
