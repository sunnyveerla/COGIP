#!/usr/bin/env Rscript

# Title: CN_sig.R
# Author: Srinivas Veerla
# Date:	14/10/24

# Description: this script is used for signature validation on each sample group.

library(NMF)
library(flexmix)
library(YAPSA)
library(tidyverse)
library(lsa)

main_script = snakemake@input[['main']]
helper_script = snakemake@input[['helper']]

# load the functions
source(main_script)
source(helper_script)

# specify the output directory
outdir <- snakemake@params[['outdir']]

# extract the required values from the snakemake workflow
sample_file = snakemake@params[['sample_info']]
print(sample_file)
sample_info <- read.table(sample_file,header=T,sep="\t") # the sample information tsv file
sample_names <- sample_info$sample_name
print(sample_names)
in_dir <- snakemake@params[['indir']] # the directory to store the segment files
binsize <- snakemake@params[['binsize']]

# generate the required input dataframes list
segment_list <- list()
for (i in sample_names) { 
  tab_df <- read.table(paste0(in_dir, i, '/05_absolute_CN/', i, '_', binsize, 'kb_seg.tsv'), sep = '\t', header = 1) %>% select('chromosome','start','end','segVal')
  a <- length(segment_list) + 1
  segment_list[[a]] <- tab_df  
}
names(segment_list) <- sample_names

# validate the signatures in the samples
CN_features <- extractCopynumberFeatures(segment_list)
sample_by_component <- generateSampleByComponentMatrix(CN_features)
print(sample_by_component)

# output the matrix files, matrix object(RDS) and a simple heatmap
output_prefix <- paste0(outdir,'CN_sig_',binsize,'kb')
# the same output for the sample-by-component matrix should also be saved for future analysis
write.table(sample_by_component, file = paste0(output_prefix, '.SCmatrix.txt'))
saveRDS(sample_by_component, file = paste0(output_prefix, '.SCmatrix.rds'), compress = FALSE)
pdf(paste0(output_prefix, 'heatmapSC.pdf'))
heatmap(sample_by_component)
dev.off()