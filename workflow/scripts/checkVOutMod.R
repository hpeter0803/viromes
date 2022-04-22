#!/usr/bin/env Rscript
# Rscript to filter the checkV output and keep only the "good quality" contigs

suppressMessages(library(tidyverse))

# Logging file
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

# args <- commandArgs(trailingOnly = TRUE)
filename <- snakemake@input[["qual_in"]]

dat <- read_tsv(filename, show_col_types = FALSE)

dat_filt <- dat %>%
        filter(!checkv_quality %in% c("Low-quality","Not-determined"))

fileOut <- gsub("quality_summary.tsv", "goodQual.tsv", filename)

write_tsv(dat_filt, snakemake@output[["qual_out"]])

# collecting only "complete" contigs
commplete <- dat %>%
  filter(checkv_quality %in% c("Complete"))

fileOut <- gsub("quality_summary.tsv", "complete_contigs.tsv", filename)
write_tsv(dat_filt, snakemake@output[["complete"]])
