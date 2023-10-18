suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

dat <- read_tsv(filename,show_col_types = FALSE)

dat_filt <- dat %>%
        filter(!checkv_quality %in% c("Low-quality","Not-determined"))

fileOut <- gsub("/quality_summary.tsv", "_goodQual.tsv", filename)

write_tsv(dat_filt, fileOut)