#!/usr/bin/env Rscript
# Rscript to get consensus output from PHIST

suppressMessages(library(tidyverse))

# Logging file
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

# filename <- "allPredictions.tsv"
filename <- snakemake@input[[""]]
tax <- read_tsv(filename,col_names = T, show_col_types = F)
tax <- tax %>%
  filter(`adj-pvalue` <= 0.05)

duplicatedValues <- tax$phage[duplicated(tax$phage)]

tax_toDo <- tax[tax$phage %in% duplicatedValues,]
tax_done <- tax[!(tax$phage %in% duplicatedValues),]

tax_class <- tax_toDo %>%
  separate(Taxa,sep = ";", remove = F, into = c("Domain","Phylum", "Class", "Order", "Family", "Genus","Species"))



finalPerc <- 0

VirTax <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(VirTax) <- c("Phage","HostTaxonomy", "NumbMatch", "Perc", "PercNoNA")
# j = "VAR_69_k141_824765"
for (j in unique(tax_class$phage)){
  temp <- tax_class %>%
    filter(phage == j)
  tempdf = ""
  line <- nrow(temp)
  for (i in c("Domain","Phylum", "Class", "Order", "Family", "Genus")){
    col <- temp[[i]]
    col <- gsub(".__$","",col, perl = T)
    col[col==""] <- NA
    col_noNA <- na.omit(col)
    linesNoNA <- length(col_noNA)
    
    col_noNA <- as.data.frame(col_noNA)
    if(nrow(col_noNA) == 0) next
    col_uni <- col_noNA %>% count(col_noNA)
    
    val_max <- max(col_uni$n)
    row_max <- which.max(col_uni$n)
    tempMax <- col_uni[row_max,1]
    
    perc <- val_max/length(col)
    percNoNa <- val_max/linesNoNA
    
    if(percNoNa >= 0.75 || (perc > 0.1 && percNoNa >= finalPerc)){
      finalPerc <- percNoNa
    }
    else{
      break
    }
    taxa <- temp$Taxa[grepl(tempMax, temp$Taxa)]
    # taxa[1]
    tax_host = gsub(paste(tempMax,".*$", sep = ""), tempMax, taxa[1], perl = T)
    
    tempdf <- c(j, tax_host, line, perc, percNoNa)
    tempdf <- as.data.frame(t(tempdf))
    
    colnames(tempdf)  <- c("Phage","HostTaxonomy", "NumbMatch", "Perc", "PercNoNA")
    #VirTax <- rbind(VirTax,tempdf)
  }
  VirTax <- rbind(VirTax,tempdf)
}

write_tsv(tax_done, snakemake@output[["UNIQ"]])
write_tsv(VirTax, snakemake@output[["MULT"]])

