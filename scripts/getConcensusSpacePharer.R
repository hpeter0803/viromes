suppressMessages(library(tidyverse))

setwd("../Desktop/Virus/spacePharer/")


filename <- "allPredictionsTax.tsv"
tax <- read_tsv(filename,col_names = T, show_col_types = F)

duplicatedValues <- tax$phage_acc[duplicated(tax$phage_acc)]

tax_toDo <- tax[tax$phage_acc %in% duplicatedValues,]
tax_done <- tax[!(tax$phage_acc %in% duplicatedValues),]

tax_class <- tax_toDo %>%
  separate(prok_tax,sep = ";", remove = F, into = c("Domain","Phylum", "Class", "Order", "Family", "Genus","Species"))



finalPerc <- 0

VirTax <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(VirTax) <- c("Phage","HostTaxonomy", "NumbMatch", "Perc", "PercNoNA")
#j = "OTE_1_k141_101569"
for (j in unique(tax_class$phage_acc)){
  temp <- tax_class %>%
    filter(phage_acc == j)
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
    taxa <- temp$prok_tax[grepl(tempMax, temp$prok_tax)]
    # taxa[1]
    tax_host = gsub(paste(tempMax,".*$", sep = ""), tempMax, taxa[1], perl = T)
    
    tempdf <- c(j, tax_host, line, perc, percNoNa)
    tempdf <- as.data.frame(t(tempdf))
    
    colnames(tempdf)  <- c("Phage","HostTaxonomy", "NumbMatch", "Perc", "PercNoNA")
    #VirTax <- rbind(VirTax,tempdf)
  }
  VirTax <- rbind(VirTax,tempdf)
}

write_tsv(tax_done, "allPredictionsUniqueMatch.tsv")
write_tsv(VirTax, "allPredictionsMultipleMatch.tsv")

