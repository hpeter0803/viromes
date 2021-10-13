suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

filename <- "VAR_61_k141_1551594.txt"
tax <- read_tsv(filename,col_names = F, show_col_types = F)
colnames(tax) <- c("Query_ID","Subject_ID","Percentage_of_identical_matches","Alignment_length","Number_of_mismatches","Number_of_gap_openings","Start_of_alignment_in_query","End_of_alignment_in_query","Start_of_alignment_in_subject","End_of_alignment_in_subject","Expected_value","Bit_score")
tax$Subject_ID <- gsub("\\|.*" , "",tax$Subject_ID,perl =TRUE)

host <- read_tsv("IMGVR_all_Host_information.tsv",show_col_types = F)
names(host) <- gsub(x = names(host), pattern = "## ", replacement = "")  

seqInfo <- read_tsv("IMGVR_all_Sequence_information.tsv",show_col_types = F)
names(seqInfo) <- gsub(x = names(seqInfo), pattern = "## ", replacement = "")  

seqInfo_subset <- seqInfo %>%
  select(c(UViG, vOTU,`Taxonomic classification`, `Host taxonomy prediction`, `Sequence origin (doi)`, `Gene content (total genes;cds;tRNA;VPF percentage)`))

tax_merge <- merge(tax, seqInfo_subset, by.x = "Subject_ID", by.y = "UViG")

tax_class <- tax_merge %>%
  select(c(Query_ID, Subject_ID, `Taxonomic classification`, `Host taxonomy prediction`))

tax_class <- tax_class %>%
  separate(`Taxonomic classification`,sep = ";", remove = F, into = c("VRealm","VKingdom", "VPhylum", "VClass", "VOrder", "VFamily", "VGenus", "VSpecies")) %>%
  separate(`Host taxonomy prediction`,sep = ";", remove = F, into = c("BDomain","BPhylum", "BClass", "BOrder", "BFamily", "BGenus", "BSpecies"))


finalPerc <- 0

bacTax <- ""
# i <- "BKingdom"
for (i in c("BDomain","BPhylum", "BClass", "BOrder", "BFamily", "BGenus", "BSpecies")){
  col <- tax_class[[i]]
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
  
  tempdf <- c(i, val_max, perc, percNoNa, tempMax)
  tempdf <- as.data.frame(t(tempdf))
  colnames(tempdf)  <- c("TaxLevel", "Value", "PercB", "PercNoNA" , "Taxa")
  bacTax <- tempdf
}

finalPerc <- 0
VirTax <- ""
#i <- "BKingdom"
for (i in c("VRealm","VKingdom", "VPhylum", "VClass", "VOrder", "VFamily", "VGenus", "VSpecies")){
  col <- tax_class[[i]]
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
  
  if(perc > 0.1 && percNoNa >= finalPerc){
    finalPerc <- percNoNa
  }
  else{
    break
  }
  
  tempdf <- c(i, val_max, perc, percNoNa, tempMax)
  tempdf <- as.data.frame(t(tempdf))
  colnames(tempdf)  <- c("TaxLevel", "Value", "PercV", "PercNoNA" , "Taxa")
  VirTax <- tempdf
}

# a <- grepl(VirTax$Taxa, x = tax_merge$`Taxonomic classification`)
# b <- unique(tax_merge$`Taxonomic classification`[a])
# VirTax$Taxa <- b
# 
# a <- grepl(bacTax$Taxa, x = tax_merge$`Host taxonomy prediction`)
# c <- unique(tax_merge$`Host taxonomy prediction`[a])
# bacTax$Taxa <- c

sample <- gsub(".txt", "", filename)
fileout <- gsub(".txt", "_temp.txt", filename)

fileout <- gsub(".*\\/","", fileout, perl = T)

fileout<- paste("Temp/", fileout, sep = "")

if(!is.data.frame(VirTax)){
  VirTax <-  data.frame(matrix(ncol = 5, nrow = 0))
  VirTax[1,] <- NA
  colnames(VirTax) <- c("TaxLevel", "Value", "PercV", "PercNoNA" , "Taxa")
}
  
  
if(!is.data.frame(bacTax)){
  bacTax <- data.frame(matrix(ncol = 5, nrow = 0))
  bacTax[1,] <- NA
  colnames(bacTax) <- c("TaxLevel", "Value", "PercB", "PercNoNA" , "Taxa")
}



datfinal <- cbind(sample,VirTax$Taxa, VirTax$PercV, VirTax$TaxLevel, bacTax$Taxa, bacTax$PercB,bacTax$TaxLevel)
datfinal <- as.data.frame(datfinal)
write_tsv(datfinal, fileout, col_names = F)
