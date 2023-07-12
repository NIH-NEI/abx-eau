# Author : Vijay Nagarajan PhD
# Institute : NEI/NIH
# This scripts reads in a tab delimited gene set file and generates the gmt format input file for gsea
# load the data science library
library(tidyverse)

#read/import the data
amigo = read.delim("/Users/Mouse_PanCancer_Immune_Profiling_ListABXEAU.txt", skip = 1, header = TRUE, check.names = FALSE)
setwd("/Users/myanalysis/")

#look at the imported data
View(amigo)
View(as.data.frame(colnames(amigo)))
#V5 is go ids, V6 is go description, V2 is gene names
#exract the unique GO terms (column V5)
#golist = unique(amigo %>% pull(colnames))
golist = colnames(amigo)
golist=golist[-c(1,2,28)]
golist

#for every unique goterm, extract relevant genes (column V2)
for (goone in golist)
{
  print(goone)
  gooneline=""
  goonegenes=amigo$`Gene Name`[amigo[[goone]]=="+"]
  #gooneline=c(goone, unique(gooneall$V5), gooneall$V3)
  
  gooneline=c(goone, goone, goonegenes)
  print(gooneline)
  gooneline=paste(gooneline, collapse="\t")
  lapply(gooneline, write, 'tcell_reiko_annotated.gmt', append=TRUE)
}
