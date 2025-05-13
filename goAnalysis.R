## Go Analysis: 

## Based on this: https://support.bioconductor.org/p/9144574/

library(AnnotationHub)

hub <- AnnotationHub()

query(hub,c("poecilia","orgdb"))

pot <- hub[["AH114800"]]
potlist <- mapIds(pot, keys(pot), "GOALL", "ENTREZID", multiVals = "list")

## Grab each gene, check go terms, look for go term enrichment. 
sig_list <- read.csv('~/Documents/Scripts/TimecourseAnalysis/all_sig.csv',header=F)
head(sig_list)

sig_list <- as.list(sig_list$V1)

foo <- gost(sig_list,organism = 'pformosa')
