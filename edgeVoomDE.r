
library(edgeR)

### From https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html


### Need to find all counts
counts <- read.delim("/data/sequencing/TimeFights/results/all_quant.sf",row.names = 1)

counts

## Create DGEList object
d0 <- DGEList(counts)

## Calculate normalization factors
d0 <- calcNormFactors(d0)

### Figure out cutoff (might want to fiddle with this)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

## Get sample names
snames <- colnames(counts) # Sample names
snames

