
library(edgeR)

### From https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html


### Need to find all counts
counts <- read.delim("/data/sequencing/TimeFights/results/all_quant.sf",row.names = 1)

counts

samples <- read.csv('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName

samples$Batch <- as.factor(samples$Batch)
samples.23 <- samples[samples$Batch != 1,]
samples.23A <- samples.23[samples.23$Time == 'A',]
samples <- samples.23A
head(samples[,c("Treatment",'Time','Batch')])

## Create DGEList object
d0 <- DGEList(counts)

## Calculate normalization factors
d0 <- calcNormFactors(d0)

d0
### Figure out cutoff (might want to fiddle with this)


## Get sample names
snames <- colnames(counts) # Sample names
snames

time <- substr(snames, 3,3)
treatment <- substr(snames, 2,2)
time
treatment[1] <- "L"
treatment

batch <- substr(snames, 4,4)
batch
groups <- interaction(treatment,time)
groups

d23 <- d0[batch != 1,]
d <- d23

if (F) {
  cutoff <- 3 ## This sets the cutoff for max count in at least one sample
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  dim(d) # number of genes left
} else if (T) {
  d <- d0 ## This does no filtering
  dim(d)
} else {
  ### Alternatively, get every gene that has > 1 for 4 of at least 1 group
  countsPM.bool <- cpm(counts) > 1
  countsPM.num <- 1*countsPM.bool
  group_sums <- rowsum(t(countsPM.num),as.numeric(groups))
  
  ## Check group sums: 
  apply(group_sums,1,FUN=max)
  
  ## What this is saying is that it's in at least 
  ## four samples for at lesat one group. 
  keep.group <- apply(group_sums,2,FUN=max) >= 2  ## This approach leaves 12000
  d <- d0[keep.group,]
  dim(d)
}

plotMDS(d,col = as.numeric(batch))

mm <- model.matrix(~0 + groups)
y <- voom(d,mm,plot=T)

fit <- lmFit(y,mm)
head(coef(fit))
contr <- makeContrasts(groupsW.A - groupsL.A,levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp,sort.by = 'P',n=Inf)
head(top.table,20)

length(which(top.table$adj.P.Val < 0.05))

mm2 <- model.matrix(~treatments*times)
colnames(mm2)
y <- voom(d,mm2,plot=T)
fit <- lmFit(y,mm2)

head(coef(fit))
tmp <- contrasts.fit(fit, coef = 3) # Directly test second coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05)) # number of DE genes
