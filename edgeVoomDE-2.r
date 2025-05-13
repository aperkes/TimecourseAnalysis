
library(edgeR)

### From https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html


### Need to find all counts
setwd("~/Documents/Scripts/TimecourseAnalysis/")
#counts0 <- read.delim("/data/sequencing/TimeFights/results/all_quant.sf",row.names = 1)
counts0 <- read.delim("~/Documents/Scripts/TimecourseAnalysis/all_quant.sf",row.names = 1)

counts0


samples <- read.csv('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName

samples <- samples[order(samples$SampleName),]
samples$Batch <- as.factor(samples$Batch)
samples.23 <- samples[samples$Batch != 1,]
samples.23A <- samples.23[samples.23$Time == 'A',]
#samples <- samples.23A
head(samples[,c("Treatment",'Time','Batch')])

## Create DGEList object
d0 <- DGEList(counts0)

## Calculate normalization factors
d0 <- calcNormFactors(d0)

d0
### Figure out cutoff (might want to fiddle with this)


## Get sample names
snames <- colnames(counts0) # Sample names
snames

time0 <- substr(snames, 3,3)
treatment0 <- substr(snames, 2,2)
time0
treatment0[1] <- "L" ## There's one typo in the sample names
treatment0
treatment0[treatment0 == "L"] <- "Z" ## Z for loZer
day0 <- time0
day0[time0 == "A"] <- 1
day0[time0 == "B"] <- 1
day0[time0 == "G"] <- 1
day0[time0 == "D"] <- 2
day0[time0 == "E"] <- 3
day0[time0 == "F"] <- 4

batch0 <- substr(snames, 4,4)
batch0

d23 <- d0[,batch0 != 1]
batch23 <- batch0[batch0 != 1]
treatment23 <- treatment0[batch0 != 1]
time23 <- time0[batch0 != 1]
counts23 <- counts0[,batch0 != 1]

d.expected <- d0[,samples$Upset != 1]
batch.expected <- batch0[samples$Upset != 1]
treatment.expected <- treatment0[samples$Upset != 1]
time.expected <- time0[samples$Upset != 1]
counts.expected <- counts0[,samples$Upset != 1]

d.contest <- d0[,samples$Upset == 0 & samples$Contest == 1]
batch.contest <- batch0[samples$Upset == 0 & samples$Contest == 1]
treatment.contest <- treatment0[samples$Upset == 0 & samples$Contest == 1]
time.contest <- time0[samples$Upset == 0 & samples$Contest == 1]
counts.contest <- counts0[,samples$Upset == 0 & samples$Contest == 1]

FILTER.DATA = "NoUpset"
FILTER.DATA = "Contests"
#FILTER.DATA = "None"
#FILTER.DATA = "Batch23"

if (FILTER.DATA == "Batch23") { 
  d1 <- d23
  
  batch <- batch23
  time <- time23
  treatment <- treatment23
  counts <- counts23
} else if (FILTER.DATA == "NoUpset") {
  d1 <- d.expected
  batch <- batch.expected
  time <- time.expected
  treatment <- treatment.expected
  counts <- counts.expected
} else if (FILTER.DATA == "Contests") {
  d1 <- d.contest
  batch <- batch.contest
  time <- time.contest
  treatment <- treatment.contest
  counts <- counts.contest
} else { 
  d1 <- d0
  
  batch <- batch0
  time <- time0
  treatment <- treatment0
  counts <- counts0
  }

groups <- interaction(treatment,time)
groups

FILTER = "GROUP"

if (FILTER == "INDV") {
  cutoff <- 1                         ## This sets the cutoff for max count 
  drop <- which(apply(cpm(d0), 1, max) < cutoff) ## in at least one sample
  d <- d1[-drop,] 
  dim(d) # number of genes left
} else if (FILTER == "NONE") {
  d <- d1 ## This does no filtering
  dim(d)
} else {
  ### Alternatively, get every gene that has > 1 for 2 of at least 1 group
  countsPM.bool <- cpm(counts) > 1
  countsPM.num <- 1*countsPM.bool
  group_sums <- rowsum(t(countsPM.num),as.numeric(groups))
  
  ## Check group sums: 
  apply(group_sums,1,FUN=max)
  
  ## What this is saying is that it's in >= n
  ## samples for at least one group. 
  keep.group <- apply(group_sums,2,FUN=max) >= 2  ## This approach leaves 12000
  d <- d1[keep.group,]
  dim(d)
}

color.key <- data.frame(treatment=c("L","S","W"),
                        color=c("blue","grey","gold"))
treatment.df <- data.frame(treatment=treatment23)
color.df <- merge(treatment.df,color.key,by="treatment")

#plotMDS(d,col = as.numeric(batch))
#plotMDS(d,col = time23)
#plotMDS(d,col = color.df$color)

GROUPS <- T
CORR <- T

### Standard group x group comparison
if (GROUPS) { 
  mm <- model.matrix(~0 + groups)
  y <- voom(d,mm,plot=T)
} else {
  mm <- model.matrix(~0 + time*treatment)
  colnames(mm) <- make.names(colnames(mm)) ## feels like bad design, but you need to rename it
  colnames(mm)
  y <- voom(d,mm,plot=T)
}

if (CORR) {
  corfit <- duplicateCorrelation(y, mm, block = batch)
  fit.corr <- lmFit(y, mm, block = batch, correlation =
                    corfit$consensus)
  fit <- fit.corr
} else { 
  fit <- lmFit(y,mm)
}

func.checkContrast <- function(fit,coef = NULL,contrasts = NULL,readout = F) { 
  tmp <- contrasts.fit(fit,coef = coef,contrasts = contrasts)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp,sort.by = "P", n= Inf)
  #head(top.table,20)
  n_genes <- length(which(top.table$adj.P.Val < 0.1))
  print(n_genes) # number of DE genes
  unsorted.table <- top.table[order(row.names(top.table)),]
  return(list(unsorted.table,n_genes))
  }

if (GROUPS == F) {
  res.a <- func.checkContrast(fit2,coef = "timeA")
  res.b <- func.checkContrast(fit2,coef = "timeB")
  res.c <- func.checkContrast(fit2,coef = "timeG")
  res.d <- func.checkContrast(fit2,coef = "timeD")
  res.e <- func.checkContrast(fit2,coef = "timeE")
  res.f <- func.checkContrast(fit2,coef = "timeF")
} else {
  ## Silly one liner to get all combinations
  comparisons <- unlist(as.list(outer(c('WC','LC','WL'),c('A','B','C','D','E','F'), paste,sep=".")))

  print('Winner vs Control')
  res.wsa <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.A - groupsW.A,levels = colnames(coef(fit))))
  res.wsb <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.B - groupsW.B,levels = colnames(coef(fit))))
  res.wsc <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.G - groupsW.G,levels = colnames(coef(fit))))
  res.wsd <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.D - groupsW.D,levels = colnames(coef(fit))))
  res.wse <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.E - groupsW.E,levels = colnames(coef(fit))))
  res.wsf <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.F - groupsW.F,levels = colnames(coef(fit))))
  
  print('Loser vs Control')
  res.zsa <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.A - groupsZ.A,levels = colnames(coef(fit))))
  res.zsb <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.B - groupsZ.B,levels = colnames(coef(fit))))
  res.zsc <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.G - groupsZ.G,levels = colnames(coef(fit))))
  res.zsd <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.D - groupsZ.D,levels = colnames(coef(fit))))
  res.zse <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.E - groupsZ.E,levels = colnames(coef(fit))))
  res.zsf <- func.checkContrast(fit,contrasts = makeContrasts(groupsS.F - groupsZ.F,levels = colnames(coef(fit))))
  
  print('Winner vs Loser')
  res.zwa <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.A - groupsW.A,levels = colnames(coef(fit))))
  res.zwb <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.B - groupsW.B,levels = colnames(coef(fit))))
  res.zwc <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.G - groupsW.G,levels = colnames(coef(fit))))
  res.zwd <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.D - groupsW.D,levels = colnames(coef(fit))))
  res.zwe <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.E - groupsW.E,levels = colnames(coef(fit))))
  res.zwf <- func.checkContrast(fit,contrasts = makeContrasts(groupsZ.F - groupsW.F,levels = colnames(coef(fit))))

  list.wc <- list(res.wsa[[1]],res.wsb[[1]],res.wsc[[1]],res.wsd[[1]],res.wse[[1]],res.wsf[[1]])
  list.wl <- list(res.zsa[[1]],res.zsb[[1]],res.zsc[[1]],res.zsd[[1]],res.zse[[1]],res.zsf[[1]])
  list.lc <- list(res.zwa[[1]],res.zwb[[1]],res.zwc[[1]],res.zwd[[1]],res.zwe[[1]],res.zwf[[1]])
  
  stacked.ws.log <- do.call("cbind",list(res.wsa[[1]]$logFC,res.wsb[[1]]$logFC,res.wsb[[1]]$logFC,res.wsc[[1]]$logFC,res.wse[[1]]$logFC,res.wsf[[1]]$logFC))
  stacked.wl.log <- do.call("cbind",list(res.zwa[[1]]$logFC,res.zwb[[1]]$logFC,res.zwb[[1]]$logFC,res.zwc[[1]]$logFC,res.zwe[[1]]$logFC,res.zwf[[1]]$logFC))
  stacked.ls.log <- do.call("cbind",list(res.zsa[[1]]$logFC,res.zsb[[1]]$logFC,res.zsb[[1]]$logFC,res.zsc[[1]]$logFC,res.zse[[1]]$logFC,res.zsf[[1]]$logFC))

  stacked.ws <- do.call("cbind",list(res.wsa[[1]]$adj.P.Val,res.wsb[[1]]$adj.P.Val,
                                     res.wsb[[1]]$adj.P.Val,res.wsc[[1]]$adj.P.Val,res.wse[[1]]$adj.P.Val,res.wsf[[1]]$adj.P.Val))
  stacked.wl <- do.call("cbind",list(res.zwa[[1]]$adj.P.Val,res.zwb[[1]]$adj.P.Val,
                                     res.zwb[[1]]$adj.P.Val,res.zwc[[1]]$adj.P.Val,res.zwe[[1]]$adj.P.Val,res.zwf[[1]]$adj.P.Val))
  
  maxes.ws <- apply(stacked.ws, 1, min, na.rm = TRUE)
  maxes.wl <- apply(stacked.wl, 1, min, na.rm = TRUE)
  
  maxes.ws.log <- apply(abs(stacked.ws.log), 1, max, na.rm = TRUE)
  maxes.wl.log <- apply(abs(stacked.wl.log), 1, max, na.rm = TRUE)
  maxes.ls.log <- apply(abs(stacked.ls.log), 1, max, na.rm = TRUE)
  
  sorted.counts <- data.frame(d$counts[order(row.names(d$counts)),])
  
  counts.bygroup <- data.frame(matrix(ncol = 18,nrow = 20708))
  counts.bygroup <- setNames(counts.bygroup,comparisons)
  row.names(counts.bygroup) <- row.names(sorted.counts)
  
  counts.bygroup[,] <- 0
  sorted.counts$WC <- 0 
  sorted.counts$LC <- 0 
  sorted.counts$WL <- 0 
  
  ### Look, I'm bad at R, ok? 
  counts.bygroup[res.zsa[[1]]$adj.P.Val < 0.1,"LC.A"] <- 1
  counts.bygroup[res.zsb[[1]]$adj.P.Val < 0.1,"LC.B"] <- 1
  counts.bygroup[res.zsc[[1]]$adj.P.Val < 0.1,"LC.C"] <- 1
  counts.bygroup[res.zsd[[1]]$adj.P.Val < 0.1,"LC.D"] <- 1
  counts.bygroup[res.zse[[1]]$adj.P.Val < 0.1,"LC.E"] <- 1
  counts.bygroup[res.zsf[[1]]$adj.P.Val < 0.1,"LC.F"] <- 1
  
  counts.bygroup[res.wsa[[1]]$adj.P.Val < 0.1,"WC.A"] <- 1
  counts.bygroup[res.wsb[[1]]$adj.P.Val < 0.1,"WC.B"] <- 1
  counts.bygroup[res.wsc[[1]]$adj.P.Val < 0.1,"WC.C"] <- 1
  counts.bygroup[res.wsd[[1]]$adj.P.Val < 0.1,"WC.D"] <- 1
  counts.bygroup[res.wse[[1]]$adj.P.Val < 0.1,"WC.E"] <- 1
  counts.bygroup[res.wsf[[1]]$adj.P.Val < 0.1,"WC.F"] <- 1
  
  counts.bygroup[res.zwa[[1]]$adj.P.Val < 0.1,"WL.A"] <- 1
  counts.bygroup[res.zwb[[1]]$adj.P.Val < 0.1,"WL.B"] <- 1
  counts.bygroup[res.zwc[[1]]$adj.P.Val < 0.1,"WL.C"] <- 1
  counts.bygroup[res.zwd[[1]]$adj.P.Val < 0.1,"WL.D"] <- 1
  counts.bygroup[res.zwe[[1]]$adj.P.Val < 0.1,"WL.E"] <- 1
  counts.bygroup[res.zwf[[1]]$adj.P.Val < 0.1,"WL.F"] <- 1
  
  ## Do the same but for overall treatment
  sorted.counts[res.zsa[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  sorted.counts[res.zsb[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  sorted.counts[res.zsc[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  sorted.counts[res.zsd[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  sorted.counts[res.zse[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  sorted.counts[res.zsf[[1]]$adj.P.Val < 0.1,"LC"] <- 1
  
  sorted.counts[res.wsa[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  sorted.counts[res.wsb[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  sorted.counts[res.wsc[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  sorted.counts[res.wsd[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  sorted.counts[res.wse[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  sorted.counts[res.wsf[[1]]$adj.P.Val < 0.1,"WC"] <- 1
  
  sorted.counts[res.zwa[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  sorted.counts[res.zwb[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  sorted.counts[res.zwc[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  sorted.counts[res.zwd[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  sorted.counts[res.zwe[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  sorted.counts[res.zwf[[1]]$adj.P.Val < 0.1,"WL"] <- 1
  }


counts.overlap <- data.frame(matrix(nrow = length(comparisons),ncol = length(comparisons)))
counts.overlap <- setNames(counts.overlap,comparisons)
row.names(counts.overlap) <- colnames(counts.overlap)

for (i in comparisons) {
  for (j in comparisons) {
    counts.overlap[i,j] <- sum(counts.bygroup[counts.bygroup[,i] == 1,j])
  }
}

trajectory.ws <- c(res.wsa[[2]],res.wsb[[2]],res.wsc[[2]],res.wsd[[2]],res.wse[[2]],res.wsf[[2]])
trajectory.wl <- c(res.zwa[[2]],res.zwb[[2]],res.zwc[[2]],res.zwd[[2]],res.zwe[[2]],res.zwf[[2]])
trajectory.ls <- c(res.zsa[[2]],res.zsb[[2]],res.zsc[[2]],res.zsd[[2]],res.zse[[2]],res.zsf[[2]])

psorted.counts <- sorted.counts[order(maxes.wl,decreasing = F),]
lsorted.counts <- sorted.counts[order(maxes.wl.log,decreasing = T),]

wl_list <- rownames(psorted.counts[psorted.counts$WL == 1,])
dim(sorted.counts[sorted.counts$WC == 1,])
dim(sorted.counts[sorted.counts$LC == 1,])
dim(sorted.counts[sorted.counts$WL == 1,])

dim(sorted.counts[sorted.counts$WC == 1 & sorted.counts$LC == 1,])
dim(sorted.counts[sorted.counts$WC == 1 & sorted.counts$WL == 1,])
dim(sorted.counts[sorted.counts$WL == 1 & sorted.counts$LC == 1,])

dim(sorted.counts[sorted.counts$WL == 1 & sorted.counts$LC == 1 & sorted.counts$WC == 1,])

dim(sorted.counts[res.wsa[[1]]$adj.P.Val < 0.1 & res.zsa[[1]]$adj.P.Val < 0.1,])

write.table(wl_list,"test.txt",row.names = F,col.names = F)

## Check for overlap: 
image(as.matrix(counts.overlap))

## Check for overlap using top log fold change
## Get order that will split negative from positive
#log_order <- order(abs(maxes.ws.log),decreasing = TRUE)
log_order <- order(abs(list.wl[[2]]$logFC),decreasing = TRUE)
resAll <- list()
for (r in list.wl) {
  resOrdered <- r[log_order,'logFC']
  i_pos <- resOrdered > 0
  i_neg <- resOrdered < 0
  resUp <- resOrdered[resOrdered > 0]
  resDown <- resOrdered[resOrdered < 0]
  resUp <- resUp[1:100]
  resDown <- resDown[1:100]
  resTop <- c(resUp,resDown)
  resAll <- c(resAll,list(resTop))
}

diff_matrix <- matrix(c(resAll[[1]],
                        resAll[[2]],
                        resAll[[3]],
                        resAll[[4]],
                        resAll[[5]],
                        resAll[[6]]),ncol = 6)

## Import the stuff, note, this might break the above code. 
library('RColorBrewer')
library('gplots')

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

#heatmap.2(m_extrema, col=coul)
heatmap.2(diff_matrix,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)


## Checking out upsets: 

## Need to compare upsets to non-upsets

## Grab unexpected wins, compare them to expected wins
rows_of_interest <- (samples$Group == "LE" & samples$Upset == 1) | 
                    (samples$Group == "LF" & samples$Upset == 1) | 
                    (samples$Group == "WE" & samples$Contest == 1 & samples$Upset == 0) | 
                    (samples$Group == "WF" & samples$Contest == 1 & samples$Upset == 0) 
d.upset <- d0[,rows_of_interest]
batch.upset  <- batch0[rows_of_interest]
treatment.upset  <- treatment0[rows_of_interest]
time.upset  <- time0[rows_of_interest]
counts.upset  <- counts0[,rows_of_interest]
snames.upset <- snames[rows_of_interest]
days.upset <- day0[rows_of_interest]
upset <- as.factor(samples[rows_of_interest,"Upset"])


cutoff <- 1                         ## This sets the cutoff for max count 
drop <- which(apply(cpm(d.upset), 1, max) < 1) ## in at least one sample
d1.upset <- d.upset[-drop,] 
dim(d1.upset) # number of genes left

mm.upset <- model.matrix(~upset)
y.upset <- voom(d1.upset,mm.upset,plot=T)

corfit.upset <- duplicateCorrelation(y.upset, mm.upset, block = days.upset)
fit.upset <- lmFit(y.upset, mm.upset, block = days.upset, correlation =
                     corfit.upset$consensus)
colnames(coef(fit.upset))
res.upset <- func.checkContrast(fit.upset,coef="upset1")

head(res.upset[[1]][order(res.upset[[1]]$adj.P.Val),],10)
sig.upset <- rownames(head(res.upset[[1]][order(res.upset[[1]]$adj.P.Val),],10))
write.table(sig.upset,"test.upset.txt",row.names = F,col.names = F)

## so there are some genes different, but we're not controlling for size here, it's sort of a mess.

library(clusterProfiler)


