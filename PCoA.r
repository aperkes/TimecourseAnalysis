
BiocManager::install('DESeq2','edgeR','limma','statmod','affycoretools','ReportingTools')

install.packages("pheatmap")
install.packages("vegan")
install.packages("ape")
install.packages("rgl")
library("edgeR")
library("statmod")
library('DESeq2')
library('vegan')
library('ape')

setwd("/data/sequencing/TimeFights/analysis")

countdata_1 <- read.delim("all_quant.sf", header=TRUE, row.names="Name") #read in the file of the count data and call it countdata, row.names tells the name in the top left cell, the gene names

countdata_1 <- round(countdata_1)

keep <- rowSums(cpm(countdata_1)>1) >=4
# keeps rows (genes) where at least 4 columns (libraries) have at least 1 count per million. This means that if a gene is only expressed in say one treatment (which has three replicates), this gene will not be thrown out of the analysis

countdata <- countdata_1[keep,] #formatting for organizing the kept rows that summed to at least 1 cpm in the step above

countdata <- round(countdata)

## Need to sort both so that they are same order
countdata <- countdata[,order(colnames(countdata))]

coldata <- read.delim("./columnInfo.txt", header=TRUE, row.names=1) #use all_v3 for all data
coldata <- coldata[order(rownames(coldata)),]

## optionally drop the first batch to remove batch effect
countdata <- countdata[,coldata$Batch != 1]
coldata <- coldata[coldata$Batch != 1,]


# Include this last row if you want to use controlA as a t0 comparison
win_rows <- coldata$PlannedTreatment == 'win' | coldata$Time == 'A0'
loss_rows <- coldata$PlannedTreatment == 'loss' | coldata$Time == 'A0'
control_rows <- coldata$PlannedTreatment == 'control' 
B_rows <- coldata$Time == 'B'
F_rows <- coldata$Time == 'F'

col_win <- coldata[win_rows,]
col_loss <- coldata[loss_rows,]
col_c <- coldata[control_rows,]
col_B <- coldata[B_rows,]
col_F <- coldata[F_rows,]

win_counts <- countdata[,win_rows]
loss_counts <- countdata[,loss_rows]
control_counts <- countdata[,control_rows]
B_counts <- countdata[,B_rows]
F_counts <- countdata[,F_rows]

## Here I can drop controls, since they are not quite balanced
#coldata <- coldata[-36:-68,]
#countdata <-countdata[,-36:-68]

#ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ POP*TREATMENT) #use POP*TREATMENT for tp 1 and 2, POP*TREATMENT*TP for all data
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ Batch + PlannedTreatment + Time) # Batch effects can't be random effect here
ddsInteractionCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ Batch + PlannedTreatment*Time) # Batch effects can't be random effect here


ddsWinTable <- DESeqDataSetFromMatrix(countData=win_counts, colData=col_win, design = ~ Batch + Time) 
ddsLossTable <- DESeqDataSetFromMatrix(countData=loss_counts, colData=col_loss, design = ~ Batch + Time) 
ddsControlTable <- DESeqDataSetFromMatrix(countData=control_counts, colData=col_c, design = ~ Batch + Time) 
ddsBTable <- DESeqDataSetFromMatrix(countData=B_counts, colData=col_B, design = ~ Batch + PlannedTreatment) 
ddsFTable <- DESeqDataSetFromMatrix(countData=F_counts, colData=col_F, design = ~ Batch + PlannedTreatment) 

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)

ddsWin <- DESeq(ddsWinTable)
ddsLoss <- DESeq(ddsLossTable)
ddsControl <- DESeq(ddsControlTable)
ddsB <- DESeq(ddsBTable)
ddsF <- DESeq(ddsFTable)

res_win <- results(ddsWin,contrast=c('Time','F','A0'))
res_win
rwinOrdered <- res_win[order(res_win$padj),]
head(rwinOrdered)
sum(res_win$padj<0.05, na.rm=TRUE)

res <- results(ddsFull,contrast=c('PlannedTreatment','win','loss'))
res
resOrdered <- res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05, na.rm=TRUE)


ddsSelect <-ddsFull

ddsSelect <- ddsWin
ddsSelect <- ddsLoss
ddsSelect <- ddsControl
ddsSelect <- ddsB
ddsSelect <- ddsF

res <- results(ddsSelect,contrast=c('PlannedTreatment','win','loss'))
#res <- results(ddsSelect,contrast=c('Time','A0','F'))
resOrdered <- res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05,na.rm=TRUE)



rld <- rlogTransformation(ddsSelect) # Not sure whether there are advantages

rld <- vst(ddsSelect) # Use if there are many samples (>30)

# assembling table of conditions to label PCoA plot:

# There's a big batch effect, so for the PCA to be useful, 
# probably best to drop the first batch. Or average treatments.

# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1 <- as.character(colData(ddsSelect)$PlannedTreatment) 
factor2 <- as.character(colData(ddsSelect)$Experience) 
factor3 <- as.character(colData(ddsSelect)$Batch) 
factor4 <- as.character(colData(ddsSelect)$Time) 



# actual PCoA analysis
vsd <- assay(rld)
dds.pcoa <- pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores <- dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
weights <- percent / sum(percent) * 100 #percent for each axes


#plot PC axis 1 and 2 for all data
xLab <- paste("PC1 (",round(weights[1],3),"%)")
yLab <- paste("PC2 (",round(weights[2],3),"%)")

## For batches:
plot(scores[,1], scores[,2], col=c("blue","red","orange")[as.numeric(as.factor(factor3))], pch=c(17, 4, 16)[as.numeric(as.factor(factor1))], xlab = xLab, ylab = yLab)
ordispider(scores,factor3, col=c("blue","red","orange"))



## For just time = B
plot(scores[,1], scores[,2], col=c("green","red","grey")[as.numeric(as.factor(factor1))], pch=c(0, 15, 1)[as.numeric(as.factor(factor1))], xlab = xLab, ylab = yLab)
ordispider(scores,factor1, col=c("darkgreen","darkred","darkgrey"))


plot(scores[,1], scores[,2], col=c("darkblue","blue","yellow","orange","red","lightblue")[as.numeric(as.factor(factor4))], pch=c(0, 15, 1, 16)[as.numeric(as.factor(factor1))], xlab = xLab, ylab = yLab)
ordispider(scores,factor4, col=c("darkblue","blue","yellow","orange","red","lightblue"))


 $#plot PC axis 3 and 4 for all data
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 15, 1, 16)[as.numeric(as.factor(factor1))], xlab = "PC1 (42.47%)", ylab = "PC2 (17.22%)")


#testing significance between group variances (need to do this because adonis will give sig p value even if groups overlap, but its sig b/c the variance is different)
#https://github.com/vegandevs/vegan/issues/233
input <- vegdist(vsd, method="manhattan")
mod <- betadisper(input, group = factor2, type = "median")
mod

oneByTwo <- paste(factor1,factor2,sep=".")
conditions <- data.frame(cbind(factor1,factor2,oneByTwo))
adonis2(t(vsd)~factor2*factor4*factor5, data = conditions, permutations = 1000000, method = "manhattan")

##test adonis order doesn't matter for my data:
adonis2(t(vsd)~factor4*factor5*factor2, data = conditions, permutations = 1000000, method = "manhattan") #almost exactly the same as above, don't need to worry about order!


#plot PC axis 1 and 2 for tp1
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (31.17%)", ylab = "PC2 (14.65%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))

oneByTwo <- paste(factor2,factor4,sep=".")
conditions <- data.frame(cbind(factor2,factor4,oneByTwo))
adonis2(t(vsd)~factor2*factor4, data = conditions, permutations = 1000000, method = "manhattan")

#plot PC axis 3 and 4 for tp1
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (11.79%)", ylab = "PC2 (9.46%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))


#plot PC1 and 2 for tp2
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(15, 16)[as.numeric(as.factor(factor4))], xlab = "PC1 (27.01%)", ylab = "PC2 (15.27%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))

oneByTwo <- paste(factor2,factor4,sep=".")
conditions <- data.frame(cbind(factor2,factor4,oneByTwo))
adonis2(t(vsd)~factor2+factor4, data = conditions, permutations = 1000000, method = "manhattan")
adonis2(t(vsd)~factor3, data = conditions, permutations = 1000000, method = "manhattan")

adonis2(t(vsd)~factor4, data = conditions, permutations = 1000000, method = "manhattan") ### for GOL tp2 only!

#plot PC axis 3 and 4 for tp2
plot(scores[,3], scores[,4], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 1)[as.numeric(as.factor(factor4))], xlab = "PC1 (12.72%)", ylab = "PC2 (9.19%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))


#plot PC1 and 2 for tp2 with phenotpyic data
PhenData <- read.table("../../lipid_protein/lipid_protein_data_forPCA.txt", header = TRUE, row.names = "library") #import and attach phenotypic data
PhenColData <- read.delim("../../lipid_protein/lipid_protein_coldata_forPCA.txt", header=TRUE, row.names=1)

vare.pca <- rda(PhenData, scale = TRUE)
vare.pca
fit2 <- envfit(vare.pca, PhenColData, perm = 999)
fit2
plot(fit2)
plot(vare.pca, scaling = 3)
biplot(vare.pca, scaling = 3)
dev.off()

plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(15, 16)[as.numeric(as.factor(factor4))], xlab = "PC1 (27.01%)", ylab = "PC2 (15.27%)")
ordispider(scores,factor3, col=c("blue", "blue", "orange", "orange"))
plot(fit2)
plot(vare.pca, scaling = 3)

ord2 <- cca(PhenData ~ protein + lipid + phospholipid + sterol + fattyacid + triacylglycerol + waxester + respiration, PhenColData)
#plot(ord2, type="p")
fit <- envfit(ord2, coldata, perm = 999, display = "lc")
fit
plot(fit, col = "red")

