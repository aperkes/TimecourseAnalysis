
BiocManager::install('DESeq2','edgeR','limma','statmod','affycoretools','ReportingTools')

install.packages("pheatmap")
install.packages("vegan")
install.packages("ape")
install.packages("rgl")

#install.packages("gplots")
#install.packages("RColorBrewer")

library("edgeR")
library("statmod")
library('DESeq2')
library('vegan')
library('ape')

#setwd("/data/sequencing/TimeFights/analysis")
setwd("~/Documents/Scripts/TimecourseAnalysis/")

countdata_1 <- read.delim("all_quant_.sf", header=TRUE, row.names="Name") #read in the file of the count data and call it countdata, row.names tells the name in the top left cell, the gene names

countdata_1 <- round(countdata_1)

keep <- rowSums(cpm(countdata_1)>1) >=4
# keeps rows (genes) where at least 4 columns (libraries) have at least 1 count per million. This means that if a gene is only expressed in say one treatment (which has three replicates), this gene will not be thrown out of the analysis
## This should be 4 per treatment, not across the whole set ###

countdata <- countdata_1[keep,] #formatting for organizing the kept rows that summed to at least 1 cpm in the step above

countdata <- round(countdata)

## Need to sort both so that they are same order
countdata <- countdata[,order(colnames(countdata))]

#coldata <- read.delim("./columnInfo.txt", header=TRUE, row.names=1) #use all_v3 for all data
coldata <- read.delim("./columnInfo.tsv", header=TRUE, row.names=1) #use tsv to drop A0

coldata <- coldata[order(rownames(coldata)),]

## optionally drop the first batch to remove batch effect
countdata_ <- countdata

countdata <- countdata[,coldata$Batch != 1]
coldata_full <- coldata
coldata <- coldata[coldata$Batch != 1,]



win_rows <- coldata$PlannedTreatment == 'win'
loss_rows <- coldata$PlannedTreatment == 'loss'
control_rows <- coldata$PlannedTreatment == 'control'

## Skip this if you're not using A0
t0_controls <- coldata$Time == 'A0' & coldata$PlannedTreatment == 'control'
## I actually need to grab the t0 control to use for DE
win_rows <- as.logical(win_rows + t0_controls)
loss_rows <- as.logical(loss_rows + t0_controls)

col_win <- coldata[win_rows,]
col_loss <- coldata[loss_rows,]
col_c <- coldata[control_rows,]

win_counts <- countdata[,win_rows]
loss_counts <- countdata[,loss_rows]
control_counts <- countdata[,control_rows]

## Here I can drop controls, since they are not quite balanced
#coldata <- coldata[-36:-68,]
#countdata <-countdata[,-36:-68]

#ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ POP*TREATMENT) #use POP*TREATMENT for tp 1 and 2, POP*TREATMENT*TP for all data

ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_full, design = ~ Batch + PlannedTreatment + Time) # Batch effects can't be random effect here
ddsInteractionCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_full, design = ~ Batch + PlannedTreatment + Time + PlannedTreatment:Time) # Batch effects can't be random effect here
ddsInteractionCountTable <- estimateSizeFactors(ddsInteractionCountTable)

ddsWinTable <- DESeqDataSetFromMatrix(countData=win_counts, colData=col_win, design = ~ Batch + Time) 
ddsLossTable <- DESeqDataSetFromMatrix(countData=loss_counts, colData=col_loss, design = ~ Batch + Time) 
ddsControlTable <- DESeqDataSetFromMatrix(countData=control_counts, colData=col_c, design = ~ Batch + Time) 

## Get WIN-LOSS contrasts for each time point:
A_rows <- coldata$Time == 'A' | coldata$Time == 'A0'
B_rows <- coldata$Time == 'B'
C_rows <- coldata$Time == 'G'
D_rows <- coldata$Time == 'D'
E_rows <- coldata$Time == 'E'
F_rows <- coldata$Time == 'F'

col_A <- coldata[A_rows,]
col_B <- coldata[B_rows,]
col_C <- coldata[C_rows,]
col_D <- coldata[D_rows,]
col_E <- coldata[E_rows,]
col_F <- coldata[F_rows,]

A_counts <- countdata[,A_rows]
B_counts <- countdata[,B_rows]
C_counts <- countdata[,C_rows]
D_counts <- countdata[,D_rows]
E_counts <- countdata[,E_rows]
F_counts <- countdata[,F_rows]

ddsATable <- DESeqDataSetFromMatrix(countData=A_counts, colData=col_A, design = ~ Batch + PlannedTreatment) 
ddsAshufTable <- DESeqDataSetFromMatrix(countData=A_counts, colData=col_A, design = ~ Batch + ShuffledTreatment) 

ddsBTable <- DESeqDataSetFromMatrix(countData=B_counts, colData=col_B, design = ~ Batch + PlannedTreatment) 
ddsCTable <- DESeqDataSetFromMatrix(countData=C_counts, colData=col_C, design = ~ Batch + PlannedTreatment) 
ddsDTable <- DESeqDataSetFromMatrix(countData=D_counts, colData=col_D, design = ~ Batch + PlannedTreatment) 
ddsETable <- DESeqDataSetFromMatrix(countData=E_counts, colData=col_E, design = ~ Batch + PlannedTreatment) 
ddsFTable <- DESeqDataSetFromMatrix(countData=F_counts, colData=col_F, design = ~ Batch + PlannedTreatment) 

ddsA <- DESeq(ddsATable)
ddsAshuf <- DESeq(ddsAshufTable)

ddsB <- DESeq(ddsBTable)
ddsC <- DESeq(ddsCTable)
ddsD <- DESeq(ddsDTable)
ddsE <- DESeq(ddsETable)
ddsF <- DESeq(ddsFTable)

res_A <- results(ddsA,contrast=c("PlannedTreatment","win","control"))
res_Ashuf <- results(ddsAshuf,contrast=c("ShuffledTreatment","win","control"))

res_B <- results(ddsB,contrast=c("PlannedTreatment","win","loss"))
res_C <- results(ddsC,contrast=c("PlannedTreatment","win","loss"))
res_D <- results(ddsD,contrast=c("PlannedTreatment","win","loss"))
res_E <- results(ddsE,contrast=c("PlannedTreatment","win","loss"))
res_F <- results(ddsF,contrast=c("PlannedTreatment","win","loss"))

rldA <- rlogTransformation(ddsA)
rldB <- rlogTransformation(ddsB)
rldC <- rlogTransformation(ddsC)
rldD <- rlogTransformation(ddsD)
rldE <- rlogTransformation(ddsE)
rldF <- rlogTransformation(ddsF)

normedA <- assay(rldA)
normedB <- assay(rldB)
normedC <- assay(rldC)
normedD <- assay(rldD)
normedE <- assay(rldE)
normedF <- assay(rldF)

A_diff <- rowMeans(normedA[,9:12]) - rowMeans(normedA[,1:4])
B_diff <- rowMeans(normedB[,9:12]) - rowMeans(normedB[,1:4])
C_diff <- rowMeans(normedC[,9:12]) - rowMeans(normedC[,1:4])
D_diff <- rowMeans(normedD[,9:12]) - rowMeans(normedD[,1:4])
E_diff <- rowMeans(normedE[,9:12]) - rowMeans(normedE[,1:4])
F_diff <- rowMeans(normedF[,9:12]) - rowMeans(normedF[,1:4])

diff_matrix <- matrix(c(A_diff,B_diff,C_diff,D_diff,E_diff,F_diff),ncol = 6)
mean_diffs <- rowMeans(diff_matrix)
head(normedF)

A_order <- order(A_diff,decreasing=TRUE)
F_order <- order(F_diff,decreasing = TRUE)
mean_order <-order(mean_diffs,decreasing = TRUE)

A_statorder <- order(res_A$padj)
F_statorder <- order(res_F$padj)
#All_statorder <- order(res$padj)
All_statorder <- order(res_full$padj)

m_sorted <- diff_matrix[A_order,]

## It works a little different if you're using the stat order
m_sorted <- head(diff_matrix[All_statorder,],n=100)
m_sorted <- head(diff_matrix[F_statorder,],n=100)
m_pos <- m_sorted[rowMeans(m_sorted) > 0,]
m_neg <- m_sorted[rowMeans(m_sorted) < 0,]
m_extrema <- rbind(head(m_pos,n=50),head(m_neg,n=50))

#m_sorted <- m_sorted[order(rowMeans(m_sorted),decreasing = TRUE),]
m_sorted <- m_sorted[order(m_sorted[,1],decreasing = TRUE),]


m_top = head(m_sorted,n=50)
m_bottom = tail(m_sorted,n=50)

m_extrema = rbind(m_top,m_bottom)

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

#heatmap.2(m_extrema, col=coul)
heatmap.2(m_extrema,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)


geneNames <- rownames(countdata)
A_genes <- head(geneNames[A_statorder],n=100)
F_genes <- head(geneNames[F_statorder],n=100)
X_genes <- head(geneNames[All_statorder],n=100)

AF_genes <- intersect(A_genes,F_genes)
AX_genes <- intersect(A_genes,X_genes)
FX_genes <- intersect(F_genes,X_genes)
## " End Here "

## Shuffled
ddsShuffledCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_full, design = ~ Batch + ShuffledTreatment + Time) # Batch effects can't be random effect here

## Real
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_full, design = ~ Batch + PlannedTreatment + Time) # Batch effects can't be random effect here

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)
res_full = results(ddsFull,contrast = c("PlannedTreatment","win","loss"))

#ddsShuffle <- DESeq(ddsShuffledCountTable)
#res_shuffled = results(ddsShuffle,contrast = c("ShuffledTreatment","win","loss"))
#res_full <- res_shuffled

rfullOrdered <- res_full[order(res_full$padj),]
head(rfullOrdered)
sum(res_full$padj<0.05,na.rm=TRUE)


ddsInteraction <-DESeq(ddsInteractionCountTable)
head(ddsInteraction)
res_int = results(ddsInteraction, tidy=TRUE, )
rintOrdered <- res_int[order(res_int$padj),]
head(rintOrdered)
sum(res_int$padj<0.05, na.rm=TRUE)


ddsWin <- DESeq(ddsWinTable)
ddsLoss <- DESeq(ddsLossTable)
ddsControl <- DESeq(ddsControlTable)

res_win <- results(ddsWin,contrast=c('Time','A','F'))
res_win
rwinOrdered <- res_win[order(res_win$padj),]
head(rwinOrdered)
sum(res_win$padj<0.05, na.rm=TRUE)

res <- results(ddsFull,contrast=c('PlannedTreatment','win','loss'))
res
resOrdered <- res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05, na.rm=TRUE)
sum(res$pvalue<0.05,na.rm=TRUE)


ddsSelect <- ddsWin
ddsSelect <- ddsLoss
ddsSelect <- ddsControl

res <- results(ddsSelect,contrast=c('Time','A','F'))
resOrdered <- res[order(res$padj),]
head(resOrdered)
sum(res$padj<0.05,na.rm=TRUE)




#rld <- rlogTransformation(ddsFull) # Not sure whether there are advantages
rld <- vst(ddsSelect) # Use if there are many samples (>30)
rld_control <- vst(ddsControl)
rld_win <- vst(ddsWin)

normed_counts = assay(rld)
normed_controls = assay(rlc_control)

## IF winner: 
means_X = rowMeans(normed_counts[,1:4])
means_A <- rowMeans(normed_counts[,5:8])
means_C <- rowMeans(normed_counts[,25:28]) ## I used G instead of C, it was dumb and I'm sorry
means_B <- rowMeans(normed_counts[,9:12])
means_D <- rowMeans(normed_counts[,13:16])
means_E <- rowMeans(normed_counts[,17:20])
means_F <- rowMeans(normed_counts[,21:24])

#IF loser: 
means_X = rowMeans(normed_counts[,25:28])

# If control:
means_X = rowMeans(normed_counts[,1:4])

means_A <- rowMeans(normed_counts[,1:4])
means_C <- rowMeans(normed_counts[,21:24]) ## I used G instead of C, it was dumb and I'm sorry
means_B <- rowMeans(normed_counts[,5:8])
means_D <- rowMeans(normed_counts[,9:12])
means_E <- rowMeans(normed_counts[,13:16])
means_F <- rowMeans(normed_counts[,17:20])


## There is absolutely an easier way to do this, 
## But I'm new to R and don't have time to figure out how to be clever today

diff_XX <- means_X - means_X
diff_AA <- means_A - means_X
diff_AB <- means_B - means_X
diff_AC <- means_C - means_X
diff_AD <- means_D - means_X
diff_AE <- means_E - means_X
diff_AF <- means_F - means_X

m <- matrix(c(diff_XX,diff_AA,diff_AB,diff_AC,diff_AD,diff_AE,diff_AF),ncol=7)

mean_all <- rowMeans(as.matrix(m))

order_by_col = order(m[,3])
order_by_mean = order(mean_all) ## What if we order it by win mean? 

m_sorted <- m[order(m[,3]),]
m_sorted <- m[order(mean_all),]

head(m_sorted)

#m_sorted <- matrix(as.numeric(m_sorted),ncol=ncol(m_sorted))

m_top = head(m_sorted,n=500)
m_bottom = tail(m_sorted,n=500)

m_extrema = rbind(m_top,m_bottom)

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

#heatmap.2(m_extrema, col=coul)
heatmap.2(m_extrema,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)

# assembling table of conditions to label PCoA plot:

# There's a big batch effect, so for the PCA to be useful, 
# probably best to drop the first batch. Or average treatments.

# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1 <- as.character(colData(ddsSelect)$PlannedTreatment) 
factor2 <- as.character(colData(ddsSelect)$Experience) 
factor3 <- as.character(colData(ddsSelect)$Batch) 
factor4 <- as.character(colData(ddsSelect)$Time) 

head(assay(rld))

hist(assay(rld))


library(gplots) 

library("RColorBrewer") 

#install.packages('genefilter')

BiocManager::install('genefilter')
library( "genefilter" )

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 40)



heatmap.2(assay(rld)[topVarGenes, ], scale="row",
          trace="none", dendrogram="both", key=TRUE, keysize = 1.5, margins =c(3,11), density.info = "density",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = heatmap.2( assay(rld)[ topVarGenes, ], # category labels
                           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),  # color key
                           lty= 1,             # line style
                           lwd = 10))          # line width
