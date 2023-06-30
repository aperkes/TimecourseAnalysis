
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
countdata_ <- countdata

countdata <- countdata[,coldata$Batch != 1]
coldata_full <- coldata
coldata <- coldata[coldata$Batch != 1,]



win_rows <- coldata$PlannedTreatment == 'win'
loss_rows <- coldata$PlannedTreatment == 'loss'
control_rows <- coldata$PlannedTreatment == 'control'

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

ddsFullCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_, design = ~ Batch + PlannedTreatment + Time) # Batch effects can't be random effect here
ddsInteractionCountTable <- DESeqDataSetFromMatrix(countData=countdata_, colData=coldata_, design = ~ Batch + PlannedTreatment*Time) # Batch effects can't be random effect here


ddsWinTable <- DESeqDataSetFromMatrix(countData=win_counts, colData=col_win, design = ~ Batch + Time) 
ddsLossTable <- DESeqDataSetFromMatrix(countData=loss_counts, colData=col_loss, design = ~ Batch + Time) 
ddsControlTable <- DESeqDataSetFromMatrix(countData=control_counts, colData=col_c, design = ~ Batch + Time) 

ddsFull <- DESeq(ddsFullCountTable) # this is the analysis!
head(ddsFull)

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
