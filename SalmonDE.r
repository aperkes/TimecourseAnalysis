## R script to process salmon data

#install.packages("BiocManager", repos = "https://cloud.r-project.org")

#library('BiocManager')
## For some reason this hasn't worked inside the script,
# Had to install via R terminal
#BiocManager::install('DESeq2')
#BiocManager::install('tximport')

#BiocManager::install('edgeR','limma','statmod','affycoretools','ReportingTools')

library('DESeq2')
library('tximport')
library('jsonlite')
#library('ggplot2')
library("edgeR")
#library("statmod")
#library('vegan')
#library('ape')


samples <- read.csv('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName
samples['BFF136','Group'] <- 'LF' ## There's a typo in one sample name

head(samples[,c("Treatment",'Time','Batch')])

script_dir = "~/Documents/Scripts/TimecourseAnalysis"
data_dir = "/data/sequencing/TimeFights/results/Salmon_quant"
data_dir = "~/Documents/Data/Salmon_quant"

files <-file.path(data_dir,samples$SampleName)
files <- paste(sep="",files,".quant/quant.sf")
names(files) <- samples$SampleName

head(files)
#tx2gene <- read_csv(tx2gene.gencode.csv)
#txi <- tximport(files,type= "salmon",tx2gene=tx2gene)
txi <- tximport(files,type= "salmon",txOut = TRUE)

## only keep 
countsPM.bool <- cpm(txi$counts) > 1
countsPM.num <- 1*countsPM.bool
group_sums <- rowsum(t(countsPM.num),samples$Group)

## What this is saying is that it's in at least 
##  four samples for every group. 
keep.group <- apply(group_sums,2,FUN=min) >= 4 

## Alternatively, I could use max to say that it's in at least 
##  4 samples for at least 1 group. 
keep.group <- apply(group_sums,2,FUN=max) >= 4 

#keep <- rowSums(cpm(ddsTxi_) > 1) >= 4
keep <- rowSums(cpm(txi$counts) > 1) >= 4 ## This approach leaves 26000


#ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Time + Batch)
ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Batch)

#ddsTxi <- ddsTxi_

ddsTxi <- ddsTxi_[keep,]
ddsTxi <- ddsTxi_[keep.group,]

dds <- DESeq(ddsTxi)
res <- results(dds,contrast = c("Treatment","W","L"))

head(res)
resOrdered <- res[order(res$padj),]
head(res[order(res$padj),])

sum(res$padj<0.05,na.rm=TRUE)

## Plot MA of all points
DESeq2::plotMA(res)

## Get just the genes that are significant
res_sig <- res[res$padj < 0.05 & !is.na(res$padj),]
dds_SigGenes <- dds[res$padj < 0.05 & !is.na(res$padj),]


### Plot heatmap of for fold change over time (that requires all those sup points)

ddsA <- dds_SigGenes[,samples$Time == 'A']
ddsB <- dds_SigGenes[,samples$Time == 'B']
ddsC <- dds_SigGenes[,samples$Time == 'G']
ddsD <- dds_SigGenes[,samples$Time == 'D']
ddsE <- dds_SigGenes[,samples$Time == 'E']
ddsF <- dds_SigGenes[,samples$Time == 'F']

resA <- results(ddsA,contrast = c("Treatment","W","L"))
resB <- results(ddsB,contrast = c("Treatment","W","L"))
resC <- results(ddsC,contrast = c("Treatment","W","L"))
resD <- results(ddsD,contrast = c("Treatment","W","L"))
resE <- results(ddsE,contrast = c("Treatment","W","L"))
resF <- results(ddsF,contrast = c("Treatment","W","L"))


#ddsAF <- dds_SigGenes[,samples$Time == 'A' | samples$Time == 'F']
dds_S <- dds_SigGenes[,samples$Treatment == 'S']
res_null <- results(dds_S,contrast = c("Time","A","G"))
sum(abs(res_null$log2FoldChange))
                    
## Get order that will split negative from positive
log_order <- order(abs(res_sig$log2FoldChange),decreasing = TRUE)
resOrdered <- res_sig[log_order,]
i_pos <- resOrdered$log2FoldChange > 0
i_neg <- resOrdered$log2FoldChange < 0
resOrdered <- rbind(resOrdered[i_pos,],resOrdered[i_neg,])

## Luckily, I can index by row name
resA <- resA[rownames(resOrdered),]
resB <- resB[rownames(resOrdered),]
resC <- resC[rownames(resOrdered),]
resD <- resD[rownames(resOrdered),]
resE <- resE[rownames(resOrdered),]
resF <- resF[rownames(resOrdered),]

diff_matrix <- matrix(c(resA$log2FoldChange,
                        resB$log2FoldChange,
                        resC$log2FoldChange,
                        resD$log2FoldChange,
                        resE$log2FoldChange,
                        resF$log2FoldChange),ncol = 6)

## Import the stuff, note, this might break the above code. 
library('RColorBrewer')
library('gplots')

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

#heatmap.2(m_extrema, col=coul)
heatmap.2(diff_matrix,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)

## Sum fold expression to see when the biggest contrast is. 
colSums(abs(diff_matrix))

## Same heatmap, but I want to sort it by A or F: 
diff_matrixA <- diff_matrix[order(diff_matrix[,1],decreasing = TRUE),]
heatmap.2(diff_matrixA,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)

diff_matrixF <- diff_matrix[order(diff_matrix[,6],decreasing = TRUE),]
heatmap.2(diff_matrixF,dendrogram = 'none', Rowv = FALSE,Colv = FALSE, col=coul)

# Get the list of genes for A, F, and see how which overlap
tableA <- ddsTxi[,samples$Time == 'A']
tableF <- ddsTxi[,samples$Time == 'F']

## Need to reassign the design now that I only have one time each
design(tableA) <- ~ Treatment + Batch 
design(tableF) <- ~ Treatment + Batch

dds_justA <- DESeq(tableA)
dds_justF <- DESeq(tableF)

res_justA <- results(dds_justA,contrast=c("Treatment","W","L"))
res_justF <- results(dds_justF,contrast=c("Treatment","W","L"))

sig_ja <- res_justA[res_justA$padj < 0.05 & !is.na(res_justA$padj),]
sig_jf <- res_justF[res_justF$padj < 0.05 & !is.na(res_justF$padj),]

shared_genes_i <- res_justA$padj < 0.05 & res_justF$padj < 0.05 & !is.na(res_justA$padj) & !is.na(res_justF$padj)


genes <- rownames(res)
## Note the '!'
unique_Ai <- res_justA$padj < 0.05 & !is.na(res_justA$padj) & !shared_genes_i
unique_Fi <- res_justF$padj < 0.05 & !is.na(res_justF$padj) & !shared_genes_i

res_uniqueA <- res_justA[unique_Ai,]
res_uniqueF <- res_justF[unique_Fi,]
res_shared <- res_justA[shared_genes_i,]



## Save all these to csv's so I can inspect gene id's
write.csv(rownames(res_uniqueA[order(res_uniqueA$padj),]),
          file.path(script_dir,'unique_A.csv'),row.names = FALSE)
write.csv(rownames(res_uniqueF[order(res_uniqueF$padj),]),
          file.path(script_dir,'unique_F.csv'),row.names = FALSE)
write.csv(rownames(res_shared[order(res_shared$padj),]),
          file.path(script_dir,'shared_AF.csv'),row.names = FALSE)
write.csv(rownames(res_sig[order(res_sig$padj),]),
          file.path(script_dir,'all_sig.csv'),row.names = FALSE)


## On to the PCA! 
library('vegan')
library('ape')

dds_23 <- dds
#dds_23 <- dds[res$padj < 0.1 & !is.na(res$padj)]
dds_23 <- dds_23[,dds$Batch != 1]
#rld <- rlogTransformation(dds) # Not sure whether there are advantages

dds_ <- dds
dds_ <- dds_23

rld <- vst(dds_) # Use if there are many samples (>30), but breaks if <1000
#rld <- varianceStabilizingTransformation((dds_23))
#rld <- rlogTransformation(dds_23)


# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
#factor1 <- as.character(colData(dds)$PlannedTreatment) 
factor2 <- as.character(colData(dds_)$Treatment) 
factor3 <- as.character(colData(dds_)$Batch) 
factor4 <- as.character(colData(dds_)$Time) 

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
plot(scores[,1], scores[,2], col=c("blue","red","orange")[as.numeric(as.factor(factor3))], pch=c(17, 4, 16)[as.numeric(as.factor(factor2))], xlab = xLab, ylab = yLab)
ordispider(scores,factor3, col=c("blue","red","orange"))

## For treatment
plot(scores[,1], scores[,2], 
     col=c("blue","grey","gold")[as.numeric(as.factor(factor2))], 
     pch=c(17, 4, 16)[as.numeric(as.factor(factor2))], 
     #xlim = c(-0.1,0.1),
     #ylim = c(-0.1,0.1),
     xlab = xLab, ylab = yLab)
ordispider(scores,factor2, col=c("blue","grey","gold"))

## For time
plot(scores[,1], scores[,2], 
     col=c("lightblue","blue","darkblue","gold","orange","red")[as.numeric(as.factor(factor4))], 
     pch=c(17, 4, 16)[as.numeric(as.factor(factor2))], 
     #xlim = c(-0.1,0.1),
     #ylim = c(-0.1,0.1),
     xlab = xLab, ylab = yLab)
ordispider(scores,factor4, col=c("lightblue","blue","darkblue","gold","orange","red"))
