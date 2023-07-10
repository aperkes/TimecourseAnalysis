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

samples <- read.csv('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName
head(samples[,c("Treatment",'Time','Batch')])

data_dir = "/data/sequencing/TimeFights/results/Salmon_quant"
data_dir = "~/Documents/Data/Salmon_quant"

files <-file.path(data_dir,samples$SampleName)
files <- paste(sep="",files,".quant/quant.sf")
names(files) <- samples$SampleName

head(files)
#tx2gene <- read_csv(tx2gene.gencode.csv)
#txi <- tximport(files,type= "salmon",tx2gene=tx2gene)
txi <- tximport(files,type= "salmon",txOut = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Time + Batch)
dds <- DESeq(ddsTxi)
res <- results(dds,contrast = c("Treatment","W","L"))

head(res)
resOrdered <- res[order(res$padj),]
head(res[order(res$padj),])

sum(res$padj<0.05,na.rm=TRUE)

## Plot MA of all points
plotMA(res)

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

shared_genes_i <- res_justA$padj < 0.05 & !res_justF$padj < 0.05 & !is.na(res_justA$padj) & !is.na(res_justF$padj)

genes <- rownames(res)
## Note the '!'
unique_Ai <- res_justA$padj < 0.05 & !is.na(res_justA$padj) & !shared_genes
unique_Fi <- res_justF$padj < 0.05 & !is.na(res_justF$padj) & !shared_genes

## Save all these to csv's so I can inspect gene id's
write.csv(genes[unique_Ai],'unique_A.csv',row.names = FALSE)
write.csv(genes[unique_Fi],'unique_F.csv',row.names = FALSE)
write.csv(genes[shared_genes_i],'unique_A.csv',row.names = FALSE)
write.csv(rownames(res_sig),'all_sig.csv',row.names = FALSE)
