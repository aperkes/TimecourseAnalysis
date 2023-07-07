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

files <-file.path('/data/sequencing/TimeFights/results/Salmon_quant',samples$SampleName)
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
