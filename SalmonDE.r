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

samples <- read.table('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName
head(samples[,c("Treatment",'Time','Batch')])

files <-file.path('/data/sequencing/TimeFights/results/Salmon_quant/',samples$SampleName,'.quant')
names(files) <- samples$SampleName

#tx2gene <- read_csv(tx2gene.gencode.csv)
#txi <- tximport(files,type= "salmon",tx2gene=tx2gene)

