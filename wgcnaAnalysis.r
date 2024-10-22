library("wgcna")
library("tidyr")

library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)

## A basic run through of WGCNA, from https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
## The start is copied from SalmonDE.r, if you're already running that, you can skip to line 73

samples <- read.csv('~/Documents/Scripts/TimecourseAnalysis/SampleInfo.csv',header=TRUE)
rownames(samples) <- samples$SampleName
#samples['BFF136','Group'] <- 'LF' ## There's are some typos in the sample names


samples$Batch <- as.factor(samples$Batch)
samples.23 <- samples[samples$Batch != 1,]
samples.23A <- samples[samples$Time == 'A',]
samples <- samples.23A
head(samples[,c("Treatment",'Time','Batch')])

script_dir = "~/Documents/Scripts/TimecourseAnalysis"
data_dir = "/data/sequencing/TimeFights/results/Salmon_quant"
#data_dir = "~/Documents/Data/Salmon_quant"

files <-file.path(data_dir,samples$SampleName)
files <- paste(sep="",files,".quant/quant.sf")
names(files) <- samples$SampleName

head(files)
#tx2gene <- read_csv(tx2gene.gencode.csv)
#txi <- tximport(files,type= "salmon",tx2gene=tx2gene)
txi <- tximport(files,type= "salmon",txOut = TRUE)


## filtering at the group level.
countsPM.bool <- cpm(txi$counts) > 1
countsPM.num <- 1*countsPM.bool
group_sums <- rowsum(t(countsPM.num),samples$Group)

## Check group sums: 
apply(group_sums,1,FUN=max)

## What this is saying is that it's in at least 
##  four samples for every group. 
keep.group <- apply(group_sums,2,FUN=min) >= 4 

## Alternatively, I could use max to say that it's in at least 
##  4 samples for at least 1 group. 
keep.group <- apply(group_sums,2,FUN=max) >= 4 

#keep <- rowSums(cpm(ddsTxi_) > 1) >= 4
## Keep genes with greater than 1 cpm for at least 4 samples. 
keep <- rowSums(cpm(txi$counts) > 1) >= 4 ## This approach leaves 26000

samples$Batch <- as.factor(samples$Batch)
ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Time + Batch)
ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Time)
#ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Batch)
#ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment)

ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Time + Batch)
ddsTxi_ <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Batch)

ddsTxi_.F <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Treatment + Batch)
ddsTxi.F <- ddsTxi_.F[,samples$Time == 'F']
dds.F <- DESeq(ddsTxi.F)
res.F <- results(dds,contrast = c("Treatment","W","S"))
DESeq2::plotMA(res.F)

ddsTxi <- ddsTxi_
ddsTxi <- ddsTxi_[keep,]
ddsTxi <- ddsTxi_[keep.group,]

dds <- DESeq(ddsTxi)

#######################
#### START HERE #######
#######################

### Some normalization and data inspection prior to WGCNA
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized[1:5,1:10]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

input_mat = t(expr_normalized)

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns

## This is the part where I copy a bunch of code I don't understand 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)

text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

## Not sure, but I don't think those plots look great...

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)


# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = samples$Group
# MEs0$treatment = samples$Treatment # I only get one of these

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


## This is probably where I stop, but here's the rest of the code: 

# pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")
