#Load the libraries
library(DESeq2)
library(ggplot2)
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(edgeR)
library(VennDiagram)

## Accept the arguments and store in the variables
args <- commandArgs(trailingOnly = TRUE)

countfile     = args[1]
control       = args[2]
treatment     = args[3]
control_rep   = args[4]
treatment_rep = args[5]
path          = args[6]

## Print the version and session-information
#header = paste(rep("#", 50), collapse = "")
header = ""

sink(file = "sessioninfo.txt")
cat(paste(header, "#Version Information", header, sep = '\n'))
cat('\n')
version
cat('\n')

cat(paste(header, "#Session Information", header, sep = '\n'))
cat('\n')
sessionInfo()
sink()

#Save DESeq2 log to file
sink(file = "DESeq2_log.txt")

##read in the data file
raw.data <- read.table(countfile, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")

if(names(raw.data)[1] != "Gene_ID") {
  stop("The genes column in counts file must be named as:Gene_ID")
}

#read the entire count matrix
rlen = length(raw.data)
count_matrix <- raw.data[2:rlen]
rownames(count_matrix) <- raw.data$Gene_ID

cat(paste(header, "#Counts file:", header, sep = '\n'))
cat('\n')
head(count_matrix)
cat('\n')

# Remove rows with zero counts
count_matrix$rowsum <- rowSums(count_matrix)
counts_nozero <- subset(count_matrix, count_matrix$rowsum >0)

#Check the dimensions of the original and non-zero counts matrix
cat(paste(header, "#Data Dimensions:", header, sep = '\n'))
cat("\n","Complete counts matrix dimensions","\n")
dim(count_matrix)
cat("\n","Nonzero counts matrix dimensions","\n")
dim(counts_nozero)

#Remove the rowsums column
countData <- subset(counts_nozero, select = -c(rowsum))

#get CPM matrix
counts_nozero_CPM<-cpm(countData)
#write.table(counts_nozero_CPM, file = "counts_CPM.txt", sep = "\t")

#Get log CPM matrix
counts_nozero_logCPM <- cpm(countData, log = TRUE)
#write.table(counts_nozero_logCPM, file = "counts_log_CPM.txt", sep = "\t", quote = FALSE)

#Replace 0 counts with 1
countData <- replace(countData, countData == 0, 1)

###################################################################################
#                           Prepare the DESeq2 Experiment design
###################################################################################

group1 = rep(control,control_rep)
group2 = rep(treatment,treatment_rep)

ExpDesign1 <- data.frame(row.names=colnames(countData),
                         condition=factor(c(group1, group2))   )

cat(paste(header, "#Experiment Design Matrix:", header, sep = '\n'))
cat('\n')
ExpDesign1
cat('\n')

dds <- DESeqDataSetFromMatrix(countData,
                              colData= ExpDesign1,
                              design= ~condition)

dds <- DESeq(dds)

cat(paste(header, "#DESeq2 Sizefactors:", header, sep = '\n'))
cat('\n')
sizeFactors(dds)
cat('\n')
cat(paste(header, "#DESeq2 ResultsNames:", header, sep = '\n'))
cat('\n')
resultsNames(dds)
cat('\n')


dds_result <- results(dds, contrast=c("condition",treatment, control) )
save(dds_result, file = "DESeq2_results.RData")

###################################################################################
#                           Unsupervised Data Exploration
###################################################################################
#rld is preferable if size factors vary a lot and mostly they do.
rld <- rlog(dds, blind= FALSE)

###################################################################################
#                           heatmap of sample-to-sample distances
###################################################################################

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cat(paste(header, "#Preparing Sample_distance_matrix.png", header, sep = '\n'))
cat('\n')

png(filename ="Sample_distance_matrix.png",
    width = 2000, height = 2000, pointsize = 6,  bg = "white", res = 300)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         show_rownames = TRUE,
         show_colnames = TRUE)

dev.off()

###################################################################################
#                           Cooks distance boxplot
###################################################################################

cat(paste(header, "#Preparing DESeq2_Cooks_distance_boxplot.png", header, sep = '\n'))
cat('\n')

png(filename ="DESeq2_Cooks_distance_boxplot.png",
    width = 2000, height = 2000, pointsize = 12,  bg = "white", res = 300)

par(mfrow=c(1,1))

par(cex.axis=1)

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, ylab="Cook's Distance Values")

dev.off()

###################################################################################
#                           DESeq2 MA plot
###################################################################################

cat(paste(header, "#Preparing DESeq2_MA.png", header, sep = '\n'))
cat('\n')

png(filename ="DESeq2_MA.png",
    width = 2000, height = 2000, pointsize = 12,  bg = "white", res = 300)

DESeq2::plotMA(dds_result ,main="DESeq2 MA", ylim=c(-2,2) )

dev.off()

###################################################################################
#                           DESeq2 cluster Dendrogram
###################################################################################

cat(paste(header, "#Preparing DESeq2_Cluster_Dendrogram.png", header, sep = '\n'))
cat('\n')

png(filename ="DESeq2_Cluster_Dendrogram.png",
    width = 2000, height = 2000, pointsize = 12,  bg = "white", res = 300)

dists1 <-dist(t(assay(rld)))
plot(hclust(dists1))

dev.off()

###################################################################################
#                           DESeq2 PCA with Label
###################################################################################

cat(paste(header, "#Preparing DESeq2_PCA_with_Labels.png", header, sep = '\n'))
cat('\n')

png(filename ="DESeq2_PCA_with_Labels.png",
    width = 2000, height = 2000, pointsize = 8,  bg = "white", res = 300)

data1 <- plotPCA(rld, intgroup="condition", returnData=TRUE)

percentVar1 <- round(100 * attr(data1, "percentVar"))

p1 <- ggplot(data1, aes(PC1, PC2,color= condition, theme(legend.position = "none"),
                    label = rownames(data1))) +
                    geom_point(size=1) +
                    xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar1[2],"% variance"))

p1 + geom_point() + geom_text(hjust = 0.5, vjust = 1, nudge_x = 0.5, size = 3) + 
  theme(legend.position = "none")

dev.off()


###################################################################################
#                           Parsing the DESeq2 results
###################################################################################

#load(file = "DESeq2_results.RData")
dds_result <- as.data.frame(dds_result)
dds_result = rownames_to_column(dds_result, var = "Gene_ID")

dds_result <- mutate(dds_result, FoldChange = 2^log2FoldChange)
dds_result <- select(dds_result, "Gene_ID", "baseMean", "FoldChange", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
dds_result <- arrange(dds_result, padj)

###################################################################################
#                           Merging the Annotations
###################################################################################

cat(paste(header, "#Reading Annotations", header, sep = '\n'))
cat('\n')
cat(paste("Annotations were read from the file:", path, sep = " "))
cat('\n')
cat("DESeq2 results table dimensions before merging: \n")
dim(dds_result)
cat('\n')

Anno <- read.table(path, sep = "\t", stringsAsFactors = F, header = TRUE, quote = "")

if(names(Anno)[1] != "Gene_ID") {
  stop("The genes column in Annotation file must be named as: Gene_ID")
}

dds_result <- left_join(dds_result, Anno, by = "Gene_ID")

cat("DESeq2 results table dimensions After merging: \n")
dim(dds_result)
cat('\n')


###################################################################################
#                           Filter DESeq2 results
###################################################################################

dds_result_FDR005 <- filter(dds_result, padj <= 0.05)
dds_result_FDR001 <- filter(dds_result, padj <= 0.01)

cat(paste(header, "#Data Dimensions", header, sep = '\n'))
cat('\n')
cat("Dimension for the entire DESeq2 results table: \n")
dim(dds_result)
cat('\n')
cat("Dimension for the DESeq2 results table at FDR cutoff 0.05: \n")
dim(dds_result_FDR005)
cat('\n')
cat("Dimension for the DESeq2 results table at FDR cutoff 0.01: \n")
dim(dds_result_FDR001)
cat('\n')

###################################################################################
#                           Write to result tables
###################################################################################

cat(paste(header, "#Writing Results", header, sep = '\n'))
cat('\n')
cat("Complete DESeq2 Results: DESeq2_All.tsv \n")
cat("Filtered results with 95% FDR confidence: DESeq2_FDR005_filtered.tsv \n")
cat("Filtered results with 99% FDR confidence: DESeq2_FDR001_filtered.tsv \n")
cat('\n')

write.table(dds_result, file="DESeq2_All.tsv", quote=FALSE, sep='\t', row.names = FALSE)
write.table(dds_result_FDR005, file="DESeq2_FDR005_filtered.tsv", quote=FALSE, sep='\t', row.names = FALSE)
write.table(dds_result_FDR001, file="DESeq2_FDR001_filtered.tsv", quote=FALSE, sep='\t', row.names = FALSE)

##########################################################################################
#  Prepare file for GSEA
#########################################################################################
cat(paste(header, "#Preparing file for GSEA:", header, sep = '\n'))
cat('\n')

GSEA = select(dds_result, "Gene_ID", "Gene.name", "log2FoldChange", "pvalue")
GSEA$Rank = sign(GSEA$log2FoldChange) * -log10(GSEA$pvalue)
GSEA_export = select(GSEA, "Gene_ID", "Gene.name", "Rank")
GSEA_export = dplyr::arrange(GSEA_export, desc(Rank)) %>% drop_na()
GSEA_export = select(GSEA_export, "Gene_ID", "Rank")
names(GSEA_export) = c("Gene_ID", "Rank")
write.table(GSEA_export, file = "GSEA.rnk", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Pre-ranked input for GSEA: GSEA.rnk \n")
cat("Dimensions of GSEA.rnk file are: \n")
dim(GSEA_export)
cat('\n')

sink()

##########################################################################################
#  Differential Gene Expresssion with edgeR
#########################################################################################

#Save DESeq2 log to file
sink(file = "edgeR_log.txt")

cat(paste(header, "#Library Sizes", header, sep = '\n'))
cat('\n')
#calculate library sizes for each replicate
library.sizes = colSums(countData[1:ncol(countData)])
library.sizes
cat('\n')

y1 = DGEList(counts=as.matrix(countData), lib.size=library.sizes, group=ExpDesign1$condition)
y1 = calcNormFactors(y1)
y1 = estimateCommonDisp(y1)
y1 = estimateTagwiseDisp(y1)
edgeR = exactTest(y1, pair=c(control, treatment)) # Exact test

edgeR_results = topTags(edgeR, n=Inf)
edgeR_results = edgeR_results$table
edgeR_results$fc <- 2^edgeR_results$logFC

cat(paste(header, "#edgeR Results Table", header, sep = '\n'))
cat('\n')
head(edgeR_results)
cat('\n')
save(edgeR_results, file = "edgeR_results.RData")
###################################################################################
#                           Merging edgeR results and Annotations
###################################################################################

edgeR_results = rownames_to_column(edgeR_results, var = "Gene_ID")

cat(paste(header, "#Merging Annotations", header, sep = '\n'))
cat('\n')
cat(paste("Annotations were read from the file:", path, sep = " "))
cat('\n')
cat("edgeR results table dimensions before merging: \n")
dim(edgeR_results)
cat('\n')

edgeR_results <- left_join(edgeR_results, Anno, by = "Gene_ID")

cat("edgeR results table dimensions After merging: \n")
dim(edgeR_results)
cat('\n')

edgeR_results <- arrange(edgeR_results, FDR)

edgeR_result_FDR005 <- dplyr::filter(edgeR_results, FDR <= 0.05)
edgeR_result_FDR001 <- dplyr::filter(edgeR_results, FDR <= 0.01)


cat(paste(header, "#Writing Results", header, sep = '\n'))
cat('\n')
cat("Complete edgeR Results: edgeR_All.tsv \n")
cat("Dimensions of edgeR_All.tsv \n")
dim(edgeR_results)
cat("Filtered results with 95% FDR confidence: edgeR_FDR005_filtered.tsv \n")
cat("Dimensions of edgeR_FDR005_filtered.tsv \n")
dim(edgeR_result_FDR005)
cat("Filtered results with 99% FDR confidence: edgeR_FDR001_filtered.tsv \n")
cat("Dimensions of edgeR_FDR001_filtered.tsv \n")
dim(edgeR_result_FDR001)
cat('\n')

write.table(edgeR_results, file="edgeR_All.tsv", quote=FALSE, sep='\t', row.names = FALSE)
write.table(edgeR_result_FDR005, file="edgeR_FDR005_filtered.tsv", quote=FALSE, sep='\t', row.names = FALSE)
write.table(edgeR_result_FDR001, file="edgeR_FDR001_filtered.tsv", quote=FALSE, sep='\t', row.names = FALSE)

###################################################################################
#                           edgeR Visualizations
###################################################################################

#smear plot
cat(paste(header, "#Writing smear plot: edgeR_smear.png", header, sep = '\n'))
cat('\n')

png(file="edgeR_smear.png", width = 1200, height = 1200, bg = "white",pointsize=10, res=300)
summary(de1 <- decideTestsDGE(edgeR),p.value=0.05)
detags1 <- rownames(y1)[as.logical(de1)]
plotSmear(edgeR, de.tags=detags1)
abline(h=c(-2, 2), col="blue")
dev.off()


##mds plot
cat(paste(header, "#Writing mds plot: edgeR_mds.png", header, sep = '\n'))
cat('\n')
png(file="edgeR_mds.png" ,width = 2000, height = 2000, bg = "white", pointsize = 10, res=300)
plotMDS(y1, method="logFC")
dev.off()

