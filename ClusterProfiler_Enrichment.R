library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Cf.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(cowplot)

##R script accepts arguments
args <- commandArgs(trailingOnly = TRUE)

filepath      = args[1]
prefix        = args[2]
organism      = args[3]
pval          = as.numeric(args[4])
qval          = as.numeric(args[5])


#define organism of interest
if(organism == "Human") {db = "org.Hs.eg.db"}
if(organism == "Mouse") {db = "org.Mm.eg.db"}
if(organism == "Dog")   {db = "org.Cf.eg.db"}
if(organism == "Yeast") {db = "org.Sc.sgd.db"}

if(organism == "Human") {kegg_code = "hsa"}
if(organism == "Mouse") {kegg_code = "mmu"}
if(organism == "Dog")   {kegg_code = "cfa"}
if(organism == "Yeast")   {kegg_code = "sce"}


x = read.table(filepath, header = TRUE, stringsAsFactors = FALSE)

# Fetch ENTREZID for the given ENSEMBL IDs
x <- bitr(x$Gene_ID, fromType="ENSEMBL", toType="ENTREZID", OrgDb=db)

# #Define Go Categories
CATEGORY <- c("BP", "CC", "MF")

plot_list = list()

for (c in CATEGORY)
{
  ego2 <- enrichGO(gene         = x$ENTREZID,
                   OrgDb         = db,
                   ont           = c,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = pval,
                   qvalueCutoff  = qval,
                   readable = TRUE)


  outfig <- paste(prefix, c , "pdf", sep = ".")
  outcsv <- paste(prefix,c,"csv", sep = ".")
  title <- paste("GO", c, "Enrichment")

  bp <- barplot(ego2, showCategory=10) + ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  dp <- dotplot(ego2, showCategory=10) + ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

  plot_list[[c]] <- dp

  write.csv(ego2, file=outcsv)

  pdf(outfig,width=16,height=10,paper='special')
  print(bp)
  print(dp)

  dev.off()

}


#kegg Enrichment
title = "KEGG Enrichment"
outcsv = paste(prefix,"KEGG","csv", sep = ".")
outfig <- paste(prefix, "KEGG" , "pdf", sep = ".")

kk <- enrichKEGG(gene         = x$ENTREZID,
                 organism     = kegg_code,
                 pvalueCutoff = pval)

kk_readable = setReadable(kk, OrgDb = db, keytype ="ENTREZID")

bp <- barplot(kk, showCategory= 10) + ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dp <- dotplot(kk, showCategory= 10) + ggtitle(title) + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

plot_list[["KEGG"]] <- dp

write.csv(kk_readable, file=outcsv, row.names = FALSE)

pdf(outfig,width=14,height=10,paper='special')

print(bp)
print(dp)

dev.off()

combine_outfile <- paste(prefix, "combined" , "pdf", sep = ".")

pdf(file = combine_outfile, width=20,height=20,paper='special')
print(plot_grid(plotlist = plot_list))
dev.off()

