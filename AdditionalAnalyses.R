library("biomaRt")
#load packages we will need
library('DESeq2')
library('enrichplot')
library('Biobase')
library('DESeq')
library('edgeR')
library('genefilter')
library('gplots')
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("biomaRt")
library('KEGGgraph')
library('pathview')

ol <- read.csv("overlap_HOM_WT_FDR05_use.csv", header=TRUE,row.names=1)
autism <- read.csv("neurogenes.csv", header=TRUE, row.names=1)

autism.ol <- ol[rownames(ol)%in%rownames(autism),]
#write.csv(autism.ol, file="autism_deg.csv")

listMarts()
#choose ensembl dataset and filters to use
ensembl=useMart("ensembl")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = 'www.ensembl.org',
                  ensemblRedirect = FALSE)
ensembl = useEnsembl(biomart = "ensembl",dataset="mmusculus_gene_ensembl")
#ensembl = useEnsembl(biomart = "ensembl",dataset="mmusculus_gene_ensembl",mirror="useast")

filters = listFilters(ensembl)

filters[1:5,]
filterOptions("ensembl_gene_id",ensembl)

annotations_all<-getBM(attributes=c('ensembl_gene_id','entrezgene_id','description','external_gene_name','chromosome_name','start_position','end_position','strand'), filters='ensembl_gene_id', values=rownames(ol), mart=ensembl)
idmap_all = merge(x = data.frame(ol), y = annotations_all, by.x="row.names",by.y='ensembl_gene_id',all.x=TRUE)
#write.csv(idmap_all,file="annotated_union_degs.csv")


#pca plot and clustering for all samples together
#load full dataset in
setwd("/depot/pccr/data/Yangyang/hr02688_SCN2A-RNASeq/RNASeq_pipeline_2/output/counts/Ranalysis_rep2/")
load("yang_July2019RNAseq.RData")

setwd("/depot/pccr/data/Yangyang/H202SC19122900_RNASeq/Analysis/2_RNASeq_pipeline/output/DE_pipeline/")

newres <- read.table('newMat.txt',sep="\t",header=TRUE)
mat <- cbind(rawdata[,c(1:9,11:15)],newres[,2:9])
samples <- data.frame(row.names=c("HETbrain1","HOMbrain1","WTbrain1","X240CC1","X240HC2","X242CC2","X242HC1","X244CC1","X244HC1","HET716","HOM717","WT718","HOM739","WT742","newWT1","newWT2","newWT3","newWT4","newHOM1","newHOM2","newHOM3","newHOM4"),
                      genotype=as.factor(c("HET","HOM","WT",rep("HOM", 2),rep("HET",2),rep("WT",2),"HET","HOM","WT","HOM","WT",rep("WT",4),rep("WT",4))),
                     batch=as.factor(c(rep(1,9),rep(3,3),rep(2,2),rep(4,8))),
                      sample=c("HETbrain1","HOMbrain1","WTbrain1","X240CC1","X240HC2","X242CC2","X242HC1","X244CC1","X244HC1","HET716","HOM717","WT718","HOM739","WT742","newWT1","newWT2","newWT3","newWT4","newHOM1","newHOM2","newHOM3","newHOM4"),
                      sex=as.factor(c(rep("male",9),rep("female",3),rep("male",10))),
                      tissueType=as.factor(c(rep("wholeBrain",3), "cortex","hipp","cortex","hipp","cortex","hipp",rep("wholeBrain",13))))
data<-newCountDataSet(mat,samples)
countData<-counts(data)
colData<-pData(data)[,c("genotype","sex","sample","batch")]
dds<-DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design =~ genotype)
dds<- DESeq(dds)
res2 <- results(dds,contrast=c("genotype","HOM","WT"),alpha=0.05)
summary(res2)
res2["ENSMUSG00000075318",]

#rld is preferable if size factors vary a lot, and some of them do.  so ill use this
rld <- rlog(dds)

#generate a heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sample)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(filename="allSamples_sampleHeatmap.png",width=6,height=4,res=300,unit="in")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
png(filename="allSamples_sampleHeatmap_colorBatch.png",width=7,height=4,res=300,unit="in")
DESeq2::plotPCA(rld, intgroup=c("batch"))
dev.off()
DESeq2::plotPCA(rld, intgroup=c("sex"), ellipse=TRUE)
png(filename="allSamples_sampleHeatmap_colorGenotype.png",width=7,height=4,res=300,unit="in")
DESeq2::plotPCA(rld, intgroup=c("genotype"))
dev.off()

#pathway diagrams

#get count matrix of all stats
load('DESeq2_results.RData')
dds.df <- data.frame(dds_result)
#get entrez ids
annotations_dds<-getBM(attributes=c('ensembl_gene_id','entrezgene_id','description','external_gene_name','chromosome_name','start_position','end_position','strand'), filters='ensembl_gene_id', values=rownames(dds.df), mart=ensembl)
idmap_dds = merge(x = data.frame(dds.df), y = annotations_dds, by.x="row.names",by.y='ensembl_gene_id',all.x=TRUE)
#write.csv(idmap_dds,file="annotated_deseq2_all.csv")
idmap_dds <- read.csv("annotated_deseq2_all.csv")

deg <- idmap_dds[idmap_dds$Row.names %in% row.names(ol),]
deg <- na.omit(deg)
#write.csv(deg,file="annotated_deseq2_degs.csv")
idmap_dds.1 <- idmap_dds[!is.na(idmap_dds$external_gene_name),]
idmap_dds.1[idmap_dds.1$external_gene_name=="Grm2",9] <-"108068"
idmap_dds.keep <- idmap_dds.1[!is.na(idmap_dds.1$entrezgene_id),]
idmap_dds.keep$kegg <- translateGeneID2KEGGID(idmap_dds.keep$entrezgene_id,organism = "mmu")
idmap_dds.keep <- idmap_dds.keep[!duplicated(idmap_dds.keep$entrezgene_id),]
rownames(idmap_dds.keep) <- idmap_dds.keep$entrezgene_id
idmap_dds.keep.deg <- idmap_dds.keep[idmap_dds.keep$Row.names %in% row.names(ol),]
idmap_dds.keep.deg <- idmap_dds.keep[idmap_dds.keep$padj < 0.05,]
use.dds <- (idmap_dds.keep.deg$log2FoldChange)
names(use.dds) <- idmap_dds.keep.deg$entrezgene_id

#use fold changes to color graphs HOM/WT
#pathways: 
#oxytocin signaling pathway 04921
#Glutamatergic synapse 04724
#GABAergic synapse 04727

dds.out <- pathview(gene.data=use.dds,pathway.id = "04921",species="mmu",
                     out.suffix="oxytocin_graph_DE",kegg.native=F,gene.idtype = "kegg")

dds.out <- pathview(gene.data=use.dds,pathway.id = "04727",species="mmu",
                    out.suffix="GABA_graph_DE",kegg.native=F,gene.idtype = "kegg")
dds.out <- pathview(gene.data=use.dds,pathway.id = "04727",species="mmu",
                    out.suffix="GABA_native_DE",kegg.native=T,gene.idtype = "kegg")

dds.out <- pathview(gene.data=use.dds,pathway.id = "04724",species="mmu",node.sum = "sum",low = "dodgerblue", high="yellow",
                    out.suffix="deseq_glutamatergic_graph_DE_lBlYl",kegg.native=F,gene.idtype = "kegg")
dds.out <- pathview(gene.data=use.dds,pathway.id = "04727",species="mmu",low = "dodgerblue", high="yellow",
                    out.suffix="deseq_GABA_graph_DE_lBlYl",kegg.native=F,gene.idtype = "kegg")
dds.out <- pathview(gene.data=use.dds,pathway.id = "04921",species="mmu",low = "dodgerblue", high="yellow",
                    out.suffix="deseq_oxytocin_graph_DE_lBlYl",kegg.native=F,gene.idtype = "kegg")
dds.out <- pathview(gene.data=use.dds,pathway.id = "04713",species="mmu",low = "dodgerblue", high="yellow",
                    out.suffix="deseq_circadian_graph_DE_lBlYl",kegg.native=F,gene.idtype = "kegg")


####### Barchart of KEGG pathways #######################
#pathway analysis diagram

pathways <- read.csv("../DE_pipeline/CP.KEGG.csv",header=TRUE)[,c(2,6)]
colnames(pathways) <- c("Pathway", "p")

pathways <- pathways[1:20,]
pathways$p.adjust <- -1*log10(pathways$p)

p<-ggplot(data=pathways, aes(x=reorder(Pathway,p.adjust), y=p.adjust)) +
  ylab("-log(padj)") + xlab("Pathways")+
  geom_bar(colour="black", fill="white",stat="identity", width = 0.5)+ coord_flip()+
  theme(axis.text = element_text(size = 28, face = "bold")) +
  theme(axis.title = element_text(size = 30, face = "bold"))

png(file="HOMvsWT_top20Pathways.png", width = 16000, height = 7000, bg = "white", pointsize=6, res=600)
p
dev.off()

############################################################
#  Heatmap of Oxytocin genes that Muriel sent
############################################################
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

oxy <- read.csv("Oxytocin_pathway.csv")
cpm <- read.table("../counts/DE_analysis/Yang_Yang_2020_RNAseq_CPM.txt")
cpm.use <- cpm[rownames(cpm) %in% oxy$ID,]
idmap = merge(x = data.frame(cpm.use), y = oxy, by.x="row.names",by.y='ID',all.x=TRUE)
idmap.use <- idmap[,2:10]
colnames(idmap.use) <- c("HOM1","HOM2","HOM3","HOM4","WT1","WT2","WT3","WT4","Gene.ID")



HM<-idmap.use
a = round(nrow(HM)/20)
hg = a * 800

# assigne new data frame X
x <- HM[,1:ncol(HM)-1]
x<-data.frame(log2(x+1))
row.names(x) <- HM$Gene.ID
newnames <- lapply(rownames(x), function(x) bquote(italic(.(x))))

hmcols<-colorRampPalette(c("dodgerblue","white","yellow"))(256)

png(filename = "oxytocin_heatmap_rowScale_CPM_clusterRows.png", width = 1000, height = 1800, pointsize = 1,  bg = "white", res =300)
pheatmap(x, show_rownames = T, scale="row",cluster_cols = F,color=hmcols,labels_row=as.expression(newnames))
dev.off()

############################################################
#  GO term analysis for membrane potential
###########################################################

#OBA0000099
#load genes annotated with this

oba <- read.table("musmusculus_OBO0000099.txt")

deg.oba <- deg[deg$external_gene_name %in% oba$V2,]
idmap.dds.oba <- idmap_dds[idmap_dds$external_gene_name %in% oba$V2,]
write.csv(deg.oba, file="membrane_potential_DEGs.csv")
write.csv(deg.oba, file="membrane_potential_allGenes_HOM_WT.csv")

############################################################
#  Split enrichment plot
############################################################
library('DOSE')
library('enrichplot')
library('clusterProfiler')
library('ggplot2')
library('dplyr')
library('stringr')
library('graphite')
library('forcats')

deg$kegg <- translateGeneID2KEGGID(deg$entrezgene_id,organism = "mmu")

original <- deg$log2FoldChange
names(original) <- deg$entrezgene_id
geneList <- na.omit(original)
geneList <- sort(geneList,decreasing=TRUE)
geneList.neg <- geneList[geneList< 0]
geneList.pos <- geneList[geneList> 0]

kegg.all <-enrichKEGG(gene=names(geneList),
                    organism="mmu",keyType="kegg",
                    pvalueCutoff=0.10)

kegg.p <-enrichKEGG(gene=names(geneList.pos),
             organism="mmu",keyType="kegg",
             pvalueCutoff=0.10)
kegg.n <-enrichKEGG(gene=names(geneList.neg),
                    organism="mmu",keyType="kegg",
                    pvalueCutoff=0.10)
kegg.test <- kegg.p
gse <-gseKEGG(gene=geneList,
                 organism="mmu",keyType="kegg",nPerm=10000,minGSSize = 3,maxGSSize = 800,
                 pvalueCutoff=0.1, pAdjustMethod = 'none',use_internal_data = FALSE,by="DOSE")
gse.f <-gseKEGG(gene=geneList,
              organism="mmu",keyType="kegg",nPerm=10000,minGSSize = 3,maxGSSize = 800,
              pvalueCutoff=0.1, pAdjustMethod = 'none',use_internal_data = FALSE,by="fgsea")
gse.dose <- gse
dotplot(gse.f,showCategory=50, title="Enriched Pathways",split=".sign")+facet_grid(.~.sign)
png(filename="GSEA_kegg.split.png",res=600,units="in",height=10,width=12)
dotplot(gse,showCategory=50, title="Enriched Pathways",split=".sign")+facet_grid(.~.sign)
dev.off()
png(filename="enrichment_kegg.pos.png",res=600,units="in",height=5,width=5)
dotplot(kegg.p,showCategory=50, title="Upregulated")
dev.off()
png(filename="enrichment_kegg.all_top10.png",res=600,units="in",height=5,width=8)
dotplot(kegg.all,showCategory=20, title="Enriched KEGG Pathways")
dev.off()
full.kegg <- kegg.all@result
full.kk <- setReadable(kegg.all,'org.Mm.eg.db',keytype="ENTREZID")
write.csv(full.kk,file="intersection_kegg.csv")
categories <- c("Circadian entrainment","Neuroactive ligand-receptor interaction",
                "GABAergic synapse","Calcium signaling pathway","cGMP-PKG signaling pathway",
                "Cushing syndrome","Oxytocin signaling pathway" ,"Salivary secretion",
                "GABAergic synapse","Aldosterone synthesis and secretion",
                "Retrograde endocannabinoid signaling")

keep.kegg<- full.kegg[full.kegg$Description %in% categories,]
keep.all.kegg <- kegg.all
keep.all.kegg@result <- keep.kegg
png(filename="enrichment_kegg.all_top10.png",res=600,units="in",height=5,width=8)
dotplot(keep.all.kegg,showCategory=20, title="Enriched KEGG Pathways")
dev.off()

png(filename="enrichment_kegg.neg.png",res=600,units="in",height=6,width=10)
dotplot(kegg.n,showCategory=50, title="Downregulated")
dev.off()

p <-dotplot(kegg.p,showCategory=50, title="Upregulated")
n <-dotplot(kegg.n,showCategory=50, title="Downregulated")
kp.df <- data.frame(kegg.p@result)
kp.df <- kp.df[kp.df$p.adjust <0.1,]
kp.df$NES <- rep(1,dim(kp.df)[1])
kn.df <- data.frame(kegg.n@result)
kn.df <- kp.df[kpn.df$.adjust <0.1,]
kn.df$NES <- rep(-1,dim(kn.df)[1])
kegg.df <- rbind(kn.df,kp.df)
kegg.test@result <- kegg.df
merged.res <- kegg.df
merged.res$type = "upregulated"
merged.res$type[merged.res$NES <0] = "downregulated"
merged.res.show <- merged.res[30,]
#p <- ggplot(merged.res, aes(x = GeneRatio, y=fct_reorder(Description, GeneRatio)))+
  geom_point(aes(size=GeneRatio,color=p.adjust))+theme_bw(base_size=14)+
  scale_color_gradient(limits=c(0,0.10),low="red")+ylab(NULL)+ggtitle("KEGG Enrichment")
#p + facet_grid(.~type)
kegg.test@result <- merged.res
png(filename="enrichment_kegg.split.png",res=600,units="in",height=7,width=10)
dotplot(kegg.test,showCategory=50,split="type")+facet_grid(.~type)
dev.off()

#make a volcano plot with most points as black, but with up as yellow and down as dodger blue
int <- read.csv("intersection_DEGs.csv")
threshold <- idmap_dds$Row.names %in% int$Gene_ID
length(which(threshold))
idmap_dds$threshold <- threshold
idmap_dds$threshold[which(idmap_dds$threshold == "TRUE" & idmap_dds$log2FoldChange > 0)] <- "Up"
idmap_dds$threshold[which(idmap_dds$threshold == "TRUE" & idmap_dds$log2FoldChange < 0)] <- "Down"
test <- idmap_dds
test[order(test$padj,decreasing=FALSE),]
idmap_dds2 <- idmap_dds
library(ggrepel)
idmap_dds2$genelabels <- ""
idmap_dds2$genelabels <- ifelse(idmap_dds2$external_gene_name %in% yangGenes, TRUE, FALSE)
idmap_dds3 <- idmap_dds
idmap_dds3$genelabels <- ""
idmap_dds3$genelabels <- ifelse(idmap_dds3$external_gene_name %in% yangGenes, TRUE, FALSE)

yangGenes <- c("Scn1a","Scn2a","Scn8a","Kcne2","Kcng4","Kcnv1","Kcna1","Kcna2","Kcnj10","Kcnk1")
png(filename="volcano_orange_v2.5_geneNames_11_5.png",res=600,units="in",height=7,width=10)
ggplot(idmap_dds2)+geom_point(aes(x=log2FoldChange,y=-log10(padj),colour=threshold),size=1)+
ggtitle("HOM vs WT") + xlab("log2 fold-change") + ylab("log10 adjusted p-value")+theme_bw(base_size=14) +
  theme(legend.position = "none", plot.title = element_text(size = rel(1.5), hjust = 0.5))+
  xlim(-3,2.5) + theme(axis.text = element_text(size=17), axis.title = element_text(size=rel(1.5)))+ geom_text_repel(aes(x=log2FoldChange,y=-log10(padj),label = ifelse(genelabels==T,as.character(idmap_dds2$external_gene_name),"")),segment.alpha=0.5,force=1,box.padding = unit(0.25,"lines"))+
  scale_color_manual(name="threshold",values = c("Up" = "orange","Down" = "dodgerblue","FALSE" = "darkgrey") )
dev.off()
