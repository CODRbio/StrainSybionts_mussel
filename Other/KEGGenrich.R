#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-k", "--kegg"), type="character", action = "store", help = "KEGG tables gene2term (pathway) tsv file"),
  make_option(c("-g", "--gene2ko"), type="character", action = "store", help = "gene2ko (kegg ortholog) tsv file"),
  make_option(c("-s", "--selected"), type="character", action = "store", help = "selected gene lists for enrichment analysis"),
  make_option(c("-o", "--output"), type="character", action = "store", default = "KEGG", help = "Output prefix name, default KEGG"),
  make_option(c("-p", "--padj"), type="double", default = 0.05, action = "store", help = "Expected padj value for the DEGs, default 0.05"),
  make_option(c("-P", "--pvalue"), default = FALSE, action = "store_true", help = "Use pvalue rather than padj, default fause"),
  make_option(c("-c", "--count"), type="character", action = "store", default = "NA", help = "Count table to generate DEGs, NA means equal threshold for each genes to represent KO, default NA"),
  make_option(c("-b", "--background"), type="character", action = "store", default = "NA", help = "Genes of the backgrounds, default all the annotations"),
  make_option(c("-S", "--startC"), type="integer", action = "store", default = 6, help = "start Column of the Count table to generate DEGs, default 6 for feature count, 1 for stringtie")
  )

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE, usage = "usage: %prog [options]"))


library("pathview")
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(ggupset)
library(dplyr)

keggAnn <- read.table("/nfs_genome/BioDB/kegg/202204/kegg/pathway/map_title.tab",header=F,sep="\t")

if(opt$background != "NA")
{
        bg <- read.delim(opt$background,header=F, sep="\t")
        Gene2kegg <- read.table(opt$kegg,header=F,sep="\t") %>% filter(.[[1]] %in% bg[[1]])
}else
{
Gene2kegg <- read.table(opt$kegg,header=F,sep="\t")
}


if(opt$count != "NA")
{
	count<-read.table(opt$count, header=T,sep="\t", row.names=1)
	threshold<-log(rowSums(count[,opt$startC:dim(count)[2]]),2)%/%1
} else
{
	threshold<-rep(1,dim(Gene2kegg)[1])
	names(threshold)<-Gene2kegg[,1]
}

diff <- read.table(opt$selected,header=T,sep="\t",row.names=1)

if(opt$pvalue)
{
        num = 5
	x =enricher(rownames(diff),TERM2GENE=Gene2kegg[,c(2,1)],TERM2NAME=data.frame(name=sprintf("%s%05d","ko",keggAnn$V1),def=keggAnn$V2),pvalueCutoff=1, qvalueCutoff=1)
} else
{
        num = 6
	x =enricher(rownames(diff),TERM2GENE=Gene2kegg[,c(2,1)],TERM2NAME=data.frame(name=sprintf("%s%05d","ko",keggAnn$V1),def=keggAnn$V2),pvalueCutoff=opt$padj, qvalueCutoff=(opt$padj*4))
}
enrichRes <- x[x[,num]<=opt$padj]

num_sig <- dim(enrichRes)[1]

head(enrichRes)
write.table(enrichRes, file = paste0(opt$output,"_keggEnriched.tsv"), quote = FALSE, sep="\t")


x2 <- pairwise_termsim(x)
if(num_sig >= 50)
{
        num_show = 50
} else
{
        num_show = num_sig
}

library(RColorBrewer)
col <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))
print("Statistics charts, dotplot")
pdf(paste0(opt$output,"_dotplot.pdf"),width=15,height=30)
dotplot(x,showCategory=num_show)
dev.off()
print("emetplot")
pdf(paste0(opt$output,"_emetplot.pdf"),width=20,height=20)
emapplot(x2,showCategory=num_show)
dev.off()
print("cnetplot")
pdf(paste0(opt$output,"_cnetplot.pdf"),width=20,height=20)


fc<-diff$log2FoldChange
#try(names(fc)<-rownames(diff))
if(is.null(fc)){
fc <- rep(1, dim(diff)[1])
names(fc)<-rownames(diff)
} else {
names(fc)<-rownames(diff)
switch = 1

}



fc2 <- fc
fc2[fc2 > 5] = 5
fc2[fc2 < -5] = -5


#if !(is.null(fc)){
cnetplot(x2,foldChange=fc2,showCategory=num_show,node_label="all") + scale_color_gradientn(colours=col)
dev.off()
pdf(paste0(opt$output,"_heatplot.pdf"),width=30,height=12)

print("heatplot")
heatplot(x2,foldChange=as.numeric(fc2),showCategory=num_show) #+ scale_fill_gradientn(colours=col)
dev.off()
pdf(paste0(opt$output,"_upsetplot.pdf"),width=30,height=12)

print("upsetplot")
upsetplot(x2,n=num_show)
dev.off()
print("Drawing pathways enriched")
library(pathview)
library(hash)
source("/nfs_genome/BioInformatics/perls/pathview-u.R")
koAnn <- read.table(opt$gene2ko,header=F,sep="\t")
gene2ko <- hash(keys=koAnn$V1, values=koAnn$V2)
genKO <- function(id){if (is.null(gene2ko[[id]])) {return()} else {return(gene2ko[[id]])}}
genTH <- function(id){if (is.null(threshold[id])) {return(1)} else {return(threshold[id])}}
genFC <- function(id){if (is.null(fc[id])) {return(1)} else {return(fc[id])}}
KOs <- levels(as.factor(koAnn$V2))
fcKO <- rep(0,length(KOs))
names(fcKO) <-KOs

for (num in seq(dim(enrichRes)[1])) 
{ 
path<-enrichRes[num,1]; 
enGenes<-unlist(strsplit(enrichRes[num,8],"/")); 
enKOs<-sapply(enGenes,genKO)

enFCs<-sapply(split(names(enKOs),enKOs), function(x) {weighted.mean(genFC(x),genTH(x))})
fcKO[names(enFCs)] <- enFCs

pv.out <- pathview_u(gene.data = fcKO, pathway.id = path ,species = "ko", out.suffix = opt$output, same.layer = T, low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd = "yellow"), na.col = "transparent", kegg.native = T, kegg.dir="/nfs_genome/BioDB/kegg/202204/kegg/xml/kgml/ko")
pv.out <- pathview(gene.data = fcKO, pathway.id = path ,species = "ko", out.suffix = opt$output, same.layer = T, kegg.native = F, kegg.dir="/nfs_genome/BioDB/kegg/202204/kegg/xml/kgml/ko")

}

