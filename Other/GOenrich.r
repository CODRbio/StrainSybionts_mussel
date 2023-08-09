#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("-g", "--gopropation"), type="character", action = "store", help = "Go propation annotation tables"),
  make_option(c("-s", "--selected"), type="character", action = "store", help = "selected gene lists for enrichment analysis"),
  make_option(c("-o", "--output"), type="character", action = "store", default = "GOenriched.tsv", help = "Output file, default GOenriched.tsv"),
  make_option(c("-p", "--padj"), type="double", default = 0.05, action = "store", help = "Expected padj value for the DEGs, default 0.05"),
  make_option(c("-P", "--pvalue"), default = FALSE, action = "store_true", help = "Use pvalue rather than padj, default fause"),
  make_option(c("-b", "--background"), default = "NA", action = "store", help = "Genes of the backgrounds, default all the annotations")  
)

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE, usage = "usage: %prog [options]"))

library(GO.db)
library(topGO)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(ggupset)
Gene2go <- read.delim(opt$gopropation,header=T,sep="\t")

if(opt$background != "NA")
{
	bg <- read.delim(opt$background,header=F, sep="\t")
	Gene2go <- read.delim(opt$gopropation,header=F,sep="\t") %>% filter(.[[1]] %in% bg[[1]]) 
}else
{

Gene2go <- read.delim(opt$gopropation,header=T,sep="\t")
}
goTer<-function(id){if (is.null(GOTERM[[id]])) {return("NA")} else {return(Term(GOTERM[[id]]))}}
goOnto<-function(id){if (is.null(GOTERM[[id]])) {return("NA")} else {return(Ontology(GOTERM[[id]]))}}

GOannotated<-levels(as.factor(Gene2go[,2]))
GOontolog<-sapply(GOannotated,goOnto)
GOter<-sapply(GOannotated,goTer)


diff <- read.table(opt$selected,header=T,sep="\t",row.names=1)

if(opt$pvalue)
{
        num = 5
	x =enricher(rownames(diff),TERM2GENE=Gene2go[,c(2,1)],TERM2NAME=data.frame(name=names(unlist(GOter)),def=unlist(GOter)),pvalueCutoff=1, qvalueCutoff=1)
} else
{
        num = 6
	x =enricher(rownames(diff),TERM2GENE=Gene2go[,c(2,1)],TERM2NAME=data.frame(name=names(unlist(GOter)),def=unlist(GOter)),pvalueCutoff=opt$padj, qvalueCutoff=(opt$padj*4))
}
enrichRes <- x[x[,num]<=opt$padj]

num_sig <- dim(enrichRes)[1]
goDEF<-function(id){if (is.null(GOTERM[[id]])) {return("NA")} else {return(Definition(GOTERM[[id]]))}}



GOs<-as.list(enrichRes$ID)[grepl("GO:",(as.list(enrichRes$ID)))]
enrichRes$term <- sapply(GOs, goTer)
enrichRes$Defination <- sapply(GOs, goDEF)
enrichRes$Ontolog <- sapply(GOs, goOnto)
head(enrichRes)
write.table(enrichRes, file = paste0(opt$output,"_goEnriched.tsv"), quote = FALSE, sep="\t")

x2 <- pairwise_termsim(x)
go_rich_BP <- x
go_rich_MF <- x
go_rich_CC <- x

cat_sig <- dim(x2[x2[,num]<=opt$padj])[1]
go_rich_BP@result <- go_rich_BP@result[names(subset(GOontolog, GOontolog=="BP")),]
go_rich_BP@ontology <- "BP"

go_rich_MF@result <- go_rich_MF@result[names(subset(GOontolog, GOontolog=="MF")),]
go_rich_MF@ontology <- "MF"

go_rich_CC@result <- go_rich_CC@result[names(subset(GOontolog, GOontolog=="CC")),]
go_rich_CC@ontology <- "CC"

if(cat_sig >= 150)
{
	catalog_show = 150
} else
{
	catalog_show = cat_sig
}

if(num_sig >= 50)
{
        num_show = 50
} else
{
        num_show = num_sig
}



num_BP<-length(which(grepl("GO:",go_rich_BP[go_rich_BP[,num]<=opt$padj]$ID))) #x2[x2[,num]<=opt$padj]
num_MF<-length(which(grepl("GO:",go_rich_MF[go_rich_MF[,num]<=opt$padj]$ID)))
num_CC<-length(which(grepl("GO:",go_rich_CC[go_rich_CC[,num]<=opt$padj]$ID)))
if(num_BP >= 50)
{
        num_show = 50
} else
{
        num_show = num_BP
}


pdf(paste0(opt$output,"_GOgraph.pdf"),width=30,height=15)
plotGOgraph(go_rich_BP,firstSigNodes=num_show,useInfo="def",.NO.CHAR = 30)


if(num_MF >= 50)
{
        num_show = 50
} else
{
        num_show = num_MF
}


plotGOgraph(go_rich_MF,firstSigNodes=num_show,useInfo="def",.NO.CHAR = 30)


if(num_CC >= 50)
{
        num_show = 50
} else
{
        num_show = num_CC
}

plotGOgraph(go_rich_CC,firstSigNodes=num_show,useInfo="def", .NO.CHAR = 30)
dev.off()




library(RColorBrewer)
col <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))
pdf(paste0(opt$output,"_dotplot.pdf"),width=15,height=30)
dotplot(x,showCategory=num_show)
dev.off()
pdf(paste0(opt$output,"_emetplot.pdf"),width=20,height=20)
emapplot(x2,showCategory=catalog_show)
dev.off()


fc<-diff$log2FoldChange

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


pdf(paste0(opt$output,"_cnetplot.pdf"),width=20,height=20)

cnetplot(x2,foldChange=fc2,showCategory=catalog_show,node_label="all")+scale_color_gradientn(colours=col)
dev.off()
pdf(paste0(opt$output,"_heatplot.pdf"),width=30,height=12)

heatplot(x2,foldChange=fc2,showCategory=catalog_show)+scale_fill_gradientn(colours=col)
dev.off()

if(num_sig >= 25)
{
        num_show = 25
} else
{
        num_show = num_sig
}

pdf(paste0(opt$output,"_upsetplot.pdf"),width=30,height=12)

upsetplot(x2,n=num_show)
dev.off()

