library(stringr)
library(vcfR)
library(ggplot2)
vcf <- read.vcfR( snakemake@input[[1]], verbose = FALSE )
tmp_f <- function(x){
  gt=str_split(x,':',simplify = T)[,1]
  ad=str_split(str_split(x,':',simplify = T)[,2],',',simplify = T)[,2]
  dp=str_split(x,':',simplify = T)[,3]
  return(gt)
  # return(data.frame(gt=gt,dp=as.numeric(dp),ad=as.numeric(ad)))
}

all_samples=colnames(vcf@gt)[2:length(colnames(vcf@gt))]
gt_mat <- do.call(cbind,lapply(all_samples, function(x){
  y=vcf@gt[,x]
  tmp_f(y)
}))
colnames(gt_mat)=all_samples
group <- read.table( snakemake@input[[2]], sep="\t", header=F, row.names=1) 
colnames(group)="Group"
dat=as.factor(gt_mat)
dat=as.numeric(dat)
dat=matrix(dat,nrow = nrow(gt_mat))
colnames(dat)=all_samples
dat=t(dat)
dat=as.data.frame(dat) 
data2 = merge(dat, group, by='row.names', all=F)
rownames(data2) <- data2$Row.names; data2$Row.names <- NULL

#pheatmap::pheatmap(dat,show_colnames = F)
print("Caculating PCA")
library("FactoMineR")
library("factoextra")
library(ade4)
library(magrittr)
library(RColorBrewer)
library(vegan)

dat.pca <- dudi.pca(df = data2[,1:(length(data2)-1)], scannf = FALSE, nf = 2)

fviz_pca_ind(dat.pca, repel =T,
             geom = c("point", "text"), # show points only (nbut not "text")
             col.ind = data2$Group, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave(snakemake@output[[1]])

