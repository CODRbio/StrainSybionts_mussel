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

group <- read.table( snakemake@input[[2]], sep="\t", header=F, row.names=1)
colnames(group)="Group"


all_samples=colnames(vcf@gt)[2:length(colnames(vcf@gt))]
gt_mat <- do.call(cbind,lapply(all_samples, function(x){
  y=vcf@gt[,x]
  tmp_f(y)
}))
colnames(gt_mat)=all_samples

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
library(ggrepel)

dat.dist <- vegdist(data2[,1:(length(data2)-1)],method='euclidean')
pcoa<- dudi.pco(dat.dist, scan = FALSE, nf=2)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

sample_site <- data.frame({pcoa$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2)) + geom_point(aes(colour = factor(data2$Group)),size=2)+
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) + geom_text_repel(aes(label=names))+ #可在这里修改点的透明度、大小
  scale_color_manual(values = brewer.pal(6,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))
pcoa_plot
ggsave(snakemake@output[[1]])

