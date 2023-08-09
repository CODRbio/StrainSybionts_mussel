#!/usr/bin/env Rscript
library(optparse)
library(pheatmap)

option_list <- list(
  make_option(c("-m", "--matrix"), type="character", default = "in.tsv", action = "store", help = "Data Matrix!"),
  make_option(c("-r", "--rowNormalization"), default = FALSE, action = "store_true", help = "Row Normalization!"),
  make_option(c("-R", "--rowCluster"), default = FALSE, action = "store_true", help = "Clustering on Row?"),
  make_option(c("-c", "--colNormalization"), default = FALSE, action = "store_true", help = "Column Normalization!"),
  make_option(c("-C", "--colCluster"), default = FALSE, action = "store_true", help = "Clustering on Column?"),
  make_option(c("-q", "--quantile"), default = FALSE, action = "store_true", help = "Sequential or quantile breaks?"),
  make_option(c("-p", "--Palette"), default = "colorBlueRed", action = "store", help = "Palette to use, colorBlueRed or colorWhiteRed?"),
  make_option(c("-o", "--output"), type = "character", default = "05.statistics/FstHeatmap.pdf",action = "store", help = "Output file name!")
)

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n),type=5)
  breaks[!duplicated(breaks)]
}


colorBlueRed = colorRampPalette(c("blue", "white", "red"))(20)
#mat_breaks <- quantile_breaks(, n = 20)
colorWhiteRed = colorRampPalette(c("white", "red"))(20)

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE, usage = "usage: %prog [options]"))

data <- read.table(file=opt$matrix, sep=",", header=T, row.names=1)
mat_breaks <- quantile_breaks(as.matrix(data), n = 20)
colorUse = get(opt$Palette)
if(opt$rowNormalization)
{
	scale="row"
} else if(opt$colNormalization)
{
	scale="column"
} else
{
	scale="none"
}

if(opt$quantile)
{
	pheatmap(data,cluster_rows=opt$rowCluster,cellheight=11, cellwidth=11, fontsize=10, display_numbers = F, number_format = "%.1f", cluster_cols=opt$colCluster,scale=scale,filename=opt$output,color=colorUse, breaks=mat_breaks)
} else
{
	pheatmap(data,cluster_rows=opt$rowCluster,cellheight=11, cellwidth=11, fontsize=10, display_numbers = F, number_format = "%.1f", cluster_cols=opt$colCluster,scale=scale,filename=opt$output,color=colorUse)
}
#colorRampPalette(c("blue", "white", "red"))(200))
