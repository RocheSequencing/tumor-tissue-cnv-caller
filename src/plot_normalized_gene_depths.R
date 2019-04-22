#!/usr/bin/env Rscript
## Generate plot showing bin-level log2 coverages and segmentation calls.
##
## Usage:
##      Rscript scatter_plot.R --cnr sample.cnr --cns sample.cns --gene gene_name --png output.png

library(argparse)

options(bitmapType='cairo')

geneChoices=c('EGFR', 'ERBB2', 'MET')


parser <- ArgumentParser(description='Generate plot showing bin-level log2 coverages and segmentation calls.')
parser$add_argument('--cnr', default=NULL, required=T,
                    help='cnr file with bin-level log2 ratios.')
parser$add_argument('--cns', default=NULL, required=T,
                    help='cns file with segmented log2 ratios.')
parser$add_argument('--gene', default=NULL, required=T, choices=geneChoices,
                    help='gene of interest.')
parser$add_argument('--png', default=NULL, required=T,
                    help='name of the output file.')

opt <- parser$parse_args()

# default args for Arg
file_cnr <- opt$cnr
file_cns <- opt$cns
gene <- opt$gene
file_png <- opt$png

#####################################################

Tcnr <- read.table(file_cnr,sep="\t",stringsAsFactors=F,header=T)
Tcns <- read.table(file_cns,sep="\t",stringsAsFactors=F,header=T)

## use grep to account for cases: gene1,gene2
idx_gene <- grep(gene, Tcnr$gene)
gchr <- unique(Tcnr$chromosome[idx_gene])
gstart <- min(Tcnr$start[idx_gene])
gend <- max(Tcnr$end[idx_gene])
sname <- gsub(".cnr","",basename(file_cnr))

## whole chromosome
idx <- which(Tcnr$chromosome==gchr)

## gene coordinates in index coordinates
gstart_idx <- idx_gene[1]-idx[1]+1
gend_idx <- gstart_idx+length(idx_gene)-1

## define color palette
gcols <- c(rep("gray66",5),rep("gray36",3),rep("gray10",2))

## bins colors encode weights
bin_weight <- Tcnr$weight
cols <- rep(NA,idx[length(idx)])

cols[idx] <- gcols[round(10*bin_weight[idx]) + 1]
cols[idx_gene] <- "red"

png(file_png, width=500, height=450, pointsize=16)
plot(Tcnr$log2[idx], col=cols[idx], cex=0.4, pch=19,
     main=sname,
     xlab=paste0("Bins on ",gchr),
     ylab="Copy ratio (log2)")

## Bin-by-bin identify corresponding segment values
bin_log2 <- Tcnr$log2[idx]
bin_chr <- Tcnr$chromosome[idx]
bin_start <- Tcnr$start[idx]
bin_end <- Tcnr$end[idx]

seg_log2 <- rep(NA,length(bin_log2))
for(i in seq_along(bin_log2)){
  seg_log2[i] <- Tcns$log2[which(Tcns$chromosome==bin_chr[i] &
                                   Tcns$start<=bin_start[i] &
                                   Tcns$end>=bin_end[i])]
}

lines(seg_log2,col="darkorange",lwd=3)
abline(h = 0, col = "black")
abline(v = gstart_idx, col = "black", lwd=0.5)
abline(v = gend_idx, col = "black", lwd=0.5)

## add gene name to the plot
tpos_x <- 0.5*(gstart_idx+gend_idx)
tpos_y <- 0.9*max(Tcnr$log2[idx])
text(tpos_x, tpos_y, gene, col="black", pos=1)
dev.off()