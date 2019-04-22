

#' Calculate coefficient of variant per position
#'
#' @import GenomicRanges
#' @import data.table 
#' @import CopyNumberAnalysis
#'
#' @param depths (depths per position per sample)
#' @param selector (selector region)
#'
#' @return res (positions ordered by CV increasingly)
#'
#' @examples
#' calculate.invariable.pos(depths,selector)

calculate.invariable.pos <- function(depths,selector){
        alldepths <- intersectTarget(depths, selector)
        coordCols = c('Chromosome', 'Position')
        SampleNames <- setdiff(colnames(alldepths), coordCols)

        mat <- alldepths[, SampleNames, with = F]
        res = mat
        copy.mat = mat

        depthSpaceMeds <- apply(mat, 2, median, na.rm = TRUE)  ### Median depth of each sample
        medNormMat <- sweep(mat, 2, depthSpaceMeds, FUN = "/") ## devided the median of corresponding sample
        mat = data.table(medNormMat)

        res$Mean = rowMeans(copy.mat)
        res$SD = apply(copy.mat[,SampleNames, with = F],1,sd)
        res$CV = res$SD/res$Mean

        res$norm.Mean = rowMeans(mat)
        res$norm.SD = apply(mat[,SampleNames, with = F],1,sd)
        res$norm.CV = res$norm.SD/res$norm.Mean

        res$Chromosome = alldepths$Chromosome
        res$Position = alldepths$Position
	
	if (ncol(selector) == 5){
                depthsGr = GenomicRanges::GRanges(seqnames = Rle(alldepths$Chromosome), ranges = IRanges::IRanges(start = alldepths$Position, end = alldepths$Position))
                targetGr = GenomicRanges::GRanges(seqnames = Rle(selector$Chromosome), ranges = IRanges::IRanges(start = selector$Start, end = selector$End))

                hits <- findOverlaps(depthsGr, targetGr)
                mcols(hits) = hits
                hits = data.frame(mcols(hits))
                res$Gene = selector[hits[,2] ,V4]
        }

        res <- res[order(norm.CV,decreasing=F), ]
        return(data.frame(res))
}

#' GC correction of depths (per position per sample) 
#'
#' @import GenomicRanges
#' @import data.table 
#' @import CopyNumberAnalysis
#'
#' @param alldepths (depths per position per sample, can be single sample)
#' @param selector (selector regions)
#' @param testRegions (regions of genes to be tested)
#' @param gcfile (GC content of each selector region)
#' @param remove (whether exclude test genes when fitting GC curves, default is T) 
#'
#' @return alldepths (depths after GC correction)
#' @return fit (the fitted model for each sample)
#'
#' @examples
#' GC.correction(alldepths, selector, testRegions, gcfile, remove = T)

GC.correction <- function(alldepths, selector, testRegions, gcfile, remove = T){

        gc <- read.table(gcfile,row.names=1)
        coordCols = c('Chromosome', 'Position')
        SampleNames <- setdiff(colnames(alldepths), coordCols)

        alldepths <- intersectTarget(alldepths, selector) # positions overlapped with selector regions

        # get the GC content of positions
        depthsGr = GenomicRanges::GRanges(seqnames = Rle(alldepths$Chromosome), ranges = IRanges::IRanges(start = alldepths$Position, end = alldepths$Position))
        targetGr = GenomicRanges::GRanges(seqnames = Rle(selector$Chromosome), ranges = IRanges::IRanges(start = selector$Start, end = selector$End))

        hits <- findOverlaps(depthsGr, targetGr)
        mcols(hits) = hits
      
	regions = data.frame(mcols(hits))[,2]
        alldepths$region = regions
        alldepths$gc = gc[regions,3]

        # Remove test genes
        if (remove == T){
                testGr = GenomicRanges::GRanges(seqnames = Rle(testRegions$chromosome), ranges = IRanges::IRanges(start = testRegions$start, end = testRegions$end))
                hits <- findOverlaps(targetGr, testGr)
                mcols(hits) = hits
                alldepthsFilter = alldepths[alldepths$region %in% data.frame(mcols(hits))[,1] ==FALSE,]
        }else{
                alldepthsFilter = alldepths
        }

        Med <- apply(alldepthsFilter[,SampleNames,with=F], 2, median)  # Median depth per sample)

        # Median depth per selector region per sample
        res = apply(alldepthsFilter[,SampleNames, with=F], 2, function(x) aggregate(val ~ region, cbind(val = x, region = alldepthsFilter$region),median))
        df <- data.frame(matrix(unlist(res), ncol=length(res), byrow=F),stringsAsFactors=FALSE)
        normed.exons <- data.frame(df[(length(res[[1]]$region)+1):nrow(df),])
        rownames(normed.exons) = res[[1]]$region
        colnames(normed.exons) = SampleNames

	# Fit loess curve to the data
        fit <- apply(normed.exons, 2, function(x) loess(val ~ region, data.frame(cbind(val = x, region = gc[rownames(normed.exons),3])), family = "symmetric"))

        # Apply the curve to predict depths based on GC content
        pred.value = do.call(cbind, lapply(fit, function(x) predict(x, gc[regions,3])))

        mat = alldepths[,SampleNames,with=F]
        NormMat = do.call(cbind,lapply(1:ncol(mat), function(x) mat[,x,with=F]/pred.value[,x] * Med[x]))
        alldepths[, SampleNames := data.frame(NormMat), with=F]

        # Remove NA positions
        alldepths[, c("region","gc"):=NULL]
        alldepths = alldepths[complete.cases(alldepths),]
        return(list(alldepths = alldepths, fit = fit))
}


