
doMC::registerDoMC(cores = 1)

#' Overload functions.
overloard <- function (FUN, ...) {
    .FUN <- FUN
    args <- list(...)
    invisible(lapply(seq_along(args), function(i) {
        formals(.FUN)[[names(args)[i]]] <<- args[[i]]
    }))
    .FUN
}

.write.table <- overloard(write.table, sep = "\t", quote = F, row.names = F)

#' Write TSV with consistent configuration for usage in RSU
#'
#' @param x Data frame to be written to file
#' @param filename File name
#'
#' @return
#' @export
#'
#' @examples
write.rsu.tsv <- function(x, filename, ...) {
    write.table(format(x, trim=TRUE, digits=10), file = filename, sep = "\t", quote = F, row.names = F, ...)
}

#' Read TSV with consistent configuration for usage in RSU
#'
#' @param filename File name
#'
#' @return Data table with contents in filename
#' @export
#'
#' @examples
read.rsu.tsv <- function(filename) {
    data.table(read.table(filename, header = T, stringsAsFactors = F, sep = "\t"))
}


#' Load depths from a list of frequency (*.freq) files. Removes characters after "Tumor" or "cfDNA" from the sample names
#'
#' @param collapsedCountFiles List of frequency files
#'
#' @param numCores Number of cores to load depths
#'
#' @return Data frame of depths for all positions by samples
#' @export
#'
#' @import data.table
#' @import plyr
loadRSUReadDepths <- function(collapsedCountFiles, numCores = NA) {
    
  if (!is.na(numCores))
    doMC::registerDoMC(cores = numCores)

    loadReadDepth <- function(collapsedCountFile) {
        if (!file.exists(collapsedCountFile)) {
            return(NA)
        }
        collapsedCount = fread(collapsedCountFile)
        collapsedCount = collapsedCount[order(CHR, POS), list(CHR, POS, DEPTH)]
        sampleName = sub("_(Tumor|cfDNA).*", "", sub("Sample_", "", basename(collapsedCountFile)))
        setnames(collapsedCount, c("CHR", "POS", "DEPTH"), c("Chromosome", "Position", 
            sampleName))
        collapsedCount
        
    }
    
    # names(collapsedCountFiles) = collapsedCountFiles
    
    collapsedCountFiles = collapsedCountFiles[!is.na(collapsedCountFiles)]
    collapsedCounts = llply(collapsedCountFiles, loadReadDepth, .parallel = T)
    
    annotatedCounts = collapsedCounts[[1]]
    for (cc in collapsedCounts[-1]) {
        annotatedCounts = merge(annotatedCounts, cc, by = c("Chromosome", "Position"))
    }
    
    annotatedCounts
}

#'  Load RSU selector/panel bed file
#'
#' @param selectorFile Selector/panel bed file
#'
#' @import data.table
loadRSUSelector <- function(selectorFile) {
    selector <- data.table(read.delim(selectorFile, as.is = T, header = F))
    setnames(selector, c("V1", "V2", "V3"), c("Chromosome", "Start", "End"))
    #selector <- selector[order(Chromosome, Start), ]
    selector[, `:=`(Length, End - Start + 1)]
    # selector <- selector[order(selector$Contig, selector$RegionStart), ] selector$Gene
    # = sub('gene_name=(\\S+);.*', '\\1', selector$FirstNote, perl=T) selector[ ,
    # c('Left', 'Mid', 'Right'):=NULL, with=F] selector[ grepl('[:;]', Gene), Gene:=NA ]
    # selector[ , FirstNote:= NULL ]
    selector
}
