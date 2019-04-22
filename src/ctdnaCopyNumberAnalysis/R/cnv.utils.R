doParallel = T

#' Syntactic sugar to concatenate strings
#'
#' @param x first string
#' @param y second string
#'
#' @return string - concatenated string
#' @export
#'
"%+%" = function(x, y) {
  paste(x, y, sep = "")
}


#' Write final CNV caller outputs to file. 
#'
#' @param sampleName Name of sample
#' @param outdir Output directory
#' @param result Full result list to be written to file
#' @param filteredReport CNV positive calls to be written to file
#' @param report All (positive and negative) CNV calls to be written to file
#'
#' @return No return value
write.final.output <-
  function(sampleName,
           outdir,
           result,
           filteredReport,
           report) {
    rdsFile = file.path(outdir, paste0(sampleName, '.cnv-internal.rds'))
    cnvCallFile = file.path(outdir, paste0(sampleName, '.cnv-call.bed'))
    whitelistFile = file.path(outdir, paste0(sampleName, '.cnv-whitelist.bed'))
    saveRDS(result, file = rdsFile)
    
    filteredSamp = filteredReport[filteredReport$sample == sampleName, -1]
    reportSamp = report[report$sample == sampleName, -1]
    
    header = paste0(
      '# track name=',
      sampleName,
      ' type=bed
#chrom\tchromStart\tchromEnd\tname\tscore(copy number gain;p-value)'
    )
    
    write(header, cnvCallFile)
    write.rsu.tsv(filteredSamp,
                  cnvCallFile,
                  col.names = F,
                  append = T)
    header = paste0(header, '\tthreshold')
    write(header, whitelistFile)
    write.rsu.tsv(reportSamp,
                  whitelistFile,
                  col.names = F,
                  append = T)
  }

#' Top level function to analyze CNV from test freq file
#'
#' @param testFreqFile freq file (usually with *.freq suffix) from freq generation command
#' @param targetFile panel target file in bed format, the usual usage are taret file from a preload data release
#' @param testRegionFile File name of TSV file with genes to analyze for CNV and their gene coordinates
#' @param normalDepths Depths of normal reference samples
#' @param sampleNames Sample names to be analyzed for CNV
#' @param outdir Output directory
#' @param numCores Number of CPU cores for analysis
#' @param fraction Top fraction of positions used in normalDepths for use to analyze CNV
#' @param gcCorrection Boolean to use GC content correction in CNV
#' @param gcfile File containing GC content for each position
#' @param remove Boolean to exclude test genes from testRegionFile when performing GC correction
#'
#' @import plyr
#' @import tools
#'
#' @return string - concatenated string
#' @export
#'
analyzeCNV <-
  function(testFreqFile,
           targetFile,
           testRegionFile,
           normalDepths,
           sampleNames,
           outdir = '.',
           numCores = NA,
           fraction = 1,
           gcCorrection = F,
           gcfile,
           remove = T) {
    if (!is.na(numCores))
      doMC::registerDoMC(cores = numCores)
    
    checkFile(testFreqFile)
    checkFile(targetFile)
    checkFile(testRegionFile)
    
    coordCols = c('Chromosome', 'Position')
    readDepth = loadRSUReadDepths(testFreqFile)
    selector <- loadRSUSelector(targetFile)
    testRegions = read.rsu.tsv(testRegionFile)
    
    if (fraction != 1) {
      CVnormalDepths = calculate.invariable.pos(normalDepths, selector)
    }
    if (gcCorrection == T) {
      normalDepths = GC.correction(normalDepths, selector, testRegions, gcfile, remove = remove)[[1]]
      readDepth = GC.correction(readDepth, selector, testRegions, gcfile, remove = remove)[[1]]
    }
    if (fraction != 1) {
      normalDepths = merge(normalDepths,
                           data.table(CVnormalDepths[1:round(nrow(CVnormalDepths) * fraction), c("Chromosome", "Position")]),
                           by = c("Chromosome", "Position"))
    }
    
    testSampleName <- setdiff(colnames(readDepth), coordCols)
    normSampleNames <- setdiff(colnames(normalDepths), coordCols)
    allSampleNames <- c(normSampleNames, testSampleName)
    alldepths <- merge(normalDepths, readDepth, by = coordCols)
    alldepths <- intersectTarget(alldepths, selector)
    alldepth_nz <- data.frame(alldepths, check.names = F)
    
    statistics <-
      dlply(testRegions, ~ gene + chromosome + start + end + pval.threshold,
            function(x) {
              plot_one (
                alldepth_nz,
                testSampleName,
                normSampleNames,
                x$chromosome,
                x$start,
                x$end,
                testRegions = testRegions,
                pagetitle = "10ng input",
                props = NULL
              )
            }, .parallel = T)
    
    report <- ldply(statistics,
                    function(x) {
                      data.frame(right.sided = x$right.sided,
                                 left.sided = x$left.sided)
                    }, .parallel = T)
    
    report$sample = sampleNames
    report = report[, c('sample',
                        'chromosome',
                        'start',
                        'end',
                        'gene',
                        'right.sided',
                        'pval.threshold')]
    
    filteredReport = subset(report, pval.threshold > right.sided)
    filteredReport$cnv = rep('gain', nrow(filteredReport), stringsAsFactors = FALSE)
    
    filteredReport = filteredReport[, c('sample',
                                        'chromosome',
                                        'start',
                                        'end',
                                        'gene',
                                        'right.sided',
                                        'cnv')]
    
    result = list(statistics = statistics,
                  report = report,
                  filteredReport = filteredReport)
    
    # check existence of outdir
    if (!file.exists(outdir)) {
      dir.create(file.path(outdir))
    }
    
    
    
    invisible(lapply(sampleNames, function(x)
      write.final.output(x, outdir, result, filteredReport, report)))
  }


#' Replace all zeros in data frame with the machine minimum
#'
#' @param df data frame of numeric in all columns
#'
#' @return a data frame with zeros replaced with the machine minimum
#' @export
#'
#' @examples
replaceZeros <- function(df) {
  replaceZero <- function(x) {
    x[x == 0] = .Machine$double.xmin
    x
  }
  
  dt = data.table(df)
  dt[, colnames(dt) := lapply(.SD, replaceZero), .SDcols = colnames(dt)]
  data.frame(dt)
}

#' Return the log2-median normalized depth from a data frame of absolute depths 
#'
#' @param depths Data frame of depths, the column names should be the sample names
#' @param sampleNames Sample names of depths to normalize
#'
#' @return data frame of normalized depths
#' @export
#'
#' @examples
log2normalize <- function(depths, sampleNames) {
  mat <- depths[, sampleNames, with = F]
  depthSpaceMeds <- apply(mat, 2, median, na.rm = TRUE)
  medNormMat <- sweep(mat, 2, depthSpaceMeds, FUN = "/")
  zeroReplacedMedNormMat <- replaceZeros(medNormMat)
  logNormMat <- log2(zeroReplacedMedNormMat)
  
  # remove total depth per sample This is subtracting a mean of logs from the log value which is the same as scaling by the geometric mean Alternatively - reverse the
  # operations - log2 and then scale mat = log2(mat) mat = scale(mat,scale=F) Not sure this is necessary but let's subtract the means again just to make sure the sample
  # mean log2 is zero - this will make the cell line and control samples 'match up' better mat = sweep(mat,2,colMeans(mat),FUN='-')
  retDepth = data.table(depths, check.names = F)
  retDepth[, `:=`(sampleNames, logNormMat), with = F]
  retDepth
}

#' Intersect data frame of depths with data frame of target genomic regions
#'
#' @param depths Data frame of depths with Chromosome and Position columns
#' @param target Data frame of genomic regions with Chromosome, Start, and End columns
#'
#' @return Data table with only positions that intersect target regions
#' @export
#'
#' @import GenomicRanges
#' @import gplots
#'
#' @examples
intersectTarget <- function(depths, target) {
  depthsGr = GenomicRanges::GRanges(
    seqnames = Rle(depths$Chromosome),
    ranges = IRanges::IRanges(start = depths$Position, end = depths$Position)
  )
  mcols(depthsGr) = depths
  targetGr = GenomicRanges::GRanges(
    seqnames = Rle(target$Chromosome),
    ranges = IRanges::IRanges(start = target$Start, end = target$End)
  )
  
  depthsGr = subsetByOverlaps(depthsGr, targetGr)
  intDepths = data.table(data.frame(mcols(depthsGr), check.names = F))
  names(intDepths) = colnames(depths)
  intDepths
}

#' Expand a data frame of genomic regions to a data frame with rows for each position
#'
#' @param testRegions Data frame of genomic regions with chromosome, start, end, and gene columns
#'
#' @return Data frame with rows for each position
#'
#' @examples
as.testRegionsByBase <- function (testRegions) {
  testRegionsByBase = copy(testRegions)
  testRegionsByBase[, lengths := end - start + 1]
  testRegionsByBase = data.table(
    Chromosome = rep(testRegionsByBase$chromosome, testRegionsByBase$lengths),
    Position = unlist(apply(testRegionsByBase[, c('start', 'end'), with = F], 1, function(x) {
      seq(x[1], x[2])
    })),
    geneName = rep(testRegionsByBase$gene, testRegionsByBase$lengths)
  )
  testRegionsByBase
}

#' Exclude specific genes from a data frame of depths
#'
#' @param depths Data frame of depths for each position
#' @param testRegionsByBase Data frame of all positions for a gene
#' @param geneNameArg Genes to exclude from depths data frame
#'
#' @return Data frame with geneNameArg genes excluded
#' @export
#'
#' @examples
excludeByBase <-
  function (depths, testRegionsByBase, geneNameArg) {
    depthsGr = GenomicRanges::GRanges(
      seqnames = Rle(depths$Chromosome),
      ranges = IRanges::IRanges(start = depths$Position, end = depths$Position)
    )
    
    mcols(depthsGr) = depths[, !names(depths) %in% c('Chromosome', 'Position')]
    
    filteredRegion = testRegionsByBase[geneName != geneNameArg]
    testRanges = GenomicRanges::GRanges(
      seqnames = Rle(filteredRegion$Chromosome),
      ranges = IRanges::IRanges(start = filteredRegion$Position, end = filteredRegion$Position)
    )
    
    depthsGr = depthsGr[is.na(match(depthsGr, testRanges))]
    depthsFiltered = data.table(Chromosome = as.character(seqnames(depthsGr)),
                                Position = start(ranges(depthsGr)))
    depthsFiltered[, (names(mcols(depthsGr))) := data.frame(mcols(depthsGr)), with =
                     F]
    depthsFiltered
    
  }

#' Analyze CNV for a data frame of depths
#'
#' @param alldepths Data frame of depths
#' @param testsamps List of names of test samples, must match a column name in alldepths
#' @param normalsamps List of names of normal samples, must match a column name in alldepths
#' @param targetchr Chromosome of test gene
#' @param targetstart Start position of test gene
#' @param targetend End position of test gene
#' @param testRegions Data frame of regions for all test genes, used to exclude from CNV test
#' @param pagetitle Deprecated: Title of plot to be rendered
#' @param props Log2 levels to include in plot for reference
#'
#' @return
#' @export
#'
#' @examples
plot_one <- 
  function(alldepths,
           testsamps,
           normalsamps,
           targetchr,
           targetstart,
           targetend,
           testRegions,
           pagetitle = "10ng input",
           props = c(0, 0.01, 0.01, 0.03, 0.05, 0.1, 0.25, 1)) {
    # normed_mat should be a dataframe of the log2 normalized average coverage per capture region. It should also contain coordinates in the 'chr','start' and 'stop' columns
    
    testRegionsByBase = as.testRegionsByBase(testRegions)
    geneName = testRegions[chromosome == targetchr &
                             start == targetstart &
                             end == targetend]$gene[1]
    testRegionRemovedDepth = excludeByBase(alldepths, testRegionsByBase, geneName)
    
    allSampleNames <- c(testsamps, normalsamps)
    normed_mat <-
      data.frame(log2normalize(data.table(testRegionRemovedDepth), allSampleNames), check.names =
                   F)
    
    dedup_met_only <-
      colMeans(normed_mat[(normed_mat$Chromosome == targetchr) &
                            (normed_mat$Position >= targetstart) &
                            (normed_mat$Position <= targetend), testsamps, drop = FALSE])
    
    dedup_met_only_normals <-
      colMeans(normed_mat[(normed_mat$Chromosome == targetchr) &
                            (normed_mat$Position >= targetstart) &
                            (normed_mat$Position <= targetend), normalsamps,
                          drop = FALSE])

        metonly_mat <-
      t(outer(dedup_met_only[testsamps], dedup_met_only_normals, FUN = "-"))
    colnames(metonly_mat) <- colnames(metonly_mat)
    # Use normal versus normal as a comparison point
    normals <-
      t(outer(dedup_met_only_normals, dedup_met_only_normals, FUN = "-"))
    # There's no special mixture proportion in the normals so pool them all
    normals <- as.numeric(normals)
    bb <- as.list(data.frame(metonly_mat))
    if (!is.null(props)) {
      mediancurve <- apply(metonly_mat, 2, median)
      expected <-
        props * 2 ^ mediancurve[length(mediancurve)] + (1 - props) * 1
    }
    aaa <- sort(sapply(dedup_met_only[testsamps], function(x) {
      t.test(x,
             dedup_met_only_normals,
             alternative = 'greater',
             var.equal = TRUE)$p.value
    }))
    bbb <- sort(sapply(dedup_met_only[testsamps], function(x) {
      t.test(x,
             dedup_met_only_normals,
             alternative = 'less',
             var.equal = TRUE)$p.value
    }))

        list(
      right.sided = aaa,
      left.sided = bbb,
      diff.matrix = metonly_mat
    )
  }
