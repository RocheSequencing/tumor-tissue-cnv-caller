
#' Check existence of input parameter
#'
#' @param x (parameter to check)
#' @param info (message to display)
#' @export
checkParam <- function(x, info) {
        if (is.null(x)){
                cat(paste(info," parameter must be provided. See script usage (-h)\n",sep=""))
                stop_quietly()
        }
}


#' Exit program without error message

stop_quietly <- function() {
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
}



#' Calculate coefficient of variant per position
#'
#' @import GenomicRanges
#' @import data.table 
#' @import CopyNumberAnalysis
#' @param depths (depths per position per sample)
#' @param selector (selector region)
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
#' @param alldepths (depths per position per sample, can be single sample)
#' @param selector (selector regions)
#' @param testRegions (regions of genes to be tested)
#' @param gcfile (GC content of each selector region)
#' @param remove (whether exclude test genes when fitting GC curves, default is T) 
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



#' Call CNVs (for leave-one-out cross-validation) 
#'
#' @import data.table 
#' @import CopyNumberAnalysis
#' @param readDepth (depths from the test sample)
#' @param targetFile (selector regions file)
#' @param testRegionFile (regions of genes to be tested)
#' @param normalDepths (depths of background samples)
#' @return res (p.values of the test sample)
#'
#' @examples
#' analyzeCNV.cv(readDepth, targetFile, testRegionFile, normalDepths)

analyzeCNV.cv <- function(readDepth, targetFile, testRegionFile, normalDepths, numCores = 1){
        coordCols = c('Chromosome', 'Position')
        selector <- loadRSUSelector( targetFile )
        testRegions = read.rsu.tsv( testRegionFile)

        testSampleName <- setdiff(colnames(readDepth), coordCols)
        normSampleNames <- setdiff(colnames(normalDepths), coordCols)
        allSampleNames <- c(normSampleNames, testSampleName)

        alldepths <- merge( normalDepths, readDepth, by=coordCols)  # get depths of same coordinates
        alldepths <- intersectTarget(alldepths, selector)   # get positions within the selector regions
        alldepth_nz <- data.frame(alldepths, check.names = F)


        statistics <- dlply(testRegions, ~ gene + chromosome + start + end + pval.threshold,
        function(x) {
                plot_one (alldepth_nz, testSampleName, normSampleNames, x$chromosome,
                    x$start, x$end, testRegions = testRegions, pagetitle = "10ng input", props = NULL )
        }, .parallel = F)

        report <- ldply( statistics,
        function(x) {
            data.frame(one.side.greater = x$right.sided,
                      one.side.less=x$left.sided)
        }, .parallel = F)

	res = c(report$one.side.greater, report$one.side.less)
	names(res) = c(unlist(lapply(report$gene, function(x) paste(x,".gain",sep =""))), unlist(lapply(report$gene, function(x) paste(x,".loss",sep =""))))
	return(res)
}


#' calculate geometric mean of p.values (for leave-one-out cross-validation) 
#'
#' @param x (p.values)
#' @param num (geometric mean of 'num' smallest p.values, default = 2)
#' @return the geometric mean
#'
#' @examples
#' calculate.geometric.mean(x,2)

calculate.geometric.mean <- function(x,num=2){
	x = x[order(x,decreasing=F)]
	return(exp(mean(log(x[1:num]))))
}


#' calculate p.value threshold from set of normals 
#'
#' @param SampleNames (names of normal samples)
#' @param normalDepths (depths of normal samples)
#' @param targetFile (selector regions file)
#' @param num 
#' @param testRegionFile (regions of genes to be tested)
#'
#' @return p.value thresholds for test genes
#'
#' @examples
#' calculate.p.value.threshold.from.normal(SampleNames, normalDepths, targetFile, testRegionFile)

calculate.p.value.threshold.from.normal <- function(SampleNames, normalDepths, targetFile, testRegionFile, num){
	cv.res = NULL
	for (sample in SampleNames){
        	readDepth = normalDepths[,c("Chromosome","Position",sample),with=F]
        	bgdepths = normalDepths[,setdiff(colnames(normalDepths),sample),with=F]
		res = analyzeCNV.cv(readDepth, targetFile, testRegionFile, bgdepths, numCores = 2)
		cv.res = rbind(cv.res,res)
	}
	rownames(cv.res) = SampleNames
	threshold = apply(cv.res,2, function(x) calculate.geometric.mean(x,num))	
	return(threshold)
}



#' wrapper for cnvPanelizer
#' 
#' @import CNVPanelizer 
#' @param testRegionFile (test genes)
#' @param targetFile (selector regions)
#' @param referenceFile (a file contains the background samples, one sample per line, BAM file or FREQ file)
#' @param sampledata (test samples, BAM file or FREQ file)
#' @param sampleNames (test sample names)
#' @param detailed output in Rdata format
#' @param type ('freq' or 'bam')
#' @return res (CNV calling of the test genes)
#'
#' @examples
#' cnvPanelizer(testRegionFile, targetFile, referenceFile, sampledata, samplenames, outputRdata)
 
cnvPanelizer <- function(testRegionFile, targetFile, referenceFile, sampledata, samplenames, outputRdata, type = "freq"){
        
	suppressMessages(library(CNVPanelizer))
        genomicRangesFromBed <- BedToGenomicRanges(targetFile, ampliconColumn = 4, split = "_")
        metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
        geneNames = metadataFromGenomicRanges["geneNames"][, 1]
        ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

        referenceFilenames = read.table(referenceFile,stringsAsFactor=F)[,1]

	# if input are freq files, change the directory to the bam files
	if (type == "freq"){
        	referenceFilenames = gsub("(^.*)/(.*?)/(.*?).dualindex-deduped.sorted.bam.snv.freq", "\\1/\\2/bams/\\3.dualindex-deduped.sorted.bam", referenceFilenames)
        	sampleFilenames = gsub("(^.*)/(.*?)/(.*?).dualindex-deduped.sorted.bam.snv.freq", "\\1/\\2/bams/\\3.dualindex-deduped.sorted.bam", sampledata)
	}else{
		sampleFilenames = sampledata
	}

        # Counting reads from BAMs
        removePcrDuplicates <- FALSE
        referenceReadCounts <- ReadCountsFromBam(referenceFilenames, genomicRangesFromBed,
                                         sampleNames = referenceFilenames,
                                         ampliconNames = ampliconNames,
                                         removeDup = removePcrDuplicates)

        sampleReadCounts <- ReadCountsFromBam(sampleFilenames, genomicRangesFromBed,
                                      sampleNames = sampleFilenames,
                                      ampliconNames = ampliconNames,
                                      removeDup = removePcrDuplicates)

	# Normalization
        normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts, referenceReadCounts,
                                                 ampliconNames = ampliconNames)
        samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
        referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

        # Bootstrap
        replicates <- 10000

        bootList <- BootList(geneNames, samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)

        # Background Estimation
        backgroundNoise <- Background(geneNames, samplesNormalizedReadCounts,
                              referenceNormalizedReadCounts,
                              bootList,
                              replicates = replicates,
                              significanceLevel = 0.1,
                              robust = TRUE)

        # Report
        reportTables <- ReportTables(geneNames, samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             backgroundNoise)

        save(reportTables, file = outputRdata)

        # Collect result
	testRegion = read.table(testRegionFile, header=T, stringsAsFactors=F)
        res = do.call(rbind, lapply(1:length(reportTables), function(x) cbind(gene = testRegion$gene, sample = samplenames[x], reportTables[[x]][testRegion$gene, c("MeanRatio","Passed")])))

        return(res)

   	# Plot
        # PlotBootstrapDistributions(bootList, reportTables)
}




########  The following sections are used for generating plots 


#' calculate median depth per selector region (used for GC - depth plot)
#' 
#' @import CopyNumberAnalysis
#' @import GenomicRanges
#' @import data.table
#' @param alldepths (input depth)
#' @param selector (selector region)
#' @param testRegions (regions for test genes)
#' @param normalize (whether to do log2 normalization of depths, default: F)
#' @return normed.exons (median depth per selector region)
#' @return hits (overlapping with test genes)
#' @return samMed (median depth per sample after excluding test genes)
#'
#' @examples
#' cal.median.depth.per.region(depths, selector, testRegions)

cal.median.depth.per.region <- function(alldepths, selector, testRegions){
        
	coordCols = c('Chromosome', 'Position','region')
        SampleNames <- setdiff(colnames(alldepths), coordCols)

        alldepths <- intersectTarget(alldepths, selector)
        
	depthsGr = GenomicRanges::GRanges(seqnames = Rle(alldepths$Chromosome), ranges = IRanges::IRanges(start = alldepths$Position, end = alldepths$Position))
        targetGr = GenomicRanges::GRanges(seqnames = Rle(selector$Chromosome), ranges = IRanges::IRanges(start = selector$Start, end = selector$End))

        hits <- findOverlaps(depthsGr, targetGr)
        mcols(hits) = hits
        
	alldepths$region = data.frame(mcols(hits))[,2]

        normed_mat <- data.frame(alldepths)
        
	normed_mat$region = alldepths$region

        # Median depth per selector region
	res = apply(data.frame(normed_mat[,3:(ncol(normed_mat)-1)]), 2, function(x) aggregate(val ~ region, cbind(val = x, region = normed_mat$region), median))
	df <- data.frame(matrix(unlist(res), ncol=length(res), byrow=F),stringsAsFactors=FALSE)
        normed.exons <- data.frame(df[(length(unique(alldepths$region))+1):nrow(df),])
        rownames(normed.exons) = res[[1]]$region

        # Get selector regions that are test regions
        depthsGr = GenomicRanges::GRanges(seqnames = Rle(testRegions$chromosome), ranges = IRanges::IRanges(start = testRegions$start, end = testRegions$end))
        hits <- findOverlaps(depthsGr, targetGr)
        mcols(hits) = hits
        hits = data.frame(mcols(hits))

        samMed = apply(data.frame(normed_mat[normed_mat$region %in% hits[,2]==FALSE,3:(ncol(normed_mat)-1)]), 2, median)
        return(list(normed.exons = normed.exons, hits=hits, samMed=samMed))
}


#' Generate GC - depth plot
#' 
#' @param normed.exon (vector of median depths per selector region of the sample)
#' @param hits (from 'cal.median.depth.per.region', record the selector regions that are overlapping with testRegions)
#' @param region.names (the index of selector regions for normed.exon)
#' @param title (title of the plot, usually the sample name)
#' @param fit (fit curve from GC correction)
#' @param testRegions (test genes)
#' @param gcfile (GC content file)
#' @return generate the plot
#'
#' @examples
#' plot.gc.depth(normed.exon, hits, region.names, title = "",fit, testRegions, gcfile)

plot.gc.depth <- function(normed.exon, hits, region.names, title = "",fit, testRegions, gcfile){
	gc <- read.table(gcfile,row.names=1)	

	names(normed.exon) = region.names
	mat = cbind(gc[region.names,3], normed.exon)
	 
	plot(gc[region.names,3], normed.exon,xlab="GC",ylab="Coverage", main = title, col= "#c6c5c4")
	lines(gc[region.names,3][order(gc[region.names,3])],predict(fit,gc[region.names,3])[order(gc[region.names,3])], lty=2, type="l")

	num = sort(unique(hits[,1]))
	colors = rainbow(length(num))
	for (i in 1:length(num)){
		points(gc[as.character(hits[hits[,1] == num[i],2]),3], normed.exon[as.character(hits[hits[,1] == num[i],2])],col= colors[i], pch=20)	
	}
	legend('topright',bty='n', pch=20, col = colors, testRegions[num,]$gene,cex = 0.9)
}



#' Generate Mappability - depth plot
#' 
#' @param normed.exon (vector of median depths per selector region of the sample)
#' @param hits (from 'cal.median.depth.per.region', record the selector regions that are overlapping with testRegions)
#' @param region.names (the index of selector regions for normed.exon)
#' @param title (title of the plot, usually the sample name)
#' @param testRegions (test genes)
#' @param map (mapability file)
#' @return generate the plot
#'
#' @examples
#' plot.mappability.depth(normed.exon, hits, region.names, title = "",testRegions, map)

plot.mappability.depth <- function(normed.exon, hits, region.names, title = "",testRegions, map){
        map <- read.table(map,row.names=1)

        names(normed.exon) = region.names
        mat = cbind(map[region.names,3], normed.exon)

        plot(map[region.names,3], normed.exon,xlab="Mapability",ylab="Coverage", main = title, col= "#c6c5c4")

        num = sort(unique(hits[,1]))
        colors = rainbow(length(num))
        for (i in 1:length(num)){
                points(map[as.character(hits[hits[,1] == num[i],2]),3], normed.exon[as.character(hits[hits[,1] == num[i],2])],col= colors[i], pch=20)
        }
        legend('topright',bty='n', pch=20, col = colors, testRegions[num,]$gene,cex = 0.9)
}



#' Generate BAF plot
#' 
#' @param file (freq)
#' @param selector (selector regions)
#' @param sampleName (sampleName) 
#' @param testRegions (test genes)
#' @return generate the plot
#'
#' @examples
#' plot.BAF(file, selector, sampleName, testRegions)

plot.BAF <- function(file, selector, sampleName, testRegions){
        data = fread(file, header=T)
        colnames(data)[1:2] = c('Chromosome', 'Position')
        alldepths <- intersectTarget(data, selector)
  	
	chr = c("chr1", "chr2","chr3", "chr4","chr5","chr6", "chr7", "chr8","chr9", "chr10","chr11","chr12", "chr13", "chr14","chr15", "chr16","chr17","chr18","chr19", "chr20","chr21", "chr22","chrX")
	alldepths$Chromosome = factor(alldepths$Chromosome, levels=chr)
	alldepths = alldepths[order(alldepths$Chromosome),]
 
	alldepths$mut  = apply(alldepths[,7:14,with=F],1,sum)
        alldepths$BAF = alldepths$mut/alldepths$DEPTH

        depthsGr = GenomicRanges::GRanges(seqnames = Rle(alldepths$Chromosome), ranges = IRanges::IRanges(start = alldepths$Position, end = alldepths$Position))
	testGr = GenomicRanges::GRanges(seqnames = Rle(testRegions$chromosome), ranges = IRanges::IRanges(start = testRegions$start, end = testRegions$end))

        hits <- findOverlaps(depthsGr, testGr)
        mcols(hits) = hits
        hits = data.frame(mcols(hits))

        plot(1:nrow(alldepths), alldepths$BAF,xlab="Position",ylab="BAF",pch=20,main=sampleName)
        wid = 0.5
        abline(h=1,lty=2,col="red",lwd=wid)
        abline(h=0,lty=2,col="red",lwd=wid)
        abline(h=0.5,lty=2,col="red",lwd=wid)
        abline(h=1/3,lty=2,col="red",lwd=wid)
        abline(h=2/3,lty=2,col="red",lwd=wid)

	num = sort(unique(hits[,2]))
        colors = rainbow(length(num))

        for (i in 1:length(num)){
		gene = testRegions[num[i],]$gene
		points(hits[hits[,2] == num[i],1] ,alldepths[hits[hits[,2] == num[i],1],]$BAF, col=colors[i],pch=20)
        }
        
	legend(x= 0,y=0.9, bty='n', pch=20, col = colors, testRegions[num,]$gene,cex = 0.9)

        #sel = alldepths[hits[,1],]
        #sel$name  = name
        #return(sel[sel$BAF > 0.1 & sel$BAF < 0.9,c(1:14,16,18,19),with=F])
}



#' Calculate log2 normalization of depths 
#' 
#' @param alldepths (depths)
#' @param selector (selector regions)
#' @return log2 transformed depths
#'
#' @examples
#' coverage_matrix(alldepths,selector)

coverage_matrix = function(alldepths,selector ){
        alldepths <- intersectTarget(alldepths, selector)  #### selector position
        fs_dt = as.data.frame(alldepths)
        medians = apply(fs_dt[,3:ncol(fs_dt)],2,median,na.rm=T)
        fs_dt[,3:ncol(fs_dt)] = sweep(fs_dt[,3:ncol(fs_dt)],2,medians,FUN='/')
        fs_dt_log = fs_dt
        fs_dt_log[,3:ncol(fs_dt_log)] =  log2(fs_dt_log[,3:ncol(fs_dt_log)])

        fs_log_mat = as.matrix(fs_dt_log[,3:ncol(fs_dt_log)])
        #fs_log_mat[is.na(fs_log_mat)] = -10
        rMns = rowMeans(fs_log_mat,na.rm=T)
        rSds = rowSds(fs_log_mat,na.rm=T)
        rQuants = rowQuantiles(fs_log_mat,probs=c(.05,.5,.95),na.rm=T)
        # Merge everything back in
        fs_log_merged = cbind(fs_dt_log[,c(1,2)],rQuants,avg=rMns,sds=rSds)
        fs_log_merged$ind = seq(1,nrow(fs_log_merged))
        fs_log_merged = cbind(fs_log_merged,fs_log_mat)
	return(fs_log_merged)
}


#' Plot coverage for each gene
#'
#' @param normalDepths (normal depths)
#' @param testDepths (test depths)
#' @param selector (selector regions)
#' @param testRegions (test genes)
#' @return coverage plot
#'
#' @examples
#' coverage.plot(normalDepths, testDepths, selector, testRegions)

coverage.plot <- function(normalDepths, testDepths, selector, testRegions){
	
	coordCols = c('Chromosome', 'Position')
	normSampleNames <- setdiff(colnames(normalDepths), coordCols)
	testSampleNames <- setdiff(colnames(testDepths), coordCols)
	alldepths <- merge( normalDepths, testDepths, by=coordCols) 

	cc = coverage_matrix(alldepths, selector)
	
	depthsGr = GenomicRanges::GRanges(seqnames = Rle(cc$Chromosome), ranges = IRanges::IRanges(start = cc$Position, end = cc$Position))
       	testGr = GenomicRanges::GRanges(seqnames = Rle(testRegions$chromosome), ranges = IRanges::IRanges(start = testRegions$start, end = testRegions$end))

        hits <- findOverlaps(depthsGr, testGr)
        mcols(hits) = hits
	hits = mcols(hits)

	num = sort(unique(hits[,2]))
	par(mfrow=c(length(num),1), mar=c(4,4,4,4))
	for (i in num)
	{
		gene = testRegions[i,]$gene
        	ss = cc[hits[hits[,2] ==i,1],]
        	names = names(ss)[9:(ncol(alldepths)+6)]
        	cols = ifelse((names %in% testSampleNames),'#FF000022','#00FF0022')
        	matplot(ss$ind,ss[,names],col=cols,pch='.',main=gene ,xlab="Position",ylab="log2(depth/median depth)", ylim=c(-4,4))
		abline(h=0, lty=2,lwd=0.5)
	}
}



#' Analyze Sensitivity and Specificity
#'
#' @param path to results
#' @param pos (positive samples, separated by ','
#' @param neg (negative samples, separated by ','
#' @param test.gene (Gene tested, one gene)
#' @param type (gain or loss)
#' @param p (p value threshold)
#' @return sensitivity & specificity
#'
#' @examples
#' analyze.result(path, pos, neg, test.gene, type, p) 

analyze.result <- function(path, pos, neg, test.gene, type, p){
	p = as.numeric(p)
	filenames <- list.files(path, pattern="*.cnv-whitelist.bed", full.names=TRUE)
	ldf <- lapply(filenames, function(x) read.table(x, stringsAsFactors=F))
	tmp = lapply(filenames, function(x) unlist(strsplit(x,"/")))
	
	
	names = do.call(rbind, lapply(tmp,function(x) get.name(x))) 
	genes = ldf[[1]]$V4

	values = cbind(names, do.call(rbind,lapply(ldf, function(x) c(x$V5, x$V6))))

	pos.data = values[values[,1] %in% pos,]
	neg.data = values[values[,1] %in% neg,]
	
	mass.pos = unique(pos.data[,2])
	sensitivity= matrix(NA, nrow=length(mass.pos), ncol = length(pos))
	rownames(sensitivity)= mass.pos
	colnames(sensitivity) = pos
	
	mass.neg = unique(neg.data[,2])
	specificity= matrix(NA, nrow=length(mass.neg), ncol = length(neg))
	rownames(specificity)= mass.neg
	colnames(specificity) = neg

	sen.spc <- function(data,masses,samples,res, test,type, tag){
		if (type == "gain"){
			colm = which(genes==test)
		}else{
			colm = which(genes==test) + length(genes)
		}
		
		for (mass in as.character(masses)){
			for (sample in samples){
				sel = as.matrix(data[data[,1]==sample & data[,2]==mass,])
				if (ncol(sel) == 1){
					sel = t(sel)
				}
				values = as.numeric(sel[,colm+3])
				if (tag == "sen"){
					res[mass,sample] = paste(sum(values < p),"/",length(values),sep="")					     }else{
					res[mass,sample] = paste(sum(values > p),"/",length(values),sep="")
				}
			}
		}
		return(res)
	}
	
	sensitivity = sen.spc(pos.data, mass.pos, pos, sensitivity, test.gene,type, "sen") 
	specificity = sen.spc(neg.data, mass.neg, neg, specificity, test.gene,type, "spe") 
	return(list(sensitivity,specificity))
}


#' Get sample and mass info from file name
get.name <- function(x){
	ele = x[length(x)]
        full = gsub(".cnv-whitelist.bed","",ele)
        res = c(unlist(strsplit(full,"_"))[1:2],full)
        res[2] = unlist(strsplit(res[2],"-"))[1]
        return(res)
}

#' Put results together
#'
#' @param path to results
#' @return compile results
#'
#' @examples
#' compile.results(path) 

compile.results <- function(path){
        filenames <- list.files(path, pattern="*.cnv-whitelist.bed", full.names=TRUE)
        ldf <- lapply(filenames, function(x) read.table(x, stringsAsFactors=F))
        tmp = lapply(filenames, function(x) unlist(strsplit(x,"/")))

        names = do.call(rbind, lapply(tmp,function(x) get.name(x)))
        genes = ldf[[1]]$V4

        values = cbind(names, do.call(rbind,lapply(ldf, function(x) c(x$V5, x$V6))))
	colnames(values) = c("sample","tag","name", unlist(lapply(genes, function(x) paste(x,".gain",sep =""))), unlist(lapply(genes, function(x) paste(x,".loss",sep =""))))
	return(values)
}



#' Calculate sensitivity and specificity for each p value threshold
roc.cal <- function(A){
        thresh = sort(unique(A[,1]))
        sen = c()
        spec = c()
        res=NULL
        for (i in 1:length(thresh)){
                sen = c(sen, sum(A[A[,1] < thresh[i],2])/sum(A[,2]))
                spec = c(spec, sum(A[A[,1] > thresh[i],2]==0)/sum(A[,2]==0))
        }
        res$sen = sen
        res$spec = spec
        return(res)
}



#' Generate ROC curve based on positive and negative samples 
#' 
#' @param Path to results
#' @param pos, positive samples, separated by ','
#' @param neg, negative samples, separated by ','
#' @param p, p value threshold, will shown as a dot on the ROC curve
#' @param gene, gene to evaluate, single gene
#' @param type, must be gain or loss
#' @return a plot
#'
#' @examples
#' plot.auc(path, pos, neg, p, gene, type)

plot.auc <- function(path, pos, neg, p, test.gene, type){
	filenames <- list.files(path, pattern="*.cnv-whitelist.bed", full.names=TRUE)
        ldf <- lapply(filenames, function(x) read.table(x, stringsAsFactors=F))
        tmp = lapply(filenames, function(x) unlist(strsplit(x,"/")))

        names = do.call(rbind, lapply(tmp,function(x) get.name(x)))
        genes = ldf[[1]]$V4

        values = cbind(names, do.call(rbind,lapply(ldf, function(x) c(x$V5, x$V6))))
	if (type == "gain"){
        	colm = which(genes==test.gene)
        }else{
        	colm = which(genes==test.gene) + length(genes)
        }
	pos.data = as.numeric(values[values[,1] %in% pos,colm+3])
	neg.data = as.numeric(values[values[,1] %in% neg, colm+3])
	res = cbind(c(pos.data, neg.data), c(rep(1, length(pos.data)),rep(0,length(neg.data))))

	roc = roc.cal(res)
	plot(1-roc$spec,roc$sen,lty=2,type="l",xlab="False Positive Rate",ylab="True Positive Rate", main = paste(test.gene,type, sep=":"))
 	
	if (!is.null(p)){
		p=as.numeric(p)
		sens = sum(res[res[,1] < p,2])/sum(res[,2])
        	spec = sum(res[res[,1] > p,2]==0)/sum(res[,2] ==0 )
       	 	points(1-spec,sens,pch=15)
	}
	simple_auc <- function(TPR, FPR){
  		# inputs already sorted, best scores first 
  		dFPR <- c(diff(FPR), 0)
  		dTPR <- c(diff(TPR), 0)
  		sum(TPR * dFPR) + sum(dTPR * dFPR)/2
	}
	return(simple_auc(roc$sen,1-roc$spec))
}


