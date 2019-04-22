#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(matrixStats))
suppressMessages(library(ctdnaCopyNumberAnalysis))

option_list <- list(
 	make_option(c("-i", "--samplename"), metavar = "sample_name", help="Sample identifier"),
 	make_option(c("-n", "--normal"), metavar = "normal_background", help="Normal background file, rds file"),
	make_option(c("-d", "--depth"), metavar = "test_depth_file", help="Input test sample depth file, could be 1 or multiple freq files, see --type"),
	make_option(c("-c", "--type"), metavar = "1/2", help="Type of test depth file, (1 - single freq file; 2 - txt file with multiple freq files, one file per line)"),
	make_option(c("-s", "--selector"), metavar = "selector_region_file", help = "Input selector region file"),
	make_option(c("-t", "--test"), metavar = "test_region_file", help = "Input test region file"),

	make_option(c("-l", "--gcCorrection"), metavar = "T/F",  default = "F", help = "Logical, whether to do GC correction [default: %default]"),
	make_option(c("-f", "--gcfile"), metavar = "GC_file",  help = "Input GC content file, [Required with -l T]"),
	make_option(c("-r", "--remove"), metavar="T/F", help = "Logical, whether to remove test regions when doing GC correct, [Required when -l T]"),	
	
	make_option(c("-p", "--fraction"), default = "1", help = "Fraction of most uniform positions to be used, 0~1 [default: %default]", type="numeric"),
	make_option(c("-u", "--cpus"), default = "2", help = "Number of threads to use [default: %default]", type="numeric"),
	
	make_option(c("-o", "--outputDir"), metavar="output_dir", default = "tmp", help = "Output directory [default: %default]")	
)

opt <- parse_args(OptionParser(option_list=option_list))


checkParam(opt$normal,"n (normal_background)") # in util.R
checkParam(opt$depth,"d (test_depth_file)") # in util.R
checkParam(opt$type,"c (test_depth_file_type)")
checkParam(opt$selector,"s (selector_region_file)")
checkParam(opt$test,"t (test_region_file)")

# check test_depth_file_type
if (opt$type != 1 & opt$type !=2){
        cat(paste("c (type)"," parameter must be 1/2. See script usage (-h)\n",sep=""))
        stop_quietly()
}

# check gcCorrection argument
if (opt$gcCorrection != "T" & opt$gcCorrection != "F"){
	cat(paste("l (gcCorrection)"," parameter must be logical. See script usage (-h)\n",sep=""))
        stop_quietly() # in util.R
}

# if gcCorrection is position, check parameters
if (opt$gcCorrection == "T"){
	checkParam(opt$gcfile,"f (GC_file)")
	checkParam(opt$remove,"r (remove)")	
	if (opt$remove != "T" & opt$remove != "F"){
		cat(paste("r (remove)"," parameter must be logical. See script usage (-h)\n",sep=""))
        	stop_quietly() # in util.R
	}
	checkFile(opt$gcfile)
	if (opt$remove != "T" & file.info(opt$test)$size == 0){
		cat("test_region_file cannot be empty when remove = T\n")
        	stop_quietly()
	}
	opt$remove = as.logical(opt$remove)
	opt$gcCorrection = as.logical(opt$gcCorrection)
}

# check fraction argument
if (opt$fraction <= 0 | opt$fraction >1){
        cat("fraction parameter must be between 0~1. See script usage (-h)\n")
        stop_quietly()
}



checkFile(opt$normal)
checkFile(opt$depth)
checkFile(opt$selector)
checkFile(opt$test)



## Main 
normalDepths = readRDS(opt$normal)
normalDepths = normalDepths[complete.cases(normalDepths),]

# Run CNV caller
if (opt$type == 1){
	if (is.null(opt$samplename)) {
		samplename = sapply(strsplit(basename(opt$depth), ".", fixed=TRUE), '[', 1)
	} else {
		samplename = opt$samplename
	}
	analyzeCNV(opt$depth, opt$selector, opt$test, normalDepths, samplename, opt$outputDir, opt$cpus, opt$fraction, opt$gcCorrection, opt$gcfile, opt$remove) 

}else if (opt$type == 2){
	print (opt$fraction)
	print (opt$gcCorrection)
	print (opt$gcfile)
	print(opt$remove)

	files = read.table(opt$depth,stringsAsFactors = F)
	for (i in 1:nrow(files)){
		samplename = files[i,1]
		analyzeCNV(files[i,2], opt$selector, opt$test, normalDepths, samplename, opt$outputDir, opt$cpus, opt$fraction, opt$gcCorrection, opt$gcfile, opt$remove) 
	}
}
