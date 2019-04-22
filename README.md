# tumor-tissue-cnv-caller

The copy number analysis module identifies copy number variation for NGS data
with capture-based sequencing. This is the R-based implementation of the CNV caller.

## Table of Contents

* [tumor-tissue-cnv-caller](#tumor-tissue-cnv-caller)
* [Repository Organization](#repository-organization)
* [Usage](#usage)

## Repository Organization

```
└── src
    ├── CNV.calling.R - Top level script for CNV calling
    ├── ctdnaCopyNumberAnalysis
    │   ├── R - CNV analysis scripts and utilities
    │   ├── DESCRIPTION - R package description file
    │   ├── NAMESPACE - R package namespace file
```

## Usage

The following command generates three output files. A result file containing only CNV calls along with the genes' coordinates, name of gene, p-value, and indication of copy number gain or loss is emitted with the \*.cnv-call.bed suffix. A whitelist result file with the \*.cnv-whitelist.bed suffix is emitted and contains all genes considered for CNV. This output contains the gene chromosome, start, end, name, p-value, and the threshold used for making the call. Finally, an R RDS file with the suffix \*.cnv-internal.rds is emitted and contains intermediate data used for CNV calling and is kept for debugging or support purposes.

```
Usage: ./CNV.calling.R [options]

Options:
        -n NORMAL_BACKGROUND, --normal=NORMAL_BACKGROUND
                Normal background file, rds file

        -d TEST_DEPTH_FILE, --depth=TEST_DEPTH_FILE
                Input test sample depth file, could be 1 or multiple freq files, see --type

        -c 1/2, --type=1/2
                Type of test depth file, (1 - single freq file; 2 - txt file with multiple freq files, one file per line)

        -s SELECTOR_REGION_FILE, --selector=SELECTOR_REGION_FILE
                Input selector region file

        -t TEST_REGION_FILE, --test=TEST_REGION_FILE
                Input test region file

        -l T/F, --gcCorrection=T/F
                Logical, whether to do GC correction [default: F]

        -f GC_FILE, --gcfile=GC_FILE
                Input GC content file, [Required with -l T]

        -r T/F, --remove=T/F
                Logical, whether to remove test regions when doing GC correct, [Required when -l T]

        -p FRACTION, --fraction=FRACTION
                Fraction of most uniform positions to be used, 0~1 [default: 1]

        -o OUTPUT_DIR, --outputDir=OUTPUT_DIR
                Output directory [default: tmp]

        -h, --help
                Show this help message and exit

```
