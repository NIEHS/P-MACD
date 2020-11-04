## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

# Quick launch
# begin an R session inside of the package top directory "pmacd"
# run these commands by copy/paste into the R interpreter (uncommented) to launch - do not uncomment in this script
# change PARENTDIR to your package top directory:
# PARENTDIR = "/home/username/pmacd/"
# change REFFASTA to the path for your reference genome FASTA file:
# REFPATH = "/home/username/reference_genome/hg19.fa"
# source(paste0(PARENTDIR, "code/control.R"))
# main()
#
# to run a custom motif,
# set:
# CUSTOMMOTIFS <- TRUE
# change in main() function call below:
# filesDir = paste0(PARENTDIR, "work/res_[custom motif]/"),
# re-source the edited script:
# source(paste0(PARENTDIR, "code/control.R"))
# and run:
# main()
# end Quick launch

#CUSTOMMOTIFS <- TRUE
CUSTOMMOTIFS <- FALSE
CUSTOMTYPE <- "yCn"
#CUSTOMTYPE <- "nCg"
#CUSTOMTYPE <- "yTt"
#CUSTOMTYPE <- "rTt"
# do not change here, will be overridden in customConfig
CUSTOMMUTLOAD <- FALSE
USESAMPLEFILE <- FALSE
# will utilize file subset filtering (for sumContext)
BATCH <- FALSE


# execute this function to launch the pipeline
main <- function(
scriptDir = paste0(PARENTDIR, "code/"),
dataDir = paste0(PARENTDIR, "data/"),
# use this option for the default APOBEC pipeline
filesDir = paste0(PARENTDIR, "work/res_tCw/"),
# use this option for the pipeline using a custom motif
# to use a different motif, a folder must be created by the user, e.g. res_rTt
#filesDir = paste0(PARENTDIR, "work/res_yCn/"),
refFasta = REFPATH,
mutMultiplier = 1.0
)
{
SCRIPTDIR <<- scriptDir
FILESDIR <<- filesDir
DATADIR <<- dataDir
# Reference Genome
REF_FA <<- refFasta
# Fraction of Genome Sequenced
MUTMULTIPLIER <<- mutMultiplier

cat(DATADIR, "\n")
setwd(SCRIPTDIR)

if (CUSTOMMOTIFS) {
	customConfigScript <- paste("customConfig_", CUSTOMTYPE, ".R", sep="")
	source(customConfigScript)
}

time00 <- Sys.time()
time10 <- Sys.time()
source("ClusterFinder_standAlone.R")
time20 <- Sys.time()
source("setupPy.R")
source("addContextWrapper.R")
time30 <- Sys.time()
source("deleteCols.R")
time40 <- Sys.time()
source("sumContextWrapper.R")
time50 <- Sys.time()
source("analyzeAutoExt.R")
time60 <- Sys.time()
source("colSumAutoMulti.R")
time70 <- Sys.time()
source("addSampleCols.R")
time80 <- Sys.time()
time90 <- Sys.time()
time100 <- Sys.time()

cat("ClusterFinder: "); print(time20-time10)
cat("addContext: "); print(time30-time20)
cat("deleteCols: "); print(time40-time30)
cat("sumContext: "); print(time50-time40)
cat("analyze: "); print(time60-time50)
cat("colSum: "); print(time70-time60)
cat("addSampleCols: "); print(time80-time70)
}
