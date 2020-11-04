## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

OLDEXT <- "_anz3.txt"
OLDEXT2 <- "_sum_all_fisher_Pcorr.txt"
NEWEXT <- "_anz5.txt"

if ((CUSTOMMOTIFS) & !(CUSTOMMUTLOAD)) {
	#newTitles <- c("[htCw_to_G+htCw_to_T]_per_mut", "htCw_to_G+htCw_to_T", "BH_Fisher_p-value_htCw", "APOBEChtCw_enrich", "htCw_to_G_enrich", "htCw_to_T_enrich", "p-value_GvT_skew", "BH_p-value_GvT_skew")
sampleTitles <- c(
paste("[", CUSTOMTYPE, "_to_G+", CUSTOMTYPE, "_to_T]_per_mut", sep=""),
paste(CUSTOMTYPE, "_to_G+", CUSTOMTYPE, "_to_T", sep=""),
paste("BH_Fisher_p-value_", CUSTOMTYPE, sep=""),
paste("APOBEC", CUSTOMTYPE, "_enrich", sep=""),
paste(CUSTOMTYPE, "_to_G_enrich",sep=""),
paste(CUSTOMTYPE, "_to_T_enrich",sep="")
#"p-value_GvT_skew",
#"BH_p-value_GvT_skew"
)
	#OLDEXT <- "_htCw_anz3.txt"
	#OLDEXT2 <- "_htCw_sum_all_fisher_Pcorr.txt"
	#NEWEXT <- "_htCw_anz5.txt"
	INSERTSAMPLEAFTERCOLNAME <- "Ga"
} else {
	if (! CUSTOMMOTIFS) {
	# revcomp added 05/15/18 
	sampleTitles <- c("[tCw_to_G+tCw_to_T+revcomp]_per_mut", "tCw_to_G+tCw_to_T+revcomp", "BH_Fisher_p-value_tCw", "APOBEC_enrich", "tCw_to_G_enrich", "tCw_to_T_enrich", "p-value_GvT_skew", "BH_p-value_GvT_skew", "APOBEC_MutLoad_MinEstimate")
	#sampleTitles <- c("[tCw_to_G+tCw_to_T]_per_mut", "tCw_to_G+tCw_to_T", "BH_Fisher_p-value_tCw", "APOBEC_enrich", "tCw_to_G_enrich", "tCw_to_T_enrich", "p-value_GvT_skew", "BH_p-value_GvT_skew", "APOBEC_MutLoad_MinEstimate")
	INSERTSAMPLEAFTERCOLNAME <- "APOBEC_mutation_to_T"
	}
}


sampleTitlesLength <- length(sampleTitles)

baseOrder <- c("A", "C", "G", "T")

INPUT_HEADER <- TRUE

files <- dir(FILESDIR)

for (file in files) { 
	if (substr(file,nchar(file)-nchar(OLDEXT)+1,nchar(file))!=OLDEXT) {
		next
	}
	sampleFileName <- paste(substr(file,1,nchar(file)-nchar(OLDEXT)), OLDEXT2, sep="")
	cat("\n", file, "\n")
	if (!(sampleFileName %in% files)) next
	sampleFile <- paste(FILESDIR, sampleFileName, sep="")
	
	table2add <- read.table(sampleFile, sep="\t", header=TRUE, as.is=TRUE)
	samples <- table2add[,"Sample"]
	sampleFileCon <- file(sampleFile, "r")
	sampleFileHeaders <- unlist(strsplit(readLines(sampleFileCon, n=1), "\t"))
	sampleFileHeaders <- gsub('"', "", sampleFileHeaders, perl=TRUE)
	
	cols2add <- match(sampleTitles, sampleFileHeaders)
	
	close(sampleFileCon)

	inputFile <- paste(FILESDIR, file, sep="")
	inputPref <- substr(inputFile, 1, nchar(inputFile)-nchar(OLDEXT))
	outputFile <- paste(inputPref, NEWEXT, sep="")
	output <- file(outputFile, "w")

	input <- file(inputFile, "r", blocking = FALSE)
	if (INPUT_HEADER) {
		line <- readLines(input, n=1)

		fieldNames <- unlist(strsplit(line, "\t"))
		INPUT_SAMPLE_FIELD <- which(fieldNames=="Tumor_Sample_Barcode")
		#INPUT_SAMPLE_FIELD <- which(fieldNames=="patient")
		INPUT_CHR_FIELD <- which(fieldNames=="Chromosome")
		INPUT_START_FIELD <- which(fieldNames=="Start_position")
		insertAfterColNum <- which(fieldNames=="\"CONTEXT(+/-20)\"")-1
		if (insertAfterColNum==length(fieldNames)) {
			newFieldNames <- c(fieldNames[1:insertAfterColNum], sampleTitles)
		} else {
			newFieldNames <- c(fieldNames[1:insertAfterColNum], sampleTitles, fieldNames[(insertAfterColNum+1):(length(fieldNames))])
		}
		writeLines(paste(newFieldNames,collapse="\t"), con=output)

		fieldNameCount <- length(fieldNames)
	}

	c <- 0
	mergePointer <- 1
	#while (FALSE) {
	while (TRUE) {
	#while (c<1000) {
		line <- readLines(input, n=1)
		lastChar <- substr(line, nchar(line), nchar(line))
		#if (lastChar=="\t") lastChar <- "TAB"
		#cat(c, lastChar, "\n")
		if (length(line)==0) {
			break
		}
		c <- c + 1
		if (c%%1000==0) {
			cat(".")
		}
		#if (DEBUG && (c%%1000!=0)) {
		#	next
		#}
		fields <- unlist(strsplit(line, "\t"))
		# fixed 07/14/15
		if (lastChar=="\t") fields <- c(fields, "")
		fieldCount <- length(fields)
		
		if (all(fields=="")) {
			break
		}

		inputSample <- fields[INPUT_SAMPLE_FIELD]
		sampleRow <- which(inputSample==samples)
		fields2add <- table2add[sampleRow, cols2add]
		
		if (insertAfterColNum==length(fieldNames)) {
			newFields <- c(fields[1:insertAfterColNum], fields2add)
		} else {
			newFields <- c(fields[1:insertAfterColNum], fields2add, fields[(insertAfterColNum+1):(length(fields))])
		}
		writeLines(paste(newFields,collapse="\t"), con=output)


	}
		
	close(input)
	close(output)

}
