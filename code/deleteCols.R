## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

OLDEXT <- "_anz4.txt"
NEWEXT <- "_anz3.txt"

INPUT_HEADER <- TRUE
CONTEXTCOUNT <- TRUE
CLUSTERCOUNT <- FALSE
if (!CUSTOMMOTIFS) {
	source("findMotifs.R")
	source("countMotifs.R")
}
files <- dir(FILESDIR)

for (file in files) { 
	if (substr(file,nchar(file)-nchar(OLDEXT)+1,nchar(file))!=OLDEXT) {
		next
	}
	cat("\n", file, "\n")

	inputFile <- paste(FILESDIR, file, sep="")
	inputPref <- substr(inputFile, 1, nchar(inputFile)-nchar(OLDEXT))
	outputFile <- paste(inputPref, NEWEXT, sep="")
	output <- file(outputFile, "w")

	#if (COUNTMOTIFS) {
	#	newTitles <- paste(c(findTitles, contextTitle, countTitles), collapse="\t")
	#} else if (FINDMOTIFS) {
	#	newTitles <- paste(c(findTitles, contextTitle), collapse="\t")
	#} else {
	#	newTitles <- contextTitle
	#}

	input <- file(inputFile, "r", blocking = FALSE)
	if (INPUT_HEADER) {
		line <- readLines(input, n=1)
		#cat("H", substr(line, nchar(line), nchar(line)), "\n")
		fieldNames <- unlist(strsplit(line, "\t"))
	
		#removeTitles <- c(match(findTitles, fieldNames), match(countTitles, fieldNames))
		#duplicatedTitles <- countTitles[duplicated(countTitles)]
		#findTitles2remove <- setdiff(findTitles, c("G", "C", "tCw", "wGa"))
		# not necessary for tC_mutation	binary column, but added for clarity 10/28/14
		findTitles2remove <- setdiff(findTitles, c("G", "C", "tC", "Ga", "tCw", "wGa"))
		countTitles2remove <- setdiff(countTitles, c("g_counts", "c_counts", "tcw_counts", "wga_counts"))
		removeTitles <- (fieldNames %in% findTitles2remove) | (fieldNames %in% countTitles2remove)
		newFieldNames <- paste(fieldNames[!removeTitles], collapse="\t")
		writeLines(paste(newFieldNames,collapse="\t"), con=output)
	}

	c <- 0
	while (TRUE) {
	#while (c<10) {
		line <- readLines(input, n=1)
		lastChar <- substr(line, nchar(line), nchar(line))
		#if (lastChar=="\t") lastChar <- "TAB"
		#cat(c, lastChar, "\n")
		if (length(line)==0) {
			break
		}
		c <- c + 1
		#if (DEBUG && (c%%1000!=0)) {
		#	next
		#}
		fields <- unlist(strsplit(line, "\t"))
		# when anz4 input is defective (misses one field alread6y)
		#if (lastChar=="\t") fields <- c(fields, "", "")
		if (lastChar=="\t") fields <- c(fields, "")
		if (c%%1000==0) {
			cat(".")
		}
		newFields <- fields[!removeTitles]
		writeLines(paste(newFields,collapse="\t"), con=output)
	}

	close(input)
	close(output)

}
