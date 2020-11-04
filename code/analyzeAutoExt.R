## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

#OLDEXT <- "_sum_all.txt"
OLDEXT <- ".txt"
NEWEXT <- "_fisher.txt"

INPUT_HEADER <- TRUE

DEBUG <- TRUE
DEBUG <- FALSE

files <- dir(FILESDIR)
#files <- files[1]

if (USESAMPLEFILE) {
if (PANCANCER) {
	sampleFile <- "panCanSample.txt"
	FIRST_NUM_COL <- 3
} else if (CUSTOMMOTIFS) {
	sampleFile <- paste("fisherSample_", CUSTOMTYPE, ".txt", sep="")
	#sampleFile <- "htcwSample.txt"
	FIRST_NUM_COL <- 2
} else {
	sampleFile <- "newcolsSample3.txt"
	FIRST_NUM_COL <- 2
}
sample <- file(sampleFile, "r", blocking = FALSE)
line <- readLines(sample, n=1)
newFieldNames <- unlist(strsplit(line, "\t"))
close(sample)
} else {
FIRST_NUM_COL <- 2
}

for (file in files) { 
	#if (substr(file,nchar(file)-nchar(OLDEXT)+1,nchar(file))!=OLDEXT) {
	if (regexpr("_sum", file)==-1) {
		next
	}

	inputPref <- substr(file, 1, nchar(file)-nchar(OLDEXT))

	# accommodate COMPLETESUBST 03/27/17
	inputFile <- paste(FILESDIR, file, sep="")
	#inputPref <- substr(inputFile, 1, nchar(inputFile)-nchar(OLDEXT))
	input <- file(inputFile, "r", blocking = FALSE)
	if (exists("COMPLETESUBST")) {
		subDir <- "tmp/"
		if ((regexpr("_fisher", file)!=-1) |
		 (substr(file, nchar(file)-nchar("_ClustSize.txt")+1, nchar(file))=="_ClustSize.txt") |
		 (substr(file, 1, 3)=="MAF") )
		{
			next
		}
		inputPref <- sub(PARENTTYPE, CUSTOMTYPE, inputPref)
	} else {
		subDir <- ""
	}
	cat("\n", file, "\n")

	outputFile <- paste(FILESDIR, subDir, inputPref, NEWEXT, sep="")
	output <- file(outputFile, "w")
	# write header w/o waiting for comment lines
	#writeLines(paste(newFieldNames,collapse="\t"), con=output)

	line <- readLines(input, n=1)
	while (substr(line,1,1)=="#") {
		# retain header comments
		writeLines(line, con=output)
		line <- readLines(input, n=1)	
	}

	fieldNames <- unlist(strsplit(line, "\t"))
if (USESAMPLEFILE) {
	writeLines(paste(newFieldNames,collapse="\t"), con=output)

	# positions of where the newFieldNames (as defined in the sample file)
	# are located in the input file
	# columns to be calculated are not present and will be NA
	fieldPosMatches <- match(newFieldNames, fieldNames)
	fieldPos2Get <- fieldPosMatches[!is.na(fieldPosMatches)]
	# positions of newFieldNames that exists in the input file
	oldFieldPos <- which(!is.na(fieldPosMatches))
	oldFieldCount <- length(oldFieldPos)
} else {

	if (!CUSTOMMOTIFS) {
		ANALYZECOLS <- c("tCw", "wGa")
		# MUTFIELDNAME should agree with custFields2
		# revcomp added 05/15/18 
		#MUTFIELDNAME <- "tCw_to_G+tCw_to_T"
		MUTFIELDNAME <- "tCw_to_G+tCw_to_T+revcomp"
		# "[nTt_to_C+aAn_to_G]_per_mut"
		perMutFieldName <- paste("[", MUTFIELDNAME, "]_per_mut", sep="")
		
		CUSTOMSKEW <- FALSE
		SKEWTYPE <- "GvT"
		skewPvalueName <- paste("p-value_", SKEWTYPE, "_skew", sep="")

		newFieldTitles <- c(
			perMutFieldName,
			MUTFIELDNAME,
			# revcomp added 05/15/18 
			#"[(C_to_G)+(C_to_T)]-[(tCw_to_G)+(tCw_to_T)]",
			"[(C_to_G)+(C_to_T)]-[(tCw_to_G)+(tCw_to_T)]+revcomp",
			#"tca+tga",
			#"c-tca",
			"tcw+revcomp",
			"c-tcw+revcomp",
			# revcomp added 05/15/18 END
			"Fisher_p-value_tCw",
			"odds_ratio_tCw",
			"95%_confidence_lower_tCw",
			"95%_confidence_upper_tCw",
			"APOBEC_enrich",
			"tCw_to_G_enrich",
			"Fisher_p-value_tCw_to_G",
			"odds_ratio_tCw_to_G",
			"tCw_to_T_enrich",
			"Fisher_p-value_tCw_to_T",
			"odds_ratio_tCw_to_T",
			skewPvalueName,
			"non-APOBEC_mutations",
			"non-APOBEC_substitutions",
			"APOBEC_to_G",
			"APOBEC_to_T"
		)

		sampleTitles <- c(
			"[tCw_to_G+tCw_to_T]_per_mut",
			"tCw_to_G+tCw_to_T", 
			"BH_Fisher_p-value_tCw", 
			"APOBEC_enrich", 
			"tCw_to_G_enrich",
			"tCw_to_T_enrich",
			"p-value_GvT_skew",
			"BH_p-value_GvT_skew"
		)

	}

	fieldPos2Get <- 1:32
	analyzeCol1 <- which(fieldNames==ANALYZECOLS[1])
	fieldPos2Get <- c(fieldPos2Get, (analyzeCol1-3):(analyzeCol1+6))
	analyzeCol2 <- which(fieldNames==ANALYZECOLS[2])
	fieldPos2Get <- c(fieldPos2Get, (analyzeCol2-3):(analyzeCol2+6))

	# new columns that are taken from the input
	oldFieldPos <- 1:(32+10+10)
	oldFieldCount <- length(oldFieldPos)
	newFieldNames <- c(fieldNames[fieldPos2Get],newFieldTitles)
	writeLines(paste(newFieldNames,collapse="\t"), con=output)
}
	
	count <- 0
	while (TRUE) {
	#while (c<200) {
		line <- readLines(input, n=1)
		if (length(line)==0) {
			break
		}
		count <- count + 1
		#if (DEBUG && (c%%1000!=0)) {
		#	next
		#}
		fields <- unlist(strsplit(line, "\t"))
		if (count%%1000==0) {
			cat(".")
		}
		oldFields <- fields[fieldPos2Get]
		eval(parse(text=paste(newFieldNames[oldFieldPos][FIRST_NUM_COL:oldFieldCount], "=", oldFields[FIRST_NUM_COL:oldFieldCount])))
		
		if (CUSTOMMOTIFS) {
			#customAnalysisScript <- paste("customAnalysis_", CUSTOMTYPE, ".R", sep="")
			#source(customAnalysisScript)
			#source("customAnalysis.R")
			outFields <- c(oldFields, getCustomFields())
		} else {
		#if (!CUSTOMMOTIFS) {

		## MODULES BEGIN HERE

		custFields <- c(
# [tCw_to_G+tCw_to_T]_per_mut
(tCw_to_G + tCw_to_T + wGa_to_C + wGa_to_A) / mutations,
# tCw_to_G+tCw_to_T
tCw_to_G + tCw_to_T + wGa_to_C + wGa_to_A,
# [(C_to_G)+(C_to_T)]-[(tCw_to_G)+(tCw_to_T)]
(C_to_G + C_to_T + G_to_C + G_to_A) - (tCw_to_G + tCw_to_T + wGa_to_C + wGa_to_A),
# tcw+wga
tcw + wga,
# c-tcw
(c + g) - (tcw + wga)
		)

		testMat <- rbind(c((tCw_to_G + tCw_to_T + wGa_to_C + wGa_to_A), (C_to_G + C_to_T + G_to_C + G_to_A - tCw_to_G - tCw_to_T - wGa_to_C - wGa_to_A)),c(tcw+wga, c+g-tcw-wga))
		ft <- fisher.test(testMat, alternative="greater")
		#ft <- fisher.test(testMat, alternative="less")
		
		testFields <- c(
# Fisher_p-value_tCw
ft$p.value,
# odds_ratio_tCw
ft$estimate,
#paste(ft$conf.int, collapse="-")
# 95%_confidence_lower_tCw
ft$conf.int[1],
# 95%_confidence_upper_tCw
ft$conf.int[2]
		)

		
		custFields2 <- c(
# APOBEC_enrich
((tCw_to_G + tCw_to_T + wGa_to_C + wGa_to_A) * (c + g))/((C_to_G + C_to_T + G_to_C + G_to_A) * (tcw + wga)))


		custFields3 <- c(
# tCw_to_G_enrich
((tCw_to_G + wGa_to_C) * (c + g))/((C_to_G + G_to_C) * (tcw + wga)))

		testMat3 <- rbind(c((tCw_to_G + wGa_to_C), (C_to_G + G_to_C - tCw_to_G - wGa_to_C)),c(tcw+wga, c+g-tcw-wga))
		ft <- fisher.test(testMat3, alternative="greater")
       
		testFields3 <- c(
# Fisher_p-value_tCw_to_G
ft$p.value,
# odds_ratio_tCw_to_G
ft$estimate
		)


		custFields4 <- c(
# tCw_to_T_enrich
((tCw_to_T + wGa_to_A) * (c + g))/((C_to_T + G_to_A) * (tcw + wga))
		)
		
		testMat4 <- rbind(c((tCw_to_T + wGa_to_A), (C_to_T + G_to_A - tCw_to_T - wGa_to_A)),c(tcw+wga, c+g-tcw-wga))
		ft <- fisher.test(testMat4, alternative="greater")
       
		testFields4 <- c(
# Fisher_p-value_tCw_to_T
ft$p.value,
# odds_ratio_tCw_to_T
ft$estimate
		)

		
		testMat5 <- rbind(c((tCw_to_G + wGa_to_C), (tCw + wGa - tCw_to_G - wGa_to_C)), c((tCw_to_T + wGa_to_A), (tCw + wGa - tCw_to_T - wGa_to_A)))
		ft <- fisher.test(testMat5, alternative="two.sided")
		
		testFields5 <- c(
# p-value_GvT_skew
ft$p.value
		)

		
		custFields6 <- c(
# non-APOBEC_mutations
mutations - tCw_to_G - tCw_to_T - wGa_to_C - wGa_to_A,
# non-APOBEC_substitutions
substitutions - tCw_to_G - tCw_to_T - wGa_to_C - wGa_to_A)

		custFields7 <- c(
# APOBEC_to_G
tCw_to_G + wGa_to_C,
# APOBEC_to_T
tCw_to_T + wGa_to_A)

		outFields <- c(oldFields,custFields,testFields,custFields2,custFields3,testFields3,custFields4,testFields4,testFields5,custFields6,custFields7)
		## MODULES END HERE
		}

		writeLines(paste(outFields,collapse="\t"), con=output)

	}

	close(input)
	close(output)

}

