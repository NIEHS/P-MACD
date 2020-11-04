## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

#OLDEXT <- "_sum_all_fisher.txt"
OLDEXT <- "_fisher.txt"
NEWEXT <- "_fisher_Pcorr.txt"

INPUT_HEADER <- TRUE

DEBUG <- TRUE
DEBUG <- FALSE

getChisqPvalue <- function(vec) {
		return(chisq.test(vec)$p.value)
}

files <- dir(FILESDIR)
#files <- files[1]

for (file in files) { 
	if (substr(file,nchar(file)-nchar(OLDEXT)+1,nchar(file))!=OLDEXT) {
	#if (substr(file,nchar(file)-8,nchar(file))!="_anz1.txt") {
	#if (substr(file,1,5)!="ME001") {
		next
	}
	cat("\n", file, "\n")

	inputFile <- paste(FILESDIR, file, sep="")
	inputPref <- substr(inputFile, 1, nchar(inputFile)-nchar(OLDEXT))
	inputTable <- read.table(inputFile, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment="#")
	outputFile <- paste(inputPref, NEWEXT, sep="")
	output <- file(outputFile, "w")
	
	input <- file(inputFile, "r", blocking = FALSE)
	line <- readLines(input, n=1)
	while (substr(line,1,1)=="#") {
		# skip header comments but retain them in output
		writeLines(line, con=output)
		line <- readLines(input, n=1)	
	}
	fields <- unlist(strsplit(line, "\t"))
	close(input)
	close(output)
	
	rowCount <- dim(inputTable)[1]
	outputTable <- inputTable
	#pValueNames <- c("Fisher_p-value_tCw", "p-value_GvT_skew")
	pValueNames <- c("Fisher_p-value_tCw", "Fisher_p-value_tCw_to_G", "Fisher_p-value_tCw_to_T", "p-value_GvT_skew")
	if (CUSTOMMOTIFS) {
		# REMOVE skew decided 10/18/16
		#pValueNames <- c(paste("Fisher_p-value_", CUSTOMTYPE, sep=""), skewPvalueName)
		pValueNames <- paste("Fisher_p-value_", CUSTOMTYPE, sep="")
	}
	for (pValCol in pValueNames) {
		inputTable <- outputTable
		colCount <- dim(inputTable)[2]
		PVALUE_FIELD <- which(fields==pValCol)
	
		pValueCol <- inputTable[1:(rowCount-1),PVALUE_FIELD]
		col1 <- c(p.adjust(pValueCol, method="BH"), NA)
		col2 <- c(p.adjust(pValueCol, method="bonferroni"), NA)
		newColNames <- paste(c("BH_","bonferroni_"), pValCol, sep="")
		colShift <- 0
	
		if (PVALUE_FIELD!=colCount) {
			outputTable <- cbind(inputTable[,1:PVALUE_FIELD],col1,col2,inputTable[,(PVALUE_FIELD+1):colCount])
			outputTable <- cbind(inputTable[,1:PVALUE_FIELD],col1,col2,inputTable[,(PVALUE_FIELD+1):colCount])
			colnames(outputTable)[1:PVALUE_FIELD] <- fields[1:PVALUE_FIELD]
			colnames(outputTable)[(PVALUE_FIELD+1):(PVALUE_FIELD+2)] <- newColNames
		
			colnames(outputTable)[(PVALUE_FIELD+3):(colCount+2)] <- fields[(PVALUE_FIELD+1):colCount]
			colShift <- 2
			fields <- c(fields[1:PVALUE_FIELD], newColNames, fields[(PVALUE_FIELD+1):colCount])
		} else {
			outputTable <- cbind(inputTable[,1:PVALUE_FIELD],col1,col2)
			colnames(outputTable)[1:PVALUE_FIELD] <- fields[1:PVALUE_FIELD]
			colnames(outputTable)[(PVALUE_FIELD+1):(PVALUE_FIELD+2)] <- newColNames
		
			colShift <- 2
			fields <- c(fields[1:PVALUE_FIELD], newColNames)	
		}
	}
	
	# additional per row operations
	if (!CUSTOMMOTIFS) {
		PVALUE_FIELD <- which(fields=="BH_Fisher_p-value_tCw")
		# revcomp added 05/15/18 
		TCW_FIELD <- which(fields=="tCw_to_G+tCw_to_T+revcomp")
		#TCW_FIELD <- which(fields=="tCw_to_G+tCw_to_T")
		APOBEC_ENRICH_FIELD <- which(fields=="APOBEC_enrich")

		nrows <- dim(outputTable)[1]
		ncols <- dim(outputTable)[2]

		mutLoadCol <- rep(0, nrows)
		signifRows <- which(outputTable[PVALUE_FIELD] <= 0.05)
		mutLoadCol[signifRows] <- round(outputTable[signifRows, TCW_FIELD] * (outputTable[signifRows, APOBEC_ENRICH_FIELD]-1)/outputTable[signifRows, APOBEC_ENRICH_FIELD], digits=0)
		mutLoadCol[nrows] <- NA
		outputTable <- cbind(outputTable, mutLoadCol)
		colnames(outputTable)[ncols+1] <- "APOBEC_MutLoad_MinEstimate"
		apobecMutloadField <- ncols+1
		
		# added 05/19/16
		PVALUE_FIELD <- which(fields=="BH_Fisher_p-value_tCw_to_G")
		TCW_FIELD <- which(fields=="APOBEC_to_G")
		apobecToGfield <- TCW_FIELD
		APOBEC_ENRICH_FIELD <- which(fields=="tCw_to_G_enrich")

		nrows <- dim(outputTable)[1]
		ncols <- dim(outputTable)[2]

		mutLoadCol <- rep(0, nrows)
		signifRows <- which(outputTable[PVALUE_FIELD] <= 0.05)
		mutLoadCol[signifRows] <- round(outputTable[signifRows, TCW_FIELD] * (outputTable[signifRows, APOBEC_ENRICH_FIELD]-1)/outputTable[signifRows, APOBEC_ENRICH_FIELD], digits=0)
		mutLoadCol[nrows] <- NA
		outputTable <- cbind(outputTable, mutLoadCol)
		colnames(outputTable)[ncols+1] <- "APOBEC_to_G_MutLoad_MinEstimate"
		apobecToGmutloadField <- ncols+1
		
		PVALUE_FIELD <- which(fields=="BH_Fisher_p-value_tCw_to_T")
		TCW_FIELD <- which(fields=="APOBEC_to_T")
		apobecToTfield <- TCW_FIELD
		APOBEC_ENRICH_FIELD <- which(fields=="tCw_to_T_enrich")

		nrows <- dim(outputTable)[1]
		ncols <- dim(outputTable)[2]

		mutLoadCol <- rep(0, nrows)
		signifRows <- which(outputTable[PVALUE_FIELD] <= 0.05)
		mutLoadCol[signifRows] <- round(outputTable[signifRows, TCW_FIELD] * (outputTable[signifRows, APOBEC_ENRICH_FIELD]-1)/outputTable[signifRows, APOBEC_ENRICH_FIELD], digits=0)
		mutLoadCol[nrows] <- NA
		outputTable <- cbind(outputTable, mutLoadCol)
		colnames(outputTable)[ncols+1] <- "APOBEC_to_T_MutLoad_MinEstimate"
		ncols <- dim(outputTable)[2]
		apobecToTmutloadField <- ncols
		
		chisqVals <- outputTable[1:nrows-1, c("APOBEC_to_G_MutLoad_MinEstimate","APOBEC_to_T_MutLoad_MinEstimate")]
		testRows <- which(apply(chisqVals,1,sum)!=0)

		chisqCol <- rep(NA, nrows)
		if (length(testRows)>0) {
			chisqCol[testRows] <- apply(chisqVals[testRows,],1,getChisqPvalue)
		}
		outputTable <- cbind(outputTable, chisqCol)	
		colnames(outputTable)[ncols+1] <- "p_value_Chi-squared_APOBEC_to_G_vs_APOBEC_to_T"
		ncols <- dim(outputTable)[2]

		adjCol <- p.adjust(chisqCol, method="BH")
		outputTable <- cbind(outputTable, adjCol)	
		colnames(outputTable)[ncols+1] <- "BH_p-value_Chi-squared_APOBEC_to_G_vs_APOBEC_to_T"
		ncols <- dim(outputTable)[2]
		newPvalueField <- ncols
		
		adjCol <- p.adjust(chisqCol, method="bonferroni")
		outputTable <- cbind(outputTable, adjCol)	
		colnames(outputTable)[ncols+1] <- "bonferroni_p-value_Chi-squared_APOBEC_to_G_vs_APOBEC_to_T"
		
		if (FALSE) {
		newCol1 <- outputTable[,apobecToGfield] - outputTable[,apobecToTfield]
		outputTable <- cbind(outputTable, newCol1)
		colnames(outputTable)[ncols+2] <- "APOBEC_to_G-APOBEC_to_T"
		
		newCol2 <- rep(0, nrows)
		signifRows <- which(outputTable[newPvalueField] <= 0.05)
		newCol2[signifRows] <- newCol1[signifRows]
		outputTable <- cbind(outputTable, newCol2)
		colnames(outputTable)[ncols+3] <- "APOBEC_to_G-APOBEC_to_T_corr"
		
		newCol3 <- rep("ND", nrows)
		newCol3[newCol2>0] <- "APOBEC_to_G"
		newCol3[newCol2<0] <- "APOBEC_to_T"
		outputTable <- cbind(outputTable, newCol3)
		colnames(outputTable)[ncols+4] <- "APOBEC_to_G_or_APOBEC_to_T"

		newCol4 <- rep("ND", nrows)
		nonZeroRows <- which(outputTable[,apobecMutloadField]!=0)
		newCol4[nonZeroRows] <- newCol1[nonZeroRows]/outputTable[nonZeroRows,apobecMutloadField]
		outputTable <- cbind(outputTable, newCol4)
		colnames(outputTable)[ncols+5] <- "relative_APOBEC_to_G_vs_T_skew"
		} # END if (FALSE) 
		
		newCol1 <- outputTable[,apobecToGmutloadField] - outputTable[,apobecToTmutloadField]
		outputTable <- cbind(outputTable, newCol1)
		colnames(outputTable)[ncols+2] <- "APOBEC_to_G_MutLoad_MinEstimate-APOBEC_to_T_MutLoad_MinEstimate"
		
		if (FALSE) {
		newCol2 <- rep(0, nrows)
		#signifRows <- which(outputTable[,apobecMutloadField] > 0)
		signifRows <- which(outputTable[newPvalueField] <= 0.05)
		newCol2[signifRows] <- newCol1[signifRows]
		outputTable <- cbind(outputTable, newCol2)
		colnames(outputTable)[ncols+3] <- "APOBEC_to_G_MutLoad_MinEstimate-APOBEC_to_T_MutLoad_MinEstimate_corr"
		}
		
		newCol2 <- rep("ND", nrows)
		signifRows <- which((outputTable[,newPvalueField] <= 0.05) & (outputTable[,apobecMutloadField]!=0))
		newCol2[signifRows] <- newCol1[signifRows]
		newCol2[nrows] <- "NA"
		outputTable <- cbind(outputTable, newCol2)
		colnames(outputTable)[ncols+3] <- "APOBEC_to_G_MutLoad_MinEstimate-APOBEC_to_T_MutLoad_MinEstimate_corr"
		
		
		
		newCol3 <- rep("ND", nrows)
		newCol3[as.numeric(newCol2)>0] <- "APOBEC_to_G"
		newCol3[as.numeric(newCol2)<0] <- "APOBEC_to_T"
		newCol3[nrows] <- "NA"
		outputTable <- cbind(outputTable, newCol3)
		colnames(outputTable)[ncols+4] <- "APOBEC_to_G_or_APOBEC_to_T"

		if (FALSE) {
		newCol4 <- rep("ND", nrows)
		nonZeroRows <- which((outputTable[,apobecToGmutloadField]+outputTable[,apobecToTmutloadField])!=0)
		newCol4[nonZeroRows] <- newCol2[nonZeroRows]/(outputTable[nonZeroRows,apobecToGmutloadField]+outputTable[nonZeroRows,apobecToTmutloadField])
		moreNDrows <- which(newCol3=="ND")
		newCol4[moreNDrows] <- "ND"
		outputTable <- cbind(outputTable, newCol4)
		colnames(outputTable)[ncols+5] <- "relative_APOBEC_to_G_vs_T_diff"
		}
		
		newCol4 <- rep("ND", nrows)
		signifRows <- which(newCol3!="ND")
		newCol4[signifRows] <- as.numeric(newCol2[signifRows])/(outputTable[signifRows,apobecToGmutloadField]+outputTable[signifRows,apobecToTmutloadField])
		newCol4[nrows] <- "NA"
		outputTable <- cbind(outputTable, newCol4)
		colnames(outputTable)[ncols+5] <- "relative_APOBEC_to_G_vs_T_diff"

	# end if (!CUSTOMMOTIFS)		
	} else {
		if (CUSTOMMUTLOAD) {
			PVALUE_FIELD <- which(fields==paste("BH_Fisher_p-value_", CUSTOMTYPE, sep=""))
			MUT_FIELD <- which(fields==MUTFIELDNAME)
			if (length(MUT_FIELD)==2) {
				MUT_FIELD <- MUT_FIELD[2]
			}
			MUT_ENRICH_FIELD <- which(fields==paste(CUSTOMTYPE, "_enrich", sep=""))

			nrows <- dim(outputTable)[1]
			ncols <- dim(outputTable)[2]

			mutLoadCol <- rep(0, nrows)
			signifRows <- which(outputTable[PVALUE_FIELD] <= 0.05)
			mutLoadCol[signifRows] <- round(outputTable[signifRows, MUT_FIELD] * (outputTable[signifRows, MUT_ENRICH_FIELD]-1)/outputTable[signifRows, MUT_ENRICH_FIELD], digits=0)
			mutLoadCol[nrows] <- NA
			outputTable <- cbind(outputTable, mutLoadCol)
			colnames(outputTable)[ncols+1] <- paste(CUSTOMTYPE, "_MutLoad_MinEstimate", sep="")
			
		}
	}
	
	write.table(outputTable, outputFile, col.names=TRUE, row.names=F, sep="\t", append=TRUE)

}

