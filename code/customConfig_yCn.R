## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "yCn"
ANALYZECOLS <- c("yCn", "nGr")

# introduced with the UV pipeline
# don't use sampleTitles from addSampleCols.R use those from the config
CUSTOMMUTLOAD <- TRUE
# MUTFIELDNAME should agree with custFields No.2
#MUTFIELDNAME <- "yCn_to_T+nGr_to_A"
#MUTFIELDNAME <- "yCn_to_T"
MUTFIELDNAME <- "yCn_to_T+revcomp"
# [yCn_to_T+nGr_to_A]_per_mut
# [yCn_to_T+revcomp]_per_mut
perMutFieldName <- paste("[", MUTFIELDNAME, "]_per_mut", sep="")

# This file should be saved as customConfig_CUSTOMTYPE.R (e.g. customConfig_htCw.R) and accompanied by a sample file with column names for analytical formulas named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) - see last section

# Motif definitions:

# for findMotifs.R
motifs2Find <- c("A", "T", "G", "C", "Cg", "cG", "[atgc]Tt", "aA[atgc]", "[ct]C[atgc]", "[atgc]G[ga]", "[atgc]Cg", "cG[atgc]")
findTitles <- c("A", "T", "G", "C", "Cg", "cG", "nTt", "aAn", "yCn", "nGr", "nCg", "cGn")

# for countMotifs.R
motifs2Count <- c("a", "t", "g", "c", "cg", "[atgc]tt", "aa[atgc]", "[ct]c[atgc]", "[atgc]g[ga]", "[atgc]cg", "cg[atgc]")
countTitles <- c("a", "t", "g", "c", "cg", "ntt", "aan", "ycn", "ngr", "ncg", "cgn")

# for analyzeAutoExt.R

# Analysis module formulas:

getCustomFields <- function() {
         custFields <- c(
# [yCn_to_T+nGr_to_A]_per_mut
(yCn_to_T + nGr_to_A) / mutations, 
# yCn_to_T + nGr_to_A 
yCn_to_T + nGr_to_A, 
# [(C_to_T)]-[(yCn_to_T)]
(C_to_T + G_to_A) - (yCn_to_T + nGr_to_A), 
#ycn + ngr 
ycn + ngr, 
# c-ycn 
(c + g) - (ycn + ngr)
         )


        testMat <- rbind(c((yCn_to_T + nGr_to_A), (C_to_T + G_to_A - yCn_to_T - nGr_to_A)),c(ycn+ngr, c+g-ycn-ngr))
        ft <- fisher.test(testMat, alternative="greater")


testFields <- c(
# Fisher_p-value_yCn
ft$p.value,
# odds_ratio_yCn
ft$estimate,
#paste(ft$conf.int, collapse="-")
# 95%_confidence_lower_yCn
ft$conf.int[1],
# 95%_confidence_upper_yCn
ft$conf.int[2]
         )


        custFields2 <- c(
# yCn_enrich
((yCn_to_T + nGr_to_A) * (c + g))/((C_to_T + G_to_A ) * (ycn + ngr))
				)


# REMOVE skew decided 10/18/16
                                					
          #testMat2 <- rbind(c((nCg_to_G+cGn_to_C), (nCg + cGn - nCg_to_G - cGn_to_C)), c((nCg_to_T+cGn_to_A), (nCg + cGn - nCg_to_T - cGn_to_A)))
					#testMat2 <- rbind(c((yCn_to_G + nGr_to_C), (yCn + nGr - yCn_to_G - nGr_to_C)), c((yCn_to_T + nGr_to_A), (yCn + nGr - yCn_to_T - nGr_to_A)))
					
          #ft <- fisher.test(testMat2, alternative="two.sided")
        #testFields2 <- c(
# p-value_GvT_skew
#ft$p.value
#        )
## EDIT END 1 ##

	#newFields <- c(custFields,testFields,custFields2,testFields2)
	newFields <- c(custFields,testFields,custFields2)
	return(newFields)
}

# REMOVE skew decided 10/18/16
## EDIT START 2 ##
# "GvT" is default, use CUSTOMSKEW <- FALSE
#CUSTOMSKEW <- FALSE
# for any other skew types
#CUSTOMSKEW <- FALSE
#SKEWTYPE <- "GvT"
## EDIT END 2 ##
#if (!CUSTOMSKEW) {
#SKEWTYPE <- "GvT"
#}
#skewPvalueName <- paste("p-value_", SKEWTYPE, "_skew", sep="")


# all the column names that appear in the above formulas as comments should be grouped here
# this will replace the fisherSample file
newFieldTitles <- c(
perMutFieldName,
MUTFIELDNAME,
## EDIT START 3 ##
"[(C_to_T)]-[(yCn_to_T)]+revcomp",
"ycn+revcomp",
"c-ycn+revcomp",
## EDIT END 3 ##
paste("Fisher_p-value_", CUSTOMTYPE, sep=""),
paste("odds_ratio_", CUSTOMTYPE, sep=""),
paste("95%_confidence_lower_", CUSTOMTYPE, sep=""),
paste("95%_confidence_upper_", CUSTOMTYPE, sep=""),
paste(CUSTOMTYPE, "_enrich", sep="")
#skewPvalueName
)


# The approach described in the comment below will be replaced 
# by newFieldTitles and a column retention procedure
# A sample file named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) must contain the following types of CUSTOMTYPE columns; this is just for illustration, do not edit here but in fisherSample


# for addSampleCols.R

sampleTitles <- c(
perMutFieldName,
MUTFIELDNAME,
paste("BH_Fisher_p-value_", CUSTOMTYPE, sep=""),
paste(CUSTOMTYPE, "_enrich", sep="")
#skewPvalueName,
#paste("BH_", skewPvalueName, sep="")
)

if (CUSTOMMUTLOAD) {
sampleTitles <- c(sampleTitles,
paste(CUSTOMTYPE, "_MutLoad_MinEstimate", sep="")
)
}



