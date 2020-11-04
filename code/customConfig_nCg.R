## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "nCg"
ANALYZECOLS <- c("nCg", "cGn")
# introduced with the UV pipeline
CUSTOMMUTLOAD <- TRUE
# MUTFIELDNAME should agree with custFields No.2
#MUTFIELDNAME <- "nCg_to_T+cGn_to_A"
#MUTFIELDNAME <- "nCg_to_T"
MUTFIELDNAME <- "nCg_to_T+revcomp"
# "[nCg_to_T+cGn_to_A]_per_mut"
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
# [nCg_to_T+cGn_to_A]_per_mut
(nCg_to_T+cGn_to_A) / mutations, 
# nCg_to_T+cGn_to_A
nCg_to_T + cGn_to_A, 
# [(C_to_T)]-[(nCg_to_T)]
(C_to_T + G_to_A) - (nCg_to_T+cGn_to_A), 
#ncg+cgn 
ncg+cgn, 
# c-ncg 
(c + g) - (ncg+cgn)
         )


       
         testMat <- rbind(c((nCg_to_T+cGn_to_A), (C_to_T + G_to_A - nCg_to_T - cGn_to_A)),c(ncg+cgn, c+g-ncg-cgn))
         ft <- fisher.test(testMat, alternative="greater")

 testFields <- c(
# Fisher_p-value_nCg
ft$p.value,
# odds_ratio_nCg
ft$estimate,
#paste(ft$conf.int, collapse="-")
# 95%_confidence_lower_nCg
ft$conf.int[1],
# 95%_confidence_upper_nCg
ft$conf.int[2]
         )

         custFields2 <- c(
# nCg_enrich
((nCg_to_T+cGn_to_A) * (c + g))/((C_to_T + G_to_A ) * (ncg+cgn)))

                                					
# REMOVE skew decided 10/18/16
                                					
          #testMat2 <- rbind(c((nCg_to_G+cGn_to_C), (nCg + cGn - nCg_to_G - cGn_to_C)), c((nCg_to_T+cGn_to_A), (nCg + cGn - nCg_to_T - cGn_to_A)))
					
          #ft <- fisher.test(testMat2, alternative="two.sided")
        #testFields2 <- c(
# p-value_TvG_skew
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
#SKEWTYPE <- "TvG"
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
"[(C_to_T)]-[(nCg_to_T)]+revcomp",
"ncg+revcomp",
"c-ncg+revcomp",
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



