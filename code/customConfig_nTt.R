## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "nTt"
ANALYZECOLS <- c("nTt", "aAn")

# instroduced with the UV pipeline
# don't use sampleTitles from addSampleCols.R use those from the config
CUSTOMMUTLOAD <- TRUE
# MUTFIELDNAME should agree with custFields No.2
#MUTFIELDNAME <- "nTt_to_C+aAn_to_G"
#MUTFIELDNAME <- "nTt_to_C"
MUTFIELDNAME <- "nTt_to_C+revcomp"
# "[nTt_to_C+aAn_to_G]_per_mut"
# "[nTt_to_C+revcomp]_per_mut"
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
# [nTt_to_C+aAn_to_G]_per_mut
# [nTt_to_C+revcomp]_per_mut
(nTt_to_C + aAn_to_G) / mutations, 
# nTt_to_C+aAn_to_G
# nTt_to_C+revcomp
nTt_to_C+aAn_to_G, 
# [(T_to_C)]-[(nTt_to_C)]
# [(T_to_C)]-[(nTt_to_C)]+revcomp
(T_to_C + A_to_G) - (nTt_to_C+aAn_to_G),
# ntt+aan
# ntt+revcompaan
ntt + aan, 
# t-ntt
# t-ntt+revcomp
(t + a) - (ntt + aan)
         )

         testMat <- rbind(c((nTt_to_C + aAn_to_G), (T_to_C + A_to_G - nTt_to_C - aAn_to_G)), c(ntt + aan, t+a-ntt-aan))
         ft <- fisher.test(testMat, alternative="greater")

         testFields <- c(
# Fisher_p-value_nTt
ft$p.value,
# odds_ratio_nTt
ft$estimate,
#paste(ft$conf.int, collapse="-")
# 95%_confidence_lower_nTt
ft$conf.int[1],
# 95%_confidence_upper_nTt
ft$conf.int[2]
         )

         custFields2 <- c(
# nTt_enrich
((nTt_to_C + aAn_to_G) * (t + a))/((T_to_C + A_to_G ) * (ntt + aan))
				 )
                                					
# REMOVE skew decided 10/18/16
                                					
          #testMat2 <- rbind(c((nTt_to_G + aAn_to_C), (nTt + aAn - nTt_to_G - aAn_to_C)), c((nTt_to_C + aAn_to_G), (nTt + aAn - nTt_to_C - aAn_to_G)))
					
          #ft <- fisher.test(testMat2, alternative="two.sided")
        #testFields2 <- c(
# p-value_GvC_skew
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
#CUSTOMSKEW <- TRUE
#SKEWTYPE <- "GvC"
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
"[(T_to_C)]-[(nTt_to_C)]+revcomp",
"ntt+revcomp",
"t-ntt+revcomp",
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



