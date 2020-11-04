## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "rTt"
ANALYZECOLS <- c("rTt", "aAy")

# introduced with the UV pipeline
CUSTOMMUTLOAD <- TRUE
# MUTFIELDNAME should agree with custFields No.2
# don't list revcomps as in tCw, to be fixed in all versions by adding "+revcomp"
#MUTFIELDNAME <- "rTt_to_C+aAy_to_G"
MUTFIELDNAME <- "rTt_to_C"
# "[rTt_to_C+aAy_to_G]_per_mut"
perMutFieldName <- paste("[", MUTFIELDNAME, "]_per_mut", sep="")

# This file should be saved as customConfig_CUSTOMTYPE.R (e.g. customConfig_htCw.R) and accompanied by a sample file with column names for analytical formulas named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) - see last section

# Motif definitions:

# for findMotifs.R
motifs2Find <- c("A", "T", "G", "C", "Cg", "cG", "[ag]Tt", "aA[tc]")
findTitles <- c("A", "T", "G", "C", "Cg", "cG", "rTt", "aAy")

# for countMotifs.R
motifs2Count <- c("a", "t", "g", "c", "cg", "[ag]tt", "aa[tc]")
countTitles <- c("a", "t", "g", "c", "cg", "rtt", "aay")

# for analyzeAutoExt.R

# Analysis module formulas:

getCustomFields <- function() {
  custFields <- c(
    # [rTt_to_C+aAy_to_G]_per_mut
    (rTt_to_C + aAy_to_G) / mutations, 
    # rTt_to_C+aAy_to_G
    rTt_to_C+aAy_to_G, 
    # [(T_to_C)]-[(rTt_to_C)]
    (T_to_C + A_to_G) - (rTt_to_C+aAy_to_G),
    #rtt+aay
    rtt + aay, 
    # t-rtt
    (t + a) - (rtt + aay)
  )
  
  testMat <- rbind(c((rTt_to_C + aAy_to_G), (T_to_C + A_to_G - rTt_to_C - aAy_to_G)), c(rtt + aay, t+a-rtt-aay))
  ft <- fisher.test(testMat, alternative="greater")
  
  testFields <- c(
    # Fisher_p-value_rTt
    ft$p.value,
    # odds_ratio_rTt
    ft$estimate,
    #paste(ft$conf.int, collapse="-")
    # 95%_confidence_lower_rTt
    ft$conf.int[1],
    # 95%_confidence_upper_rTt
    ft$conf.int[2]
  )
  
  custFields2 <- c(
    # rTt_enrich
    ((rTt_to_C + aAy_to_G) * (t + a))/((T_to_C + A_to_G ) * (rtt + aay))
  )
  
# 09/07/18 correction
# file pre-dates
# ls -l customConfig_rTt.R
# -rw-r--r-- 1 klimczakl 12140 4179 Aug 26  2016 customConfig_rTt.R
# REMOVE skew decided 10/18/16

  #testMat2 <- rbind(c((rTt_to_G + aAy_to_C), (rTt + aAy - rTt_to_G - aAy_to_C)), c((rTt_to_C + aAy_to_G), (rTt + aAy - rTt_to_C - aAy_to_G)))
  
  #ft <- fisher.test(testMat2, alternative="two.sided")
  #testFields2 <- c(
    # p-value_GvC_skew
    #ft$p.value
  #)
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
#  SKEWTYPE <- "GvT"
#}
#skewPvalueName <- paste("p-value_", SKEWTYPE, "_skew", sep="")


# all the column names that appear in the above formulas as comments should be grouped here
# this will replace the fisherSample file
newFieldTitles <- c(
  perMutFieldName,
  MUTFIELDNAME,
  ## EDIT START 3 ##
  "[(T_to_C)]-[(rTt_to_C)]",
  "rtt+aay",
  "t-rtt",
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


