## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "yTt"
ANALYZECOLS <- c("yTt", "aAr")

# introduced with the UV pipeline
CUSTOMMUTLOAD <- TRUE
# MUTFIELDNAME should agree with custFields No.2
# don't list revcomps as in tCw, to be fixed in all versions by adding "+revcomp"
#MUTFIELDNAME <- "yTt_to_C+aAr_to_G"
MUTFIELDNAME <- "yTt_to_C+revcomp"
# "[yTt_to_C+aAr_to_G]_per_mut"
perMutFieldName <- paste("[", MUTFIELDNAME, "]_per_mut", sep="")

# This file should be saved as customConfig_CUSTOMTYPE.R (e.g. customConfig_htCw.R) and accompanied by a sample file with column names for analytical formulas named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) - see last section

# Motif definitions:

# for findMotifs.R
motifs2Find <- c("A", "T", "G", "C", "Cg", "cG", "[tc]Tt", "aA[ag]")
findTitles <- c("A", "T", "G", "C", "Cg", "cG", "yTt", "aAr")

# for countMotifs.R
motifs2Count <- c("a", "t", "g", "c", "cg", "[tc]tt", "aa[ag]")
countTitles <- c("a", "t", "g", "c", "cg", "ytt", "aar")

# for analyzeAutoExt.R

# Analysis module formulas:

getCustomFields <- function() {
  custFields <- c(
    # [yTt_to_C+aAr_to_G]_per_mut
    (yTt_to_C + aAr_to_G) / mutations, 
    # yTt_to_C+aAr_to_G
    yTt_to_C+aAr_to_G, 
    # [(T_to_C)]-[(yTt_to_C)]
    (T_to_C + A_to_G) - (yTt_to_C+aAr_to_G),
    #ytt+aar
    ytt + aar, 
    # t-ytt
    (t + a) - (ytt + aar)
  )
  
  testMat <- rbind(c((yTt_to_C + aAr_to_G), (T_to_C + A_to_G - yTt_to_C - aAr_to_G)), c(ytt + aar, t+a-ytt-aar))
  ft <- fisher.test(testMat, alternative="greater")
  
  testFields <- c(
    # Fisher_p-value_yTt
    ft$p.value,
    # odds_ratio_yTt
    ft$estimate,
    #paste(ft$conf.int, collapse="-")
    # 95%_confidence_lower_yTt
    ft$conf.int[1],
    # 95%_confidence_upper_yTt
    ft$conf.int[2]
  )
  
  custFields2 <- c(
    # yTt_enrich
    ((yTt_to_C + aAr_to_G) * (t + a))/((T_to_C + A_to_G ) * (ytt + aar))
  )
  
  testMat2 <- rbind(c((yTt_to_G + aAr_to_C), (yTt + aAr - yTt_to_G - aAr_to_C)), c((yTt_to_C + aAr_to_G), (yTt + aAr - yTt_to_C - aAr_to_G)))
  
  ft <- fisher.test(testMat2, alternative="two.sided")
  testFields2 <- c(
    # p-value_GvC_skew
    ft$p.value
  )
  ## EDIT END 1 ##
  
  newFields <- c(custFields,testFields,custFields2,testFields2)
  return(newFields)
}

## EDIT START 2 ##
# "GvT" is default, use CUSTOMSKEW <- FALSE
#CUSTOMSKEW <- FALSE
# for any other skew types
CUSTOMSKEW <- TRUE
SKEWTYPE <- "GvC"
## EDIT END 2 ##
if (!CUSTOMSKEW) {
  SKEWTYPE <- "GvT"
}
skewPvalueName <- paste("p-value_", SKEWTYPE, "_skew", sep="")


# all the column names that appear in the above formulas as comments should be grouped here
# this will replace the fisherSample file
newFieldTitles <- c(
  perMutFieldName,
  MUTFIELDNAME,
  ## EDIT START 3 ##
  "[(T_to_C)]-[(yTt_to_C)]",
  "ytt+aar",
  "t-ytt",
  ## EDIT END 3 ##
  paste("Fisher_p-value_", CUSTOMTYPE, sep=""),
  paste("odds_ratio_", CUSTOMTYPE, sep=""),
  paste("95%_confidence_lower_", CUSTOMTYPE, sep=""),
  paste("95%_confidence_upper_", CUSTOMTYPE, sep=""),
  paste(CUSTOMTYPE, "_enrich", sep=""),
  skewPvalueName
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


