## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

## EDIT START 1 ##

CUSTOMTYPE <- "ytCa"
ANALYZECOLS <- c("ytCa", "tGar")
# instroduced with the UV pipeline
#CUSTOMMUTLOAD <- TRUE
CUSTOMMUTLOAD <- FALSE
# MUTFIELDNAME should agree with custFields2
MUTFIELDNAME <- "ytCa_to_G+ytCa_to_T"
# "[nTt_to_C+aAn_to_G]_per_mut"
perMutFieldName <- paste("[", MUTFIELDNAME, "]_per_mut", sep="")

# This file should be saved as customConfig_CUSTOMTYPE.R (e.g. customConfig_htCw.R) and accompanied by a sample file with column names for analytical formulas named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) - see last section

# Motif definitions:

# for findMotifs.R
motifs2Find <- c("A", "T", "G", "C", "Cg", "cG", "[ct]tCa", "tGa[ag]", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tC[atc]", "[atg]Ga", "cC", "Gg", "[at][ag]C", "G[ct][at]", "Cc", "gG", "[at]A", "T[at]")
findTitles <- c("A", "T", "G", "C", "Cg", "cG", "ytCa", "tGar", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tCh", "dGa", "cC", "Gg", "wrC", "Gyw", "Cc", "gG", "wA", "Tw")
#rbind(motifs2Find, findTitles)

# for countMotifs.R
motifs2Count <- c("a", "t", "g", "c", "cg", "[ct]tca", "tga[ag]", "tca", "tga", "tct", "aga", "tc", "ga", "tc[atc]", "[atg]ga", "cc", "gg", "[at][ag]c", "g[ct][at]", "cc", "gg", "[at]a", "t[at]")
countTitles <- c("a", "t", "g", "c", "cg", "ytca", "tgar", "tca", "tga", "tct", "aga", "tc", "ga", "tch", "dga", "cc", "gg", "wrc", "gyw", "cc", "gg", "wa", "tw")
#rbind(motifs2Count,countTitles)

# for analyzeAutoExt.R

# Analysis module formulas:

getCustomFields <- function() {

                                custFields <- c(
(ytCa_to_G + ytCa_to_T + tGar_to_C + tGar_to_A) / mutations,
ytCa_to_G + ytCa_to_T + tGar_to_C + tGar_to_A,
(C_to_G + C_to_T + G_to_C + G_to_A) - (ytCa_to_G + ytCa_to_T + tGar_to_C + tGar_to_A),
ytca + tgar,
(c + g) - (ytca + tgar)
                                )

                                testMat <- rbind(c((ytCa_to_G + ytCa_to_T + tGar_to_C + tGar_to_A), (C_to_G + C_to_T + G_to_C + G_to_A - ytCa_to_G - ytCa_to_T - tGar_to_C - tGar_to_A)),c(ytca+tgar, c+g-ytca-tgar))
                                ft <- fisher.test(testMat, alternative="greater")
                                
                                testFields <- c(
ft$p.value,
ft$estimate,
ft$conf.int[1],
ft$conf.int[2]
                                )
                                
                                custFields2 <- c(
((ytCa_to_G + ytCa_to_T + tGar_to_C + tGar_to_A) * (c + g))/((C_to_G + C_to_T + G_to_C + G_to_A) * (ytca + tgar)),
((ytCa_to_G + tGar_to_C) * (c + g))/((C_to_G + G_to_C) * (ytca + tgar)),
((ytCa_to_T + tGar_to_A) * (c + g))/((C_to_T + G_to_A) * (ytca + tgar))
                                )
# REMOVE skew decided 10/18/16																
                                #testMat2 <- rbind(c((ytCa_to_G + tGar_to_C), (ytCa + tGar - ytCa_to_G - tGar_to_C)), c((ytCa_to_T + tGar_to_A), (ytCa + tGar - ytCa_to_T - tGar_to_A)))
                                #ft <- fisher.test(testMat2, alternative="two.sided")                          
                                #testFields2 <- c(
#ft$p.value
                                #)
                                
# non-APOBEC_mutations
# non-APOBEC_substitutions                             
                                custFields3 <- c(mutations - ytCa_to_G - ytCa_to_T - tGar_to_C - tGar_to_A, substitutions - ytCa_to_G - ytCa_to_T - tGar_to_C - tGar_to_A)
## EDIT END 1 ##

	#newFields <- c(custFields,testFields,custFields2,testFields2,custFields3)
	newFields <- c(custFields,testFields,custFields2,custFields3)
	return(newFields)
}


# REMOVE skew decided 10/18/16
## EDIT START 2 ##
# "GvT" is default, use CUSTOMSKEW <- FALSE
#CUSTOMSKEW <- FALSE
# for any other skew types
#CUSTOMSKEW <- TRUE
#SKEWTYPE <- "GvT"
## EDIT END 2 ##
#if (CUSTOMSKEW) {
#SKEWTYPE <- "GvC"
#}
#skewPvalueName <- paste("p-value_", SKEWTYPE, "_skew", sep="")


# all the column names that appear in the above formulas as comments should be grouped here
# this will replace the fisherSample file
newFieldTitles <- c(
perMutFieldName,
MUTFIELDNAME,
## EDIT START 3 ##
"[(C_to_G)+(C_to_T)]-[(ytCa_to_G)+(ytCa_to_T)]",
"ytca+tgar",
"c-ytca",
## EDIT END 3 ##
paste("Fisher_p-value_", CUSTOMTYPE, sep=""),
paste("odds_ratio_", CUSTOMTYPE, sep=""),
paste("95%_confidence_lower_", CUSTOMTYPE, sep=""),
paste("95%_confidence_upper_", CUSTOMTYPE, sep=""),
# APOBECtCa_enrich
paste("APOBEC", CUSTOMTYPE, "_enrich", sep=""),
paste(CUSTOMTYPE, "_to_G_enrich", sep=""),
paste(CUSTOMTYPE, "_to_T_enrich", sep=""),
#skewPvalueName,
paste("non-APOBEC", CUSTOMTYPE, "_mutations", sep=""),
paste("non-APOBEC", CUSTOMTYPE, "_substitutions", sep="")
)


# The approach described in the comment below will be replaced 
# by newFieldTitles and a column retention procedure
# A sample file named fisherSample_CUSTOMTYPE.txt (e.g. fisherSample_htCw.txt) must contain the following types of CUSTOMTYPE columns; this is just for illustration, do not edit here but in fisherSample

# for addSampleCols.R

#sampleTitles <- c("[tCa_to_G+tCa_to_T]_per_mut",
#"tCa_to_G+tCa_to_T", 
#"BH_Fisher_p-value_tCa", 
#"APOBECtCa_enrich", 
#"tCa_to_G_enrich",
#"tCa_to_T_enrich",
#"p-value_GvT_skew",
#"BH_p-value_GvT_skew")

sampleTitles <- c(
perMutFieldName,
MUTFIELDNAME,
paste("BH_Fisher_p-value_", CUSTOMTYPE, sep=""),
paste("APOBEC", CUSTOMTYPE, "_enrich", sep=""),
paste(CUSTOMTYPE, "_to_G__enrich", sep=""),
paste(CUSTOMTYPE, "_to_T_enrich", sep="")
#skewPvalueName,
#paste("BH_", skewPvalueName, sep="")
)

if (CUSTOMMUTLOAD) {
sampleTitles <- c(sampleTitles,
paste(CUSTOMTYPE, "_MutLoad_MinEstimate", sep="")
)
}



