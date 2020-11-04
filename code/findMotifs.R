## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

if (!CUSTOMMOTIFS) {
motifs2Find <- c("A", "T", "G", "C", "Cg", "cG", "tC[at]", "[at]Ga", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tC[atc]", "[atg]Ga", "cC", "Gg", "[at][ag]C", "G[ct][at]", "Cc", "gG", "[at]A", "T[at]")
findTitles <- c("A", "T", "G", "C", "Cg", "cG", "tCw", "wGa", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tCh", "dGa", "cC", "Gg", "wrC", "Gyw", "Cc", "gG", "wA", "Tw")
#rbind(motifs2Find, findTitles)
}

#if (CUSTOMMOTIFS) {
if (FALSE) {
	customFindScript <- paste("findMotifs_", CUSTOMTYPE, ".R", sep="")
	source(customFindScript)
	#source("findMotifsCustom.R")
}

apobecTitles <- c("APOBEC_mutation", "APOBEC_mutation_to_G", "APOBEC_mutation_to_T")
# modified 10/28/14
apobecTitles <- c("tC_mutation", "tC_mutation_to_G", "tC_mutation_to_T", "APOBEC_mutation", "APOBEC_mutation_to_G", "APOBEC_mutation_to_T")

motifPresent <- rep(0, length(motifs2Find))

findMotifs <- function(seq) {
	c <- 0
	for (motif in motifs2Find) {
		c <- c + 1
		#for (i in (1:2)) {
		#}
		if (motif=="+") {
			motifPresent[c] <- motifPresent[c-1] + motifPresent[c-2]
		} else {
			match <- grep(motif, seq, perl=TRUE)
			if (length(match)!=0) {
				motifPresent[c] <- 1
				next
			}
		}
	}

	return(paste(motifPresent, collapse="\t"))
}

tCwPos <- which(findTitles=="tCw")
wGaPos <- which(findTitles=="wGa")
# modified 10/28/14 to add 3 tC_mutation columns
tCPos <- which(findTitles=="tC")
GaPos <- which(findTitles=="Ga")
isApobec <- function(findString) {
	apobecBits <- c("0", "0", "0", "0", "0", "0")
	findBits <- unlist(strsplit(findString, "\t"))
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	) {
		apobecBits[5] <- "1"
	}
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	) {
		apobecBits[6] <- "1"
	}
	if((apobecBits[5]=="1") | (apobecBits[6]=="1")) apobecBits[4] <- "1"
	
	if((apobecBits[5]=="1") |
	((findBits[tCPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[GaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	) {
		apobecBits[2] <- "1"
	}
	if((apobecBits[6]=="1") |
	((findBits[tCPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[GaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	) {
		apobecBits[3] <- "1"
	}
	if((apobecBits[2]=="1") | (apobecBits[3]=="1")) apobecBits[1] <- "1"
	return(paste(apobecBits, collapse="\t"))
}
# version before 10/28/14
isApobecOld <- function(findString) {
	apobecBits <- c("0", "0", "0")
	findBits <- unlist(strsplit(findString, "\t"))
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	) {
		apobecBits[2] <- "1"
	}
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	) {
		apobecBits[3] <- "1"
	}
	if((apobecBits[2]=="1") | (apobecBits[3]=="1")) apobecBits[1] <- "1"
	return(paste(apobecBits, collapse="\t"))
}

