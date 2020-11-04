## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

if (not CUSTOMMOTIFS):
	motifs2Find = ("A", "T", "G", "C", "Cg", "cG", "tC[at]", "[at]Ga", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tC[atc]", "[atg]Ga", "cC", "Gg", "[at][ag]C", "G[ct][at]", "Cc", "gG", "[at]A", "T[at]")
	findTitles = ("A", "T", "G", "C", "Cg", "cG", "tCw", "wGa", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tCh", "dGa", "cC", "Gg", "wrC", "Gyw", "Cc", "gG", "wA", "Tw")

# modified 10/28/14
apobecTitles = ("tC_mutation", "tC_mutation_to_G", "tC_mutation_to_T", "APOBEC_mutation", "APOBEC_mutation_to_G", "APOBEC_mutation_to_T")

def findMotifs(seq):
	motifPresent = ["0"] * len(motifs2Find)
	c = 0
	for motif in motifs2Find:
		if (motif=="+"):
			motifPresent[c] = motifPresent[c-1] + motifPresent[c-2]
		else:
			if re.search(motif, seq):
				motifPresent[c] = "1"
				next
		c = c + 1

	return reduce(lambda x, y: x+'\t'+y, motifPresent)


tCwPos = which([x=="tCw" for x in findTitles])[0]
wGaPos = which([x=="wGa" for x in findTitles])[0]
# modified 10/28/14 to add 3 tC_mutation columns
tCPos = which([x=="tC" for x in findTitles])[0]
GaPos = which([x=="Ga" for x in findTitles])[0]
def isApobec(findString):
	apobecBits = ["0", "0", "0", "0", "0", "0"]
	findBits = findString.split("\t")
	if (
	((findBits[tCwPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[wGaPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	):
		apobecBits[4] = "1"

	if (
	((findBits[tCwPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[wGaPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	):
		apobecBits[5] = "1"

	if((apobecBits[4]=="1") or (apobecBits[5]=="1")):
		apobecBits[3] = "1"
	
	if((apobecBits[4]=="1") or
	((findBits[tCPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[GaPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	):
		apobecBits[1] = "1"

	if((apobecBits[5]=="1") or
	((findBits[tCPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[GaPos]=="1") and fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	):
		apobecBits[2] = "1"

	if((apobecBits[1]=="1") | (apobecBits[2]=="1")):
		apobecBits[0] = "1"

	return list2tabString(apobecBits)

# version before 10/28/14
'''
isApobecOld = function(findString) {
	apobecBits = c("0", "0", "0")
	findBits = unlist(strsplit(findString, "\t"))
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="G") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="C")
	) {
		apobecBits[2] = "1"
	}
	if(
	((findBits[tCwPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="T") |
	((findBits[wGaPos]=="1") & fields[TUMOR_SEQ_ALLELE2_FIELD]=="A")
	) {
		apobecBits[3] = "1"
	}
	if((apobecBits[2]=="1") | (apobecBits[3]=="1")) apobecBits[1] = "1"
	return(paste(apobecBits, collapse="\t"))
}
'''
