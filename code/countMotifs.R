## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

if (!CUSTOMMOTIFS) {
motifs2Count <- c("a", "t", "g", "c", "cg", "tc[at]", "[at]ga", "tca", "tga", "tct", "aga", "tc", "ga", "tc[atc]", "[atg]ga", "cc", "gg", "[at][ag]c", "g[ct][at]", "cc", "gg", "[at]a", "t[at]")
countTitles <- c("a", "t", "g", "c", "cg", "tcw", "wga", "tca", "tga", "tct", "aga", "tc", "ga", "tch", "dga", "cc", "gg", "wrc", "gyw", "cc", "gg", "wa", "tw")
}

#if (CUSTOMMOTIFS) {
if (FALSE) {
	customCountScript <- paste("countMotifs_", CUSTOMTYPE, ".R", sep="")
	source(customCountScript)
	#source("countMotifsCustom.R")
}
#rbind(motifs2Count, countTitles)
#lapply(lapply(list(motifs2Find=motifs2Find,findTitles=findTitles,motifs2Count=motifs2Count,countTitles=countTitles), duplicated), which)

if (CONTEXTCOUNT) {
	countTitles <<- paste(countTitles, "_counts", sep="")
}
if (CLUSTERCOUNT) {
	countTitles <<- paste(countTitles, "_in cluster", sep="")
}

#cbind(motifTitles, motifs)

motifCount <- rep(0, length(motifs2Count))

countMotifs <- function(seq) {
	c <- 0
	for (motif in motifs2Count) {
		c <- c + 1
		#for (i in (1:2)) {
		#}
		if (motif=="+") {
			motifPresent[c] <- motifPresent[c-1] + motifPresent[c-2]
		} else {
			# Old Style with non-greedy search before 12/11/14
			#match <- gregexpr(motif, seq, perl=TRUE)
			# New Style with greedy search
			patternStr <- paste("(?=", motif, ")", sep="")
			match <- gregexpr(patternStr, seq, perl=TRUE)
			# end New Style
			if (unlist(match)[1]==-1) {
				motifCount[c] <- 0
			} else {
				motifCount[c] <- length(unlist(match))
			}
		}
	}

	return(paste(motifCount, collapse="\t"))
}
