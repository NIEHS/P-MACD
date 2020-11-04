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

def countMotifs(seq):
	motifCount = [0] * len(motifs2Count)
	c = 0
	for motif in motifs2Count:
		if (motif=="+"):
			motifPresent[c] = motifPresent[c-1] + motifPresent[c-2]
		else:
			# Old Style with non-greedy search before 12/11/14
			# R: match = gregexpr(motif, seq, perl=TRUE)
			# Python: 
			#match = re.findall(motif, seq) 
			# New Style with greedy search (overlapping counts)
			# R: patternStr = paste("(?=", motif, ")", sep="")
			# R: match = gregexpr(patternStr, seq, perl=TRUE)
			# Python: 
			patternStr = "(?=" + motif + ")"
			match = re.findall(patternStr, seq)
			# end New Style
			motifCount[c] = "%d" % len(match)
		c = c + 1

	return list2tabString(motifCount)

