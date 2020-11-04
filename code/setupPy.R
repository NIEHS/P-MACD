## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

Sys.setenv(CUSTOMMOTIFS=CUSTOMMOTIFS, CUSTOMTYPE=CUSTOMTYPE, REF_FA=REF_FA)
if (CUSTOMMOTIFS) {
	pythonConfig <- paste("customConfig_", CUSTOMTYPE, ".py", sep="")
	if (file.exists(pythonConfig)) {
		file.copy(pythonConfig, paste(FILESDIR, "findMotifs.py", sep=""), overwrite=TRUE)
	} else {
		configFile <- file(pythonConfig, "w")
		configLine <- paste('motifs2Find = ("', paste(motifs2Find, collapse='", "'), '")', sep="")
		writeLines(configLine, configFile)
		configLine <- paste('findTitles = ("', paste(findTitles, collapse='", "'), '")', sep="")
		writeLines(configLine, configFile)		
		configLine <- paste('motifs2Count = ("', paste(motifs2Count, collapse='", "'), '")', sep="")
		writeLines(configLine, configFile)
		configLine <- paste('countTitles = ("', paste(countTitles, collapse='", "'), '",)', sep="")
		cat(configLine, "\n")
		writeLines(configLine, configFile)	
		# why commented out ?	
		#configLine <- 'countTitles = map(lambda x: x+"_counts", countTitles)'
		configLine <- 'countTitles = tuple([x+"_counts" for x in countTitles])'
		writeLines(configLine, configFile)		
			
		close(configFile)
		file.copy(pythonConfig, paste(FILESDIR, "findMotifs.py", sep=""), overwrite=TRUE)
	}
} else {
	file.copy("findMotifs_std.py", paste(FILESDIR, "findMotifs.py", sep=""), overwrite=TRUE)
	#file.copy("findMotifs.py", paste(FILESDIR, "findMotifs.py", sep=""), overwrite=TRUE)
}
