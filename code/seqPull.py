## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

# based on generic_seq_pull.pl donated by Sara Grimm, 08/2011

REF_FA  = os.environ['REF_FA']

import re

def getSeq(chr, start, end):
	chrFirstPos = chr2x[chr]
	if (end>chrSizes[chr]):
		return 'ENDOFCHROM'
	lineCt1 = start // lineSegSize
	if (start % lineSegSize == 0):
		lineCt1 = lineCt1-1
	lineCt2 = end // lineSegSize
	if (end % lineSegSize == 0):
		lineCt2 = lineCt2-1
	readPt = chrFirstPos + start + lineCt1 - 1   # first position + query position + line breaks - 1
	readSize = end - start + 1 + (lineCt2-lineCt1)
	faCon.seek(readPt)
	readSeq = faCon.read(readSize)
	readSeq = re.sub('\n', '', readSeq).lower()
	return readSeq

faiFile = REF_FA + ".fai"
faiCon = open(faiFile, "r")

chr2x = dict()
chrSizes = dict()
while (True):
	line = faiCon.readline()
	if (len(line)==0):
		break
	faiFields = line.split("\t")
	chr2x[faiFields[0]] = int(faiFields[2])
	# 01/07/19 chr boundary check
	chrSizes[faiFields[0]] = int(faiFields[1])

lineSegSize = int(faiFields[3])

faiCon.close()

faCon = open(REF_FA, "r")

# getSeq("chr2", 220413922, 220413982)
#'agcaggaagatcgcaaggagaaaaggcggaaaaagagacctcctcgggctcccctcagagg'

#close(faCon)

