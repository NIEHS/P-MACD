## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

import os, sys

FILESDIR = sys.argv[1]
SCRIPTDIR = sys.argv[2]
CUSTOMMOTIFS = False
CUSTOMMOTIFS = os.environ['CUSTOMMOTIFS'] == 'TRUE'
CUSTOMTYPE = os.environ['CUSTOMTYPE']
print('Running from package')

def which(tf):
	return filter(lambda x: tf[x], range(len(tf)))
	
def list2tabString(l):
	return reduce(lambda x, y: x+'\t'+y, l)

execfile(SCRIPTDIR + "seqPull.py")
#import seqPull

#PANCANCER = True
PANCANCER = False
if (PANCANCER):
	chr2x["chr23"] = chr2x.pop("chrX")
	chr2x["chr24"] = chr2x.pop("chrY")

if (not CUSTOMMOTIFS):
	execfile(SCRIPTDIR + "findMotifs.py")
execfile(SCRIPTDIR + "setupPy.py")
execfile(FILESDIR + "findMotifs.py")

#import findMotifs

if (not CUSTOMMOTIFS):
	blankFindMotifs = list2tabString([""] * (len(motifs2Find)+len(apobecTitles)))
else:
	blankFindMotifs = list2tabString([""] * len(motifs2Find))

# why commented out ?
#execfile(SCRIPTDIR + "countMotifs.py")
# SPECIAL commented out for CGI 
blankCountMotifs = list2tabString([""] * len(motifs2Count))
#blankCountMotifs = ''

BEFORELENGTH = 20
AFTERLENGTH = 20
seq2read = BEFORELENGTH + AFTERLENGTH + 2
contextTitle = '\"CONTEXT(+/-' + ("%d" % BEFORELENGTH) + ')\"'
OLDEXT = "_anz2.txt"
NEWEXT = "_anz4.txt"

# 06/21/19 boundary check output
boundaryChecks = []

if (CUSTOMMOTIFS):
	NEWEXT = "_" + CUSTOMTYPE + "_anz4.txt"

def capsubstr(x, start, stop):
	x = x[:start-1] + x[start-1:stop].upper() + x[stop:]
	return x

UCSC_COORD = True
UCSCCORR = 0

INPUT_CHR_FIELD = 0 
INPUT_START_FIELD = 0
INPUT_END_FIELD = 0
VARIANT_TYPE_FIELD = 0
TUMOR_SEQ_ALLELE2_FIELD = 0
COMPLEX_ID_FIELD = 0

INSERTAFTERCOLNAME = "Tumor_Seq_Allele2"

def getFieldPosFromName():
	global INPUT_CHR_FIELD, INPUT_START_FIELD, INPUT_END_FIELD, VARIANT_TYPE_FIELD, TUMOR_SEQ_ALLELE2_FIELD, COMPLEX_ID_FIELD, INSERTAFTERCOLNAME
	
	INPUT_CHR_FIELD = which([x=="Chromosome" for x in fieldNames])[0]
	INPUT_START_FIELD = which([x=="Start_Position" for x in fieldNames])
	if (len(INPUT_START_FIELD)==0):
		# non-standard name in many Broad MAFs (violates NCI MAF definition)
		INPUT_START_FIELD = which([x=="Start_position" for x in fieldNames])[0]
	else:
		INPUT_START_FIELD = INPUT_START_FIELD[0]

	INPUT_END_FIELD = INPUT_START_FIELD
	#INPUT_END_FIELD = which(fieldNames=="End_position")
	VARIANT_TYPE_FIELD = which([x=="Variant_Type" for x in fieldNames])[0]
	TUMOR_SEQ_ALLELE2_FIELD = which([x=="Tumor_Seq_Allele2" for x in fieldNames])[0]
	COMPLEX_ID_FIELD = which([x=="Complex_ID" for x in fieldNames])[0]

	INSERTAFTERCOLNAME = "Tumor_Seq_Allele2"


files = os.listdir(FILESDIR)
files.sort()
#file = files[0]
#[x[-8:]=='anz2.txt' for x in os.listdir(FILESDIR)]

for file in files:
	if (file[-len(OLDEXT):]!=OLDEXT):
		continue
	#if (file[0:4]!="CESC"):
	#if (file[17:21]!="HNSC"):
	#	continue
	print(file)

	inputFile = FILESDIR + file
	inputPref = inputFile[:-len(OLDEXT)]

	outputFile = inputPref + NEWEXT
	
	# seq file names do not contain motif types, they are motif-independent
	# inputPref doesn't contain motif - motif is part of NEWEXT
	seqFile = inputPref + ".seq"
	
	output = open(outputFile, "w")
	
	print(inputPref, seqFile)
	if (os.path.exists(seqFile)):
		print('Using existing sequence file...')
		seqFromFile = True
		seq = open(seqFile, "r")
	else:
		seqFromFile = False
		print('Creating a new sequence file...')
		seq = open(seqFile, "w")

	if (not CUSTOMMOTIFS):
		newTitles = list2tabString(findTitles + apobecTitles + (contextTitle,) + countTitles)
	else:
		newTitles = list2tabString(findTitles + (contextTitle,) + countTitles)

	input = open(inputFile, "r")
	
	line = input.readline()[:-1]
	headPos = input.tell()
	fieldNames = line.split("\t")
	#print fieldNames

	if (PANCANCER):
			INPUT_CHR_FIELD = which(fieldNames=="chr")
			INPUT_START_FIELD = which(fieldNames=="pos")
			#INPUT_START_FIELD = which(fieldNames=="Start_Position")
			INPUT_END_FIELD = INPUT_START_FIELD
			#INPUT_END_FIELD = which(fieldNames=="End_position")
			VARIANT_TYPE_FIELD = which(fieldNames=="classification")
			TUMOR_SEQ_ALLELE2_FIELD = which(fieldNames=="newbase")
			COMPLEX_ID_FIELD = which(fieldNames=="Complex_ID")

			INSERTAFTERCOLNAME = "newbase"
		
	else:
		getFieldPosFromName()
	
	line = input.readline()[:-1]
	fields = line.split("\t")
	if (fields[INPUT_CHR_FIELD][:3]=='chr'):
		CHRPREF = ''
	else:
		CHRPREF = 'chr'
	print(CHRPREF)
	input.seek(headPos)
	
	insertAfterColNum = which([x==INSERTAFTERCOLNAME for x in fieldNames])[0]
	if (insertAfterColNum==len(fieldNames)):
		newFieldNames = fieldNames[:insertAfterColNum+1] + [newTitles,]
	else:
		newFieldNames = fieldNames[:insertAfterColNum+1] + [newTitles,] + fieldNames[(insertAfterColNum+1):]
	
	output.write(list2tabString(newFieldNames)+'\n')


	c = 0
	while True:
	#while (c<2000):
		line = input.readline()[:-1]
		if (not line):
			break

		lastChar = line[-1:]
		c = c + 1

		fields = line.split("\t")
		# not needed in Python
		#if (lastChar=="\t"):
		#	fields.append("")
		if (c % 1000==0):
			sys.stdout.write('.')

		# don't check for big files with known chromosome format
		# uncomment for original MAF, comment for processed anz
		# ClusterFinder stripped chr before modification on 01/28/19
		#fields[INPUT_CHR_FIELD] = "chr" + fields[INPUT_CHR_FIELD]
		# auto-discover 11/22/19 
		fields[INPUT_CHR_FIELD] = CHRPREF + fields[INPUT_CHR_FIELD]
		if (not chr2x.has_key(fields[INPUT_CHR_FIELD])):
			#writeLines("X", con=output)
			continue

		inputStart = int(fields[INPUT_START_FIELD])
		inputEnd = int(fields[INPUT_END_FIELD])
		overOneNt = inputEnd - inputStart
		seqStart = inputStart+UCSCCORR - BEFORELENGTH
		seqEnd = inputEnd + AFTERLENGTH
		
		# 01/07/19 chr boundary check
		if (seqStart<1):
			# 06/21/19 boundary check output
			print '\n\n\n\nBOUNDARY CHECK FAILED in line ' + ("%d" % c) + ' \n\n\n\n'
			boundaryChecks.append(c)
			if (not seqFromFile):
				seq.write('\n')
			else:
				readSeq = seq.readline()[:-1]
			continue

		if (not seqFromFile):
			readSeq = getSeq(fields[INPUT_CHR_FIELD], seqStart, seqEnd)
			if readSeq=='ENDOFCHROM':
				# 06/21/19 boundary check output
				print '\n\n\n\nBOUNDARY CHECK FAILED in line ' + ("%d" % c) + ' \n\n\n\n'
				boundaryChecks.append(c)
				seq.write('\n')
				continue
		else:
			#readSeq = seq.readline()[:-1]
			# 06/21/19 boundary check output
			#readSeq = seq.read(seq2read)[:-1]
			readSeq = seq.readline()[:-1]
			if (len(readSeq)==0):
				print '\n\n\n\nBOUNDARY CHECK FAILED in line ' + ("%d" % c) + ' \n\n\n\n'
				boundaryChecks.append(c)
				continue

		capReadSeq = capsubstr(readSeq, BEFORELENGTH+1, BEFORELENGTH+1+overOneNt)
		
		if (fields[VARIANT_TYPE_FIELD]=="INS"):
					#capReadSeq = paste(substr(readSeq, 1, BEFORELENGTH+1), toupper(fields[TUMOR_SEQ_ALLELE2_FIELD]), substr(readSeq, BEFORELENGTH+2,BEFORELENGTH+2+AFTERLENGTH), sep="")
					#capReadSeq = paste(substr(readSeq, 1+1, BEFORELENGTH+1), "INS", substr(readSeq, BEFORELENGTH+2,BEFORELENGTH+2-1+AFTERLENGTH), sep="")
					capReadSeq = readSeq[1:BEFORELENGTH+1] + "INS" + readSeq[(BEFORELENGTH+1):(BEFORELENGTH+2-1+AFTERLENGTH)]
					print(readSeq[:-1] + " - " + capReadSeq)

		if (not seqFromFile):
			seq.write(readSeq+'\n')
			
		#print "X"+fields[COMPLEX_ID_FIELD]+"X"
		
		if ((fields[VARIANT_TYPE_FIELD]=="SNP") and (fields[COMPLEX_ID_FIELD]=="")):
			findString = findMotifs(capReadSeq)
			if (not CUSTOMMOTIFS):
				apobecString = isApobec(findString)
				findString = findString + "\t" + apobecString

			# SPECIAL commented out for CGI 
			countString = countMotifs(readSeq)
			#countString = ''

		else:
			findString = blankFindMotifs
			countString = blankCountMotifs

		# combine find and count strings with sequence			
		motifString = findString + "\t" + capReadSeq + "\t" + countString

		if (insertAfterColNum==len(fieldNames)):
			newFields = fields[:insertAfterColNum+1] + [motifString,]
		else:
			newFields = fields[:insertAfterColNum+1] + [motifString,] + fields[(insertAfterColNum+1):]

		output.write(list2tabString(newFields)+'\n')
	
	input.close()
	output.close()
	seq.close()
	
	# 06/21/19 boundary check output
	if len(boundaryChecks)>0:
		bdck = open(inputPref + '.bdck', "w")
		for bd in boundaryChecks:
			bdck.write("%d\n" % bd)
		bdck.close()

faCon.close()

