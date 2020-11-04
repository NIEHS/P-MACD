## This code was developed and authored by Les Klimczak, Ph.D. 
## Unauthorized commercial reuse of the code and removal of this notice
## are prohibited.
## This research was supported by the Intramural Research Program of the NIH,
## National Institute of Environmental Health Sciences.

import os, re, sys, copy, numpy, datetime, __builtin__

FILESDIR = sys.argv[1]

execfile(sys.argv[1] + "findMotifs.py")

time40 = datetime.datetime.now()

def which(tf):
	return filter(lambda x: tf[x], range(len(tf)))

#def all(tf):
#	return reduce(lambda x, y: x&y, tf)
	
def duplicated(list):
	seen = []
	duplicated = []
	for el in list:
		if el in seen:
			duplicated.append(True)
		else:
			duplicated.append(False)
		seen.append(el)
	return duplicated

def ifelse(test, yes, no):
	if test:
		return yes
	else:
		return no

OLDEXT = "_anz4.txt"
BATCH = True
BATCH = False

#CLUSTERANZ = True
CLUSTERANZ = False

if CLUSTERANZ:
	OLDEXT = "_cluster.txt"
	ORIGDIR = "/data/PCAWG_12oct_passonly/MAF_2583/histology_split/A3A_A3B/res_ytCa/"
	motifString = re.sub('.*/res', '', ORIGDIR)[:-1]

sumTitles = ["Sample","A_T_coord_clusters", "G_C_coord_clusters", "Non_coord_clusters", "clusters","mutations", "complex","insertions", "deletions", "indels", "substitutions", "bases"]
substTo = {'A':("T","G","C"),'T':("A","C","G"),'G':("C","T","A"),'C':("G","A","T")}

colTitles = []
for title in findTitles:
	if (len(title)==1):
		mutBase = title
	else:
		mutBase = title[which([x<90 for x in [ord(x) for x in title]])[0]]
	for base in substTo[mutBase]:
		if (base==mutBase):
			continue
		colTitle = title + "_to_" + base
		colTitles.append(colTitle)

	colTitles.extend((title, title.lower()))
	if (len(title)>1):
		colTitles.extend((title + "_per_mut", title + "_per_" + mutBase, title.lower() + "_per_" + mutBase.lower(), "enrich_" + title, "freq_" + title, "reliable>=30"))

headers = sumTitles + colTitles
numCols = len(headers)
totals = [0] * numCols
numFixedCols = 32

def newCounter(type, subtype, sumTitle, sumSubtitle, typeColumn, typeValue, subtypeColumn, subtypeValue, subtypeComparison):
	type = type
	subtype = subtype
	sumTitle = sumTitle
	sumSubtitle = sumSubtitle
	typeColumn = typeColumn
	typeValue = typeValue
	if (typeValue[1:2]==","):
		typeValue = typeValue.split(",")
	subtypeColumn = subtypeColumn
	subtypeValue = subtypeValue
	subtypeComparison = subtypeComparison
	if ((subtypeComparison=="equal") or (subtypeComparison=="greater")):
		subtypeValue = int(subtypeValue)
	typeField = [0]
	subtypeField = [0]
	
	# 07/15/16 counters replaced with variants
	variants = [{}]
	complex = [{}]
	clusters = [{}]
	mutations = [0]

	mutList = [{}]
	baseList = [{}]
	baseCountList = [{}]

	totals = [0] * numCols
	output = [""]

	def initOutput(inputPref):
		if (type=='99'):
			typeString = "_" + sumSubtitle
			outputFile = inputPref + "_sum" + typeString + ".txt"
		else:
			typeString = "%02d" % int(type)
			outputFile = inputPref + "_sum" + typeString + subtype + ".txt"
		
		#if (file.exists(outputFile)) next

		output[0] = open(outputFile, "w")
		if (type!='99'):
			output[0].write("#" + sumTitle + " " + sumSubtitle + "\n")
		output[0].write('\t'.join(headers) + "\n")

	def closeOutput():
		output[0].close()

	def initFieldNames(fieldNames):
		typeField[0] = which([x==typeColumn for x in fieldNames])
		if (typeField[0]):
			typeField[0] = typeField[0][0]
		subtypeField[0] = which([x==subtypeColumn for x in fieldNames])
		if (subtypeField[0]):
			subtypeField[0] = subtypeField[0][0]
		#print type, subtype, typeColumn, typeField[0], subtypeField[0]

	def initCounters(mutList0, baseList0, baseCountList0):
		variants[0] = {}
		complex[0] = {}
		clusters[0] = {}
		mutations[0] = 0

		mutList[0] = copy.deepcopy(mutList0)
		baseList[0] = copy.deepcopy(baseList0)
		baseCountList[0] = copy.deepcopy(baseCountList0)

	def count(fields):
		if (sumTitle!="All Mutations"):

			#if (not any([x==typeValue for x in fields])):
			if (__builtin__.type(typeValue) is list):
				if (not any([x==fields[typeField[0]] for x in typeValue])):
					return
			else:
				if (not fields[typeField[0]]==typeValue):
					return
			if subtypeField[0]!=[]:
				if ((subtypeComparison=="equal") or (subtypeComparison=="greater")):
					subtypeFieldValue = int(fields[subtypeField[0]])
				else:
					subtypeFieldValue = fields[subtypeField[0]]
				if (not compare(subtypeFieldValue, subtypeValue, subtypeComparison)):
					return

			#if (True):
			#if (type=="10"):
			#if ((type!="10")&(fields[typeField]!="N")) {
			#	print lnum, type, subtype, subtypeColumn, subtypeField[0]
			# typeField, subtypeField
			#}

		if (not variants[0].has_key(fields[VARIANT_TYPE_FIELD]) and (fields[COMPLEX_ID_FIELD]=="")):
			(variants[0])[fields[VARIANT_TYPE_FIELD]] = 1
		else:
			if (fields[COMPLEX_ID_FIELD]==""):
				(variants[0])[fields[VARIANT_TYPE_FIELD]] = (variants[0])[fields[VARIANT_TYPE_FIELD]] + 1

		if (fields[CLUSTER_ID_FIELD]!=""):
			(clusters[0])[fields[CLUSTER_ID_FIELD]] = fields[CLUSTER_COORD_FIELD]

		if (fields[COMPLEX_ID_FIELD]!=""):
			(complex[0])[fields[COMPLEX_ID_FIELD]] = fields[COMPLEX_ID_FIELD]
			#cat("c", fields[COMPLEX_ID_FIELD], "c")

		# 04/15/14 change to all non-complex rows
		if (fields[COMPLEX_ID_FIELD]==""):
			mutations[0] = mutations[0]+1

		if ((fields[COMPLEX_ID_FIELD]=="") and (fields[VARIANT_TYPE_FIELD]=="SNP")):
			# 04/15/14 change to all non-complex rows
			#mutations <- mutations+1
			# for (title in findTitles[as.logical(as.numeric(fields[fieldNames %in% findTitles]))]) {
			for title in [findTitles[x] for x in which([bool(x) for x in [int(x) for x in [fields[x] for x in which([x in findTitles for x in fieldNames])]]])]:
				mutationCount = (mutList[0])[title][fields[TUMOR_SEQ_ALLELE2_FIELD]]
				(mutList[0])[title][fields[TUMOR_SEQ_ALLELE2_FIELD]] = mutationCount + 1

			#baseCounts <- as.numeric(fields[fieldNames %in% countTitles])
			baseCounts = [int(x) for x in [fields[x] for x in uniqueFieldNumbers]]
			for i in range(len(baseCounts)):
				(baseCountList[0])[uniqueFieldNames[i]] = (baseCountList[0])[uniqueFieldNames[i]] + baseCounts[i]
		else:
			sys.stdout.write("X")

	def writeSampleSum():
		print " writing sample in row", lnum
		complex[0] = len((complex[0]).keys())
		# translated ifelse won't work b/c it's trying to evaluate the False option
		if (variants[0].has_key("INS")):
			insertions = variants[0]["INS"]
		else:
			insertions = 0

		if (variants[0].has_key("DEL")):
			deletions = variants[0]["DEL"]
		else:
			deletions = 0

		if (variants[0].has_key("SNP")):
			substitutions = variants[0]["SNP"]
		else:
			substitutions = 0

		# 04/15/14 change to all rows
		#mutations <- mutations + complex + insertions + deletions
		mutations[0] = mutations[0] + complex[0]
		coordBases = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
		# iterates over keys, not values like R
		for key in clusters[0]:
			coordBases[clusters[0][key]] = coordBases[clusters[0][key]] + 1

		genCounts = (coordBases["A"] + coordBases["T"], coordBases["G"] + coordBases["C"], coordBases["N"], len(clusters[0].keys()), mutations[0], complex[0], insertions, deletions, insertions+deletions, substitutions, substitutions*41)
		#print genCounts
		substCounts = [sampleID]
		substCounts.extend(genCounts)

		for title in findTitles:
			motifCounts = []
			#print(title)
			if (len(title)==1):
				mutBase = title
			else:
				mutBase = title[which([x<90 for x in [ord(x) for x in title]])[0]]

			mutBaseSum = 0
			#for (base in c("A", "T", "G", "C")) {
			for base in substTo[mutBase]:
				if (base==mutBase):
					continue
				baseSubstCount = mutList[0][title][base]
				motifCounts.extend([baseSubstCount])
				mutBaseSum = mutBaseSum + baseSubstCount 

			motifCountName = title.lower() + "_counts"
			totMotifCount = baseCountList[0][motifCountName]
			motifCounts.extend([mutBaseSum, totMotifCount])
			if (len(title)==1):
				baseList[0][mutBase] = mutBaseSum
			else:
				mutMotifperMut = numpy.float64(mutBaseSum)/baseList[0][mutBase]
				totMotifperBase = numpy.float64(totMotifCount)/baseCountList[0][mutBase.lower() + "_counts"]
				reliable = (mutBaseSum >= 30)*1
				#print mutations[0], mutBaseSum, totMotifperBase, totMotifCount
				motifCounts.extend((numpy.float64(mutBaseSum)/mutations[0], mutMotifperMut, totMotifperBase, numpy.float64(mutMotifperMut)/totMotifperBase, numpy.float64(mutBaseSum)/totMotifCount, reliable))

			#print(motifCounts)
			substCounts.extend(motifCounts)

		totals[1:(numCols)] = [x+y for x,y in zip(totals[1:(numCols)], (substCounts[1:(numCols)]))]
		outLine = substCounts[0] + '\t' + '\t'.join([x.__str__() for x in substCounts[1:(numCols)]]) + "\n"
		output[0].write(outLine.replace('nan', 'NaN'))
		#cat(substCounts, "\n")

	fixedCols = 32

	def writeTotals():
		t = 0
 		for title in findTitles:
			if (len(title)==1):
				continue
			else:
				mutBase = title[which([x<90 for x in [ord(x) for x in title]])[0]]

			colOffset = fixedCols + t*11
			titleTotals = [0]*10
			titleTotals[0:5] = totals[(colOffset):(colOffset+5)]
	 		titleTotals[5] = numpy.float64(titleTotals[3])/totals[which([x=="mutations" for x in headers])[0]]
	 		titleTotals[6] = numpy.float64(titleTotals[3])/totals[which([x==mutBase for x in headers])[0]]
	 		titleTotals[7] = numpy.float64(titleTotals[4])/totals[which([x==mutBase.lower() for x in headers])[0]]
			titleTotals[8] = numpy.float64(titleTotals[6])/titleTotals[7]
			titleTotals[9] = numpy.float64(titleTotals[3])/titleTotals[4]
			totals[(colOffset):(colOffset+10)] = titleTotals
			t = t + 1

		totals[0] = "Totals"
		outLine = "\t".join([x.__str__() for x in totals]) + "\n"
		output[0].write(outLine.replace('nan', 'NaN'))
		totals[0:(numCols)] = [0] * numCols

	def getValues():
		return (type, subtype, sumTitle, typeValue)
	
	return {'type':type, 'getValues':getValues, 'initOutput':initOutput, 'closeOutput':closeOutput,  'initFieldNames':initFieldNames, 'initCounters':initCounters, 'count':count, 'mutList':mutList, 'writeSampleSum':writeSampleSum, 'writeTotals':writeTotals}

countersList = []
#countersList[2]['getValues']()

rulesFile = "SummaryRulesIntegr1.txt"
# special 12/13/17 
#rulesFile = "SummaryRules5d.txt"
#rulesTable = read.table(rulesFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")

rInput = open(rulesFile, "r")
line = rInput.readline()

#	switch(type, equal=x==y, greater=x>y, isNonBlank=(x!=""), isAny=TRUE)
def compare(x, y, type):
	return {'equal':x==y, 'greater':x>y, 'isNonBlank':x!="", 'isAny':True}[type]

sumCount = 0
while (True):
	line = rInput.readline().rstrip('\r\n')

	if (not line):
		break
	fields = line.split("\t")

	countersList.append(newCounter(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8]))
	sumCount = sumCount + 1

rInput.close()


files = os.listdir(FILESDIR)
#files = [os.listdir(FILESDIR)[3]]
#file = "2014_Fredriksson_HNSC_27_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz4.txt"

for file in files:
	if (file[-len(OLDEXT):]!=OLDEXT):
		continue

	if (BATCH):
		batchNames = ("BLCA", "BRCA", "HNSC", "LUAD", "LUSC")
		batchSubset = 2
		if (file.split("_")[2] not in batchNames[batchSubset]):
			continue

	print(file)
	#print("\n" + file + "\n")

	inputFile = FILESDIR + file
	inputPref = inputFile[:-len(OLDEXT)]

	if CLUSTERANZ:
		sys.stderr.write("Using clusterDEF file\n")
		#origAnz4 = ORIGDIR + re.sub('_anz2.*', '', file) + '_anz4.txt'
		origAnz4 = ORIGDIR + re.sub('_anz2.*', '', file) + motifString + '_anz4.txt'
		inputPref = inputPref + motifString
		#input = os.popen('cut -f1-77 ' + origAnz4)
		# 68 - last anz4 col; 2 - first clusterDef col; 9 - first+9-1
		#input = os.popen('bash -c "paste <(cut -f1-68 ' + origAnz4 + ')' + ' <(cut -f2-9 ' + inputFile + ')"')
		# for tCw
		#input = os.popen('bash -c "paste <(cut -f1-69 ' + origAnz4 + ')' + ' <(cut -f3-11 ' + inputFile + ')"')
		# for tCa: 62, rtCa: 63
		input = os.popen('bash -c "paste <(cut -f1-63 ' + origAnz4 + ')' + ' <(cut -f3-11 ' + inputFile + ')"')
	else:
		input = open(inputFile, "r")

	line = input.readline()
	firstChar = line[:1]
	while firstChar=="#":
		line = input.readline()
		firstChar = line[:1]

	fieldNames = line[:-1].split("\t")

	VARIANT_TYPE_FIELD = which(map(lambda(x): x=='Variant_Type', fieldNames))[0]
	TUMOR_SEQ_ALLELE2_FIELD = which(map(lambda(x): x=='Tumor_Seq_Allele2', fieldNames))[0]

	COMPLEX_ID_FIELD = which(map(lambda(x): x=='Complex_ID', fieldNames))[0]
	CLUSTER_ID_FIELD = which(map(lambda(x): x=='Dataset_Cluster_ID', fieldNames))[0]
	CLUSTER_COORD_FIELD = which(map(lambda(x): x=='Cluster_Coordination', fieldNames))[0]
	INPUT_SAMPLE_FIELD = which(map(lambda(x): x=='Tumor_Sample_Barcode', fieldNames))[0]


	#uniqueFieldNumbers <- (fieldNames %in% countTitles) & !duplicated(fieldNames)
	#uniqueFieldNames <- fieldNames[uniqueFieldNumbers]
	tf1 = map(lambda x: x in countTitles, fieldNames)
	tf2 = map(lambda x: not x, duplicated(fieldNames))
	uniqueFieldNumbers = which(map(lambda x: tf1[x] and tf2[x], range(len(fieldNames))))
	uniqueFieldNames = map(lambda(x):  fieldNames[x], uniqueFieldNumbers)
	
	mutList = dict(zip(findTitles, [0]*len(findTitles)))
	baseList = {'A': 0, 'T':0, 'G':0, 'C':0}
	for key in mutList.keys():
		mutList[key] = dict(baseList)
	baseCountList = dict(zip(uniqueFieldNames, [0]*len(uniqueFieldNames)))

	print("Initializing...")

	for sumNum in range(sumCount):
			countersList[sumNum]['initOutput'](inputPref)
			countersList[sumNum]['initFieldNames'](fieldNames)
			countersList[sumNum]['initCounters'](mutList,baseList,baseCountList)

	print("Counting...")

	sampleID = ""
	firstSample = True

	lnum = 0
	while (True):
	#while (lnum<1000):
		line = input.readline()

		if (not line):
			for sumNum in range(sumCount):
				countersList[sumNum]['writeSampleSum']()
				countersList[sumNum]['writeTotals']()
				countersList[sumNum]['closeOutput']()
			print "End: not line"
			break

		lnum = lnum+1

		if (lnum%1000==0):
			sys.stdout.write('.')

		fields = line.split("\t")
		
		if (all(map(lambda(x): x=="", fields))):
			for sumNum in range(sumCount):
				countersList[sumNum]['writeSampleSum']()
				countersList[sumNum]['writeTotals']()
				countersList[sumNum]['closeOutput']()
			print "End: empty fields"
			break
		
		if (fields[INPUT_SAMPLE_FIELD]!=sampleID):
			if (not firstSample):
				for sumNum in range(sumCount):
					print("Writing %d" % sumNum)
					countersList[sumNum]['writeSampleSum']()
					countersList[sumNum]['initCounters'](mutList,baseList,baseCountList)
			firstSample = False
			sampleID = fields[INPUT_SAMPLE_FIELD]
			print(sampleID)
		
		for sumNum in range(sumCount):
			countersList[sumNum]['count'](fields)

	input.close()
	
	print
	for sumNum in range(sumCount):
		print countersList[sumNum]['getValues']()


os.remove(sys.argv[1] + "findMotifs.py")
time50 =  datetime.datetime.now()
print(time50 - time40)
	


