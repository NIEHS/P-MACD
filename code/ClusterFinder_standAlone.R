## This code was developed and authored by Deepak Mav, Shawn Harris and Arpit Tandon. Unauthorized commercial reuse of the code and removal of this notice are prohibited. This research was supported by the Intramural Research Program of the NIH, National Institute of Environmental Health Sciences.
## Update: Can skip lines with "#" before the header
## Keep the non-standard chromosome name as it is
options(warn=-1)
library('data.table')

## trimming function for proper naming convetnion
trim <- function (s, recode.factor = TRUE, ...) {
    s <- sub(pattern = "^ +", replacement = "", x = s)
    s <- sub(pattern = " +$", replacement = "", x = s)
    s
}

##loads the MAF file in data.table data structure and performs filtering and other steps. It prepares required objects for cluster finding
Load.data <- function(mutFile, headerRowIdxMut, headerRowIdxChr, chromFile, genome.pos,chrom.pos,var.type, ref.pos, complex.threshold, ChromName.pos, ChromSize.pos, filtering, conservativeApproach,filterInfo ) {
    headerRowIdxMut <- headerRowIdxMut - 1 
    headerRowIdxChr <- headerRowIdxChr - 1
    if((complex.threshold>100)|(complex.threshold<0)) stop("complex.threshold must be in [0,100]");
	
	##Reading file before hand and finding any lines starting with "#" before header and then skip that line
	top 	<- "^#"
	con  <- file(mutFile, open = "r")
	topLine = readLines(con, n = 1, warn = FALSE)
	if(grepl(top, topLine)){
		headerRowIdxMut = headerRowIdxMut + 1
		while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
			if(grepl(top, oneLine)){
	        		headerRowIdxMut= headerRowIdxMut + 1
      			}
	      		else{
	        		break
      			}
    		}
  	} 
    
  close(con)
    indata  <- read.csv(mutFile, sep="\t", skip =headerRowIdxMut, header = TRUE, comment.char = "",check.names = FALSE, quote="" )
    time1   <- proc.time()
    print(paste("Processing file ", basename(mutFile)))
    
    ##removing any duplicate mutations.
		cat(genome.pos,chrom.pos,ref.pos,"\n")
    indataUnique <- indata[!duplicated(indata[c(genome.pos,chrom.pos,ref.pos)]),]
    indata <- indataUnique    
    
    nchroms  <- length(chroms<-levels(indata[,chrom.pos]));
    ngenomes <- length(genomes<-levels(indata[,genome.pos]));
	originalMutations <- vector("list",length(genomes))
    names(originalMutations) <- genomes
	
    for (i in 1:ngenomes) {
      originalMutations[genomes[i]] <- length(which(indata[,genome.pos] == genomes[i]))
    } 
    if(filtering) {
        infos <- unlist(strsplit(filterInfo,split="]"))
        for (i in 1:length(infos)) {
            infos[i]    <- gsub("\\[","",infos[i])
            splitInfos  <- unlist(strsplit(infos[i],split="="))
            colName     <- as.character(splitInfos[1])
            values      <- unlist(strsplit(splitInfos[2],split="\\|"))
            valBuff     <- c()
            if(!as.logical(match(infos[i], "(None)=", nomatch=FALSE))) {
				indata[colName][is.na(indata[colName])]= ""
                for(j in 1:length(values)) {
                    if(is.na(values[j])) {
                            values[j] = ""
                    }
                    valBuff <- c(valBuff, as.character(values[j]))
                }
                tempIndata  <- subset(indata,(get(colName) %in% valBuff))
                indata      <- tempIndata
            }
        }
    }
	
	##The available genomes, after filtering
    genomesAvai <- unique(indata[,genome.pos])
    genomesRest <- genomes[!(genomes %in% genomesAvai)] 

	if(nrow(indata) == 0) {
        stop("      No data is available after filtering.
        Please check the filtering criteria and the mutation file.")
        warning("No data available")
    }
    chromosomes                 <- read.csv(chromFile, sep="\t",skip =headerRowIdxChr, header = TRUE, comment.char = "",check.names = FALSE) 
    #chromosomes[,ChromName.pos] <- trim(gsub("Chr|chr","",(sapply(strsplit(as.character(chromosomes[,ChromName.pos]),split="_"),function(x) x[1]))));
    chrom.class                 <- class(indata[,chrom.pos]);
    genome.class                <- class(indata[,genome.pos]);
    indata[,genome.pos]         <- as.factor(indata[,genome.pos]);
    indata[,chrom.pos]          <- as.factor(indata[,chrom.pos]);
    #levels(indata[,chrom.pos])  <- gsub("Chr|chr","",sapply(strsplit(levels(indata[,chrom.pos]),split="_"),function(x) x[1]))
	chrOrder 					<- intersect(c((1:100),"X","Y","M"),levels(indata[,chrom.pos]))
	chrOrder 					<- c(chrOrder,sort(levels(indata[,chrom.pos])[!levels(indata[,chrom.pos])%in% chrOrder]))
	indata[,chrom.pos]			<- factor(indata[,chrom.pos], chrOrder,ordered=TRUE)
	indata[,genome.pos]         <- as.factor(indata[,genome.pos]);
    indata[,chrom.pos]          <- as.factor(indata[,chrom.pos]);
		# 09/22/16 LJK
		# sorting MAF not necessary if input already sorted
    indata                      <- indata[do.call(order,indata[,c(genome.pos, chrom.pos, ref.pos)]),];
    genome.size                 <- sum(as.numeric(chromosomes[,ChromSize.pos]));
    ngenomes                    <- length(genomes<-levels(indata[,genome.pos]));
    nchroms                     <- length(chroms<-levels(indata[,chrom.pos]));    
    nMutations                  <- nrow(indata);
    
    ##Proceed only on those chromosomes available in chromosome file. Ex- most likely chrM is not present in chromosome file
    availableChrs               <- intersect(trim(unique(indata[,chrom.pos])),chromosomes[,ChromName.pos])
    navailableChrs              <- length(availableChrs)
    DistBetween.Mutations 		<- Dataset.Mutation.ID <- RT <- Complex.ID <- Complex.Size  <- rep("",nMutations);
    Mutation.ID 				<- Complex.Start <- Complex.End <- rep(0,nMutations); 
	Adjusted_Complex_number     <- vector("list",length(genomes))
    names(Adjusted_Complex_number)<- genomes 
	
    tryCatch(for(j in 1:navailableChrs) {
        size    <- chromosomes[which(chromosomes[,ChromName.pos]==trim(availableChrs[j])),ChromSize.pos];
        w       <- which(indata[,chrom.pos]==availableChrs[j])
        RT[w]   <- size - indata[w,ref.pos];
    },error = function(e){print("missing chromosomes")
        stop("      missing chromosomes: Some chromosomes in mutation file are missing from chromosome file.
        Please ensure that you have all the chromosome information.")
        warning("missing chromosomes")})

    g2 <- 0;
    s1 <- 0;
    for(i in 1:ngenomes) {
        Adjusted_Complex_number_temp <- 0 
        g1 <- 0; 
        w1 <- which(indata[,genome.pos]==genomes[i]);
        for(j in 1:nchroms) {
            w2 <- w1[which(indata[w1,chrom.pos]==chroms[j])];
            n <- length(w2)
            if(n==1) {
                  Complex.Start[w2] <- Complex.End[w2] <- 1;
            }
            else if(n>1) {
                DistBetween.Mutations[w2[2:length(w2)]] <- indata[w2[2:length(w2)],ref.pos] - indata[w2[1:length(w2)-1],ref.pos]
                distance <- diff(indata[w2 ,ref.pos]);
                S <- 1; E <- c(); 
                for(k in 2:n) {
                    if(distance[k-1] > complex.threshold) { 
                        E <- c(E,k-1); 
                        S <- c(S,k)
                    }
                }
                E <- c(E,n);
                for(k in 1:length(S)) {
                    Complex.Start[w2[S[k]]] <- 1;
                    Complex.End[w2[E[k]]] <- 1;
                    
                    if(S[k]<E[k]) {
                        s1 <- s1+1;
                        Complex.ID[w2[S[k]]:w2[E[k]]]   <- s1;
                        Complex.Size[w2[S[k]]:w2[E[k]]] <- E[k]-S[k]+1;
                        Adjusted_Complex_number_temp    <- Adjusted_Complex_number_temp + E[k]-S[k]
                    }  
                }
            }                 
            Mutation.ID[w2] <- g1+(1:n);
            g1 <- g1+n;  
        }
        m <- length(w1)
        Dataset.Mutation.ID[w1] <- g2+(1:m);
        g2 <- g2+m;
        Adjusted_Complex_number[i] <- Adjusted_Complex_number_temp
    }

    indata$DistBetween_Mutations 	<- DistBetween.Mutations;
    indata$"Distance_to_LT_end" 	<- indata[,ref.pos];
    indata$"Distance_to_RT_end" 	<- RT;
    indata$Strain_Mutation_ID 		<- Mutation.ID; indata$Dataset_Mutation_ID <- Dataset.Mutation.ID;  indata$Complex_ID <- Complex.ID; indata$Complex_Size <- Complex.Size; indata$Complex.Start <- Complex.Start; indata$Complex.End <- Complex.End;
    invisible(list(indata=indata,genome.class=genome.class,originalMutations=originalMutations, chrom.class=chrom.class, genome.size =genome.size, genomes=genomes, ngenomes=ngenomes, chroms=chroms, nchroms=nchroms, nMutations=nMutations, genome.pos=genome.pos,chrom.pos=chrom.pos, ref.pos=ref.pos, Adjusted_Complex_number=Adjusted_Complex_number, genomesAvai=genomesAvai, genomesRest=genomesRest));
}

##Searches the clusters
Create.Clusters <- function(objects, dist.threshold, clusterFile, ref.Allele, var.Type, minClusterSize, pValue, conservativeApproach,preFilteredMul=1, fileIndex=0, ComplexAdjustment=ComplexAdjustment) {
    temp	<-data.table(objects$indata)
    rep(as.integer(NA),objects$nMutations) 	-> temp$StrainCluster_ID -> temp$Dataset_Cluster_ID   -> temp$Distance_Between_Clusters -> temp$Cluster_Start 
    rep(as.integer(NA),objects$nMutations) 	-> temp$Cluster_End -> temp$Cluster_Size_Mutations  -> temp$Cluster_Size_Complexes -> temp$Cluster_Length  
    rep("",objects$nMutations) 				-> temp$Cluster_Coordination -> temp$Content_of_non_coordinated_cluster  -> temp$Cluster_InternalCoordination -> temp$Non_Coordinated_Internalcluster;
    rep(as.double(NA),objects$nMutations) 	-> temp$Cluster_Pvalue
    
    s1 = 0;
    timeT <- proc.time()
    for(i in 1:objects$ngenomes) {
        nMutations <- (as.numeric(objects$originalMutations[objects$genomes[i]]) - as.numeric(ComplexAdjustment[objects$genomes[i]]))/preFilteredMul
        prob.mutation <- nMutations/objects$genome.size;
        time2 <- proc.time()
        print(paste("Processing file ", fileIndex, "  genome ", objects$genomes[i]))
        strain <- 0;
		if(objects$genomes[i] %in% objects$genomesAvai){
			w1 <- which(temp[,objects$genome.pos,with=FALSE]==objects$genomes[i]);
			for(j in 1:objects$nchroms) {
				print(paste("Processing file ", fileIndex, "  genome ", objects$genomes[i], "  chromosome ", objects$chroms[j]))
				w.start <- w1[which((temp[w1,objects$chrom.pos,with=FALSE]==objects$chroms[j]) & (temp[w1,Complex.Start]==1))];
				w.end   <-  w1[which((temp[w1,objects$chrom.pos,with=FALSE]==objects$chroms[j]) & (temp[w1,Complex.End]==1))];        
				n <- length(w.start)
				if(n>1) {
					distance <- temp[w.start[2:n],objects$ref.pos,with=FALSE]-temp[w.end[1:(n-1)],objects$ref.pos,with=FALSE];        
					S <- 1; E <- c(); 
					for(k in 2:n) {
						if(distance[k-1]>dist.threshold) { 
							E <- c(E,k-1); 
							S <- c(S,k)
						}
					}
					E <- c(E,n);
					mutArrayE <- mutArrayS <- refAllele <- clusterLoca <- clusterStartLoca <- c();
					alreadyFoundOneCluster <- FALSE
					for(k in 1:length(S)) {
						if(S[k]<E[k]) {
							pVal = pnbinom(as.numeric(temp[w.end[E[k]],objects$ref.pos,with=FALSE])-as.numeric(temp[w.start[S[k]],objects$ref.pos,with=FALSE])-(E[k]-S[k]),E[k]-S[k],prob.mutation)
							if(E[k]-S[k]+1 >= minClusterSize && pVal <= pValue ) {
								if(alreadyFoundOneCluster) {
									clusterDistance = temp[w.start[[S[k]]],objects$ref.pos,with=FALSE] - temp[w.end[mutArrayE[1]],objects$ref.pos,with=FALSE]
									set(temp,w.start[S[k]]:w.end[E[k]],'Distance_Between_Clusters',clusterDistance)
								}
								mutArrayE <- c(E[k])
								refAllele <- c(temp[w.start[S[k]]:w.end[E[k]],ref.Allele,with=FALSE])
								varType   <- c(temp[w.start[S[k]]:w.end[E[k]],var.Type,with=FALSE])
								varTypeUni<- unique(unlist(varType))
								clusterLoca <- c(w.start[S[k]]:w.end[E[k]])
								lenB <- length(refAllele)
								refAlleleTemp <- rev(sort(unique(unlist(refAllele))))
								lenA <- length(refAlleleTemp)
								buffer <- refAlleleTemp[1]
								if(lenA > 1) {
									for( jj in 2:lenA) {
										buffer <- paste(buffer, refAlleleTemp[jj], sep="_")
									}
								}
								##should check if any complex exists in this cluster,then
								##assign 'N' to Cluster_Coordinaton and Cluster_InternalCoordination
								complexExist = FALSE
								for (ii in 1:length(clusterLoca)) {
									if(temp[clusterLoca[ii],Complex_Size] != "") {
										complexExist = TRUE
									}
								}                                
								for (ii in 1:length(clusterLoca)) {
									if(complexExist) {
										set(temp,as.integer(clusterLoca[ii]),'Cluster_Coordination',as.character('N'))
										set(temp,as.integer(clusterLoca[ii]),'Content_of_non_coordinated_cluster', as.character(buffer))
										set(temp,as.integer(clusterLoca[1]),'Cluster_InternalCoordination',as.character('N'))
										set(temp,as.integer(clusterLoca[1]),'Non_Coordinated_Internalcluster',as.character(buffer))
									}
									else if (lenA != 1 || length(varTypeUni) != 1 || !any(varTypeUni == "SNP" )) {
										set(temp,as.integer(clusterLoca[ii]),'Cluster_Coordination',as.character('N'))
										set(temp,as.integer(clusterLoca[ii]),'Content_of_non_coordinated_cluster', as.character(buffer))
										set(temp,as.integer(clusterLoca[1]),'Cluster_InternalCoordination',as.character('N'))
										set(temp,as.integer(clusterLoca[1]),'Non_Coordinated_Internalcluster',as.character(buffer))
									}
									else {
										set(temp,as.integer(clusterLoca[ii]),'Cluster_Coordination',as.character(refAlleleTemp))
										set(temp,as.integer(clusterLoca[1]),'Cluster_InternalCoordination', as.character(refAlleleTemp))
									}
								}
								s1 <- s1+1;
								strain <- strain+1;
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Dataset_Cluster_ID',as.integer(s1))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_Size_Complexes',as.integer(E[k]-S[k]+1))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_Size_Mutations',as.integer(length(w.end[E[k]]:w.start[S[k]])))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_Length',as.integer(temp[w.end[E[k]],objects$ref.pos,with=FALSE]-temp[w.start[S[k]],objects$ref.pos,with=FALSE]+1))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_Pvalue',as.double(pnbinom(as.numeric(temp[w.end[E[k]],objects$ref.pos,with=FALSE])-as.numeric(temp[w.start[S[k]],objects$ref.pos,with=FALSE])-(E[k]-S[k]),E[k]-S[k],prob.mutation)))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_Start',as.integer(0))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'Cluster_End',as.integer(0))
								set(temp,as.integer(w.start[S[k]]),'Cluster_Start',as.integer(1))
								set(temp,as.integer(w.start[S[k]]),'Cluster_End',as.integer(1))
								set(temp,as.integer(w.start[S[k]]:w.end[E[k]]),'StrainCluster_ID',as.integer(strain))
								alreadyFoundOneCluster <- TRUE
							}
						}
					} 
					if(length(mutArrayS) >= 2)        {
						clusterDistance = temp[w.start[mutArrayS[length(mutArrayS)]],objects$ref.pos,with=FALSE] - temp[w.end[mutArrayE[1]],objects$ref.pos,with=FALSE]
						set(temp,as.integer(w.start[temp[1]]:w.end[temp[2]]),'Distance_Between_Clusters',as.integer(clusterDistance))
					}
				}        
			}
			time3 <- proc.time() - time2
			print(paste("Time to process file ", fileIndex, " genome ", objects$genomes[i], " = ", time3["elapsed"]))
		}
	}
    timeC <- proc.time() - timeT
    print(paste("Time to process file ", fileIndex, " = ", timeC["elapsed"]))
    if(objects$genome.class!="factor") {
        #tmp <- as.character(temp[,objects$genome.pos,with=FALSE]);
        #class(tmp) <- "double";
        #temp[,objects$genome.pos] <- tmp;
        #remove(tmp);
    }   
    if(objects$chrom.class!="factor") {
		#tmp <- as.character(temp[,objects$chrom.pos]);
        #class(tmp) <- "double";
        #temp[,objects$chrom.pos] <- tmp;
        #remove(tmp);
    }

	##replacing NAs with "" in newly introduced columns
	newColumns = c("Cluster_Coordination","Content_of_non_coordinated_cluster","Cluster_InternalCoordination","Dataset_Cluster_ID","Cluster_Size_Complexes","Cluster_Length","Cluster_Pvalue","Cluster_Start","Cluster_End","StrainCluster_ID","Distance_Between_Clusters","Cluster_Size_Mutations")
    newColumnsNum = c("Dataset_Cluster_ID","Cluster_Size_Complexes","Cluster_Length","Cluster_Pvalue","Cluster_Start","Cluster_End","StrainCluster_ID","Distance_Between_Clusters","Cluster_Size_Mutations")
	genomeIndata <- as.data.frame(temp) 
	genomeIndata[,newColumns] <- apply(genomeIndata[,newColumns] ,2, function(x){replace (x, is.na(x) ,"" )})
	#genomeIndata[,newColumns] <- apply(genomeIndata[,newColumns] ,2, function(x){ gsub("^\\s+|\\s+$", "", x)})
	genomeIndata[,newColumnsNum] <- sapply(genomeIndata[,newColumnsNum], as.numeric)
	
	genomeIndata[] <- lapply(genomeIndata, as.character)
	write.table(genomeIndata[,setdiff(names(genomeIndata),c("Original.nMutations","Cluster_InternalCoordination","Non_Coordinated_Internalcluster","Complex.Start","Complex.End","Cluster_Start","Cluster_End"))], file=clusterFile,row.names=FALSE,quote=FALSE,sep="\t",na="")
	Sys.chmod(clusterFile, (file.info(clusterFile)$mode | "775"))
    invisible(list(genomeIndata=genomeIndata,Adjusted_Complex_number=objects$Adjusted_Complex_number,genome.pos=objects$genome.pos, genomesRest=objects$genomesRest));        
}

##Creates Summary file
create.SummaryFile <- function(summaryFile = summaryFile,clusterFile=clusterFile, objects,complex.threshold=complex.threshold, minClusterSize=minClusterSize, dist.threshold=dist.threshold, pValue=pValue, conservativeApproach =conservativeApproach, preFilteredMul=preFilteredMul, filterInfo =filterInfo) {
    #If summary file exists, then just append the summary on it
    #else, create the new one
	print(summaryFile)
    if(!file.exists(summaryFile)) {
        buffer                        <- c("Description","Complex_Threshold","Min _#_mutations_in_a_cluster","Distance_threshold","Maximal_p-value","Calculation_Method", "Prefiltered_Multiplier","Filter_Criteria", "Output_File_Name")
        buffer                        <- c("Description","Complex_Threshold","Min _#_mutations_in_a_cluster","Distance_threshold","Maximal_p-value","Calculation_Method", "Prefiltered_Multiplier","Filter_Criteria", "Output_File_Name")
        buffer                        <- c(buffer,"Detected_number_of_clusters","#_of_A-coordinated_cluster","#_of_A_N-coordinated","#_of_C-coordinated_cluster","#_of_C_N-coordinated","#_of_G-coordinated_cluster","#_of_G_N-coordinated")
        buffer                        <- c(buffer,"#_of_T-coordinated_cluster","#_of_T_N-coordinated","#_of_Non_uniform_cluster","Min_cluster_length_(nt)_detected","Max_cluster_length_(nt)_detected")
        buffer                        <- c(buffer,"Minimal_#_of_mutation_complexes_in_a_cluster_detected","Maximal_#_of_mutation_complexes_in_a_cluster_detected")
        buffer                        <- c(buffer,"Minimal_#_of_mutations_in_a_cluster_detected","Maximal_#_of_mutations_in_a_cluster_detected")
        buffer                        <- c(buffer,"Mimimum_p-value_detected","Maximum_p-value_detected","Complex_adjustment","Total_number_of_mutations","Total_number_of_mutation_complexes")
        bufferTemp                    <- ""
        bufferTemp                    <- (paste(unlist(buffer),collapse="\t"))                                                        
        lapply(bufferTemp, write, summaryFile, ncolumns = length(buffer), append=TRUE, sep="")
    }
    filterInfo <- gsub("\\[\\(None\\)=\\]", "" ,filterInfo)
    filterInfo <- gsub("=\\|","='(Blanks)'|",filterInfo)
    filterInfo <- gsub("\\|]","|'(Blanks)']",filterInfo)
    filterInfo <- gsub("\\|\\|","|'(Blanks)'|",filterInfo)
    filterInfo <- gsub("=]","='(Blanks)']",filterInfo)
    
    countA <- countG <- countT <- countC <- countN <- 0
	countA_non <- countG_non <- countT_non <- countC_non <- countN_non <- 0
    ##Count occurrence of A,C,G,T and N coordinate clusters
    print(objects$genome.pos)
    genomeIndata  <- split(objects$genomeIndata, objects$genomeIndata[,objects$genome.pos])
    jj = 0
    totalClusters <- totalA <- totalG <- totalT <- totalC <- totalN <- totalComplexAdj <- totalnoMutations <- totalnoMutCom <- 0
	totalA_non <- totalG_non <- totalT_non <- totalC_non <- totalN_non <- 0
    totalClusterLen <- totalNumClusters <- totalNumComplexes <- totalmaxpValue <- totalminpValue <- c()

    approach = "LiberalApproach"
    if(conservativeApproach) {
        approach = "conservativeApproach"
    }
    for(i in names(genomeIndata)) {
        jj = jj + 1
        buffer                     <- c()
        splitGenomes               <- as.data.frame(genomeIndata[[i]])
        coordinates                <- table(splitGenomes$Cluster_InternalCoordination)
        non_coordinates            <- table(splitGenomes$Non_Coordinated_Internalcluster)
        coordinatesNam             <- names(coordinates)
        non_coordinatesNam         <- names(non_coordinates)
        coordinatesVal             <- as.vector(coordinates)
        non_coordinatesVal         <- as.vector(non_coordinates)
		
		##for coordinates
        countA                     <- coordinatesVal[which(coordinatesNam == "A")][1]
        countG                     <- coordinatesVal[which(coordinatesNam == "G")][1]
        countT                     <- coordinatesVal[which(coordinatesNam == "T")][1]
        countC                     <- coordinatesVal[which(coordinatesNam == "C")][1]
        countN                     <- coordinatesVal[which(coordinatesNam == "N")][1]
		
        if(is.na(countA)) { countA <- 0}
        if(is.na(countG)) { countG <- 0}
        if(is.na(countT)) { countT <- 0}
        if(is.na(countC)) { countC <- 0}
        if(is.na(countN)) { countN <- 0}
        totalA                     <- totalA + countA
        totalG                 	   <- totalG + countG
        totalT                	   <- totalT + countT
        totalC                	   <- totalC + countC
		
		##for non-coordinates
        countA_non                 <- non_coordinatesVal[which(non_coordinatesNam == "A")][1]
        countG_non                 <- non_coordinatesVal[which(non_coordinatesNam == "G")][1]
        countT_non                 <- non_coordinatesVal[which(non_coordinatesNam == "T")][1]
        countC_non                 <- non_coordinatesVal[which(non_coordinatesNam == "C")][1]
         
        if(is.na(countA_non)) { countA_non <- 0}
        if(is.na(countG_non)) { countG_non <- 0}
        if(is.na(countT_non)) { countT_non <- 0}
        if(is.na(countC_non)) { countC_non <- 0}
        if(is.na(countN_non)) { countN_non <- 0}
        countN_non                  <- countN - (countA_non + countT_non + countG_non + countC_non)
 
        totalA_non                  <- totalA_non + countA_non
        totalG_non                  <- totalG_non + countG_non
        totalT_non                  <- totalT_non + countT_non
        totalC_non                  <- totalC_non + countC_non
        totalN_non                  <- totalN_non + countN_non
                
        numOfClusters                 <- max(as.integer(splitGenomes$StrainCluster_ID),na.rm = TRUE)
        if(is.infinite(numOfClusters)){
                numOfClusters   <- 0
                countA          <- countG <- countT <- countC <- countN <- NA
				countA_non         <- countG_non <- countT_non <- countC_non <- countN_non <- NA
        }
        minClusterLen                <- min(as.integer(splitGenomes$Cluster_Length),na.rm = TRUE)
        if(is.infinite(minClusterLen)){
                minClusterLen   <- NA
        }
        maxClusterLen                 <- max(as.integer(splitGenomes$Cluster_Length),na.rm = TRUE)
        if(is.infinite(maxClusterLen)){
                maxClusterLen   <- NA
        }
        minMutComplex                 <- min(as.integer(splitGenomes$Cluster_Size_Complexes),na.rm = TRUE)
        if(is.infinite(minMutComplex)){
                minMutComplex   <- NA
        }
        maxMutComplex                <- max(as.integer(splitGenomes$Cluster_Size_Complexes),na.rm = TRUE)
        if(is.infinite(maxMutComplex)){
                maxMutComplex   <- NA
        }
        minMutCluster                 <- min(as.integer(splitGenomes$Cluster_Size_Mutations),na.rm = TRUE)
        if(is.infinite(minMutCluster)){
                minMutCluster   <- NA
        }
        maxMutCluster                 <- max(as.integer(splitGenomes$Cluster_Size_Mutations),na.rm = TRUE)
        if(is.infinite(maxMutCluster)){
                maxMutCluster   <- NA
        }
        minPvalue                <- min(as.double(splitGenomes$Cluster_Pvalue),na.rm = TRUE)
        if(is.infinite(minPvalue)){
                minPvalue       <- NA
        }
        minPvalue               <- format(minPvalue, digits = 5)
        maxPvalue               <- max(as.double(splitGenomes$Cluster_Pvalue),na.rm = TRUE)
        if(is.infinite(maxPvalue)){
                maxPvalue       <- NA
        }
        maxPvalue               <- format(maxPvalue, digits = 5)
        temp                    <- as.integer(splitGenomes$Complex_ID)
        complexAdju             <- as.numeric(objects$Adjusted_Complex_number[names(genomeIndata)[jj]])
        noOfMutations           <- length(splitGenomes$Strain_Mutation_ID)
        noOfMutComplexes        <- noOfMutations - complexAdju
        buffer                  <- c(i, complex.threshold, minClusterSize, dist.threshold, pValue, approach, preFilteredMul,filterInfo, clusterFile, numOfClusters, countA, countA_non, countC, countC_non, countG, countG_non, countT, countT_non, countN_non, minClusterLen, maxClusterLen, minMutComplex, maxMutComplex)
        buffer                  <- c(buffer, minMutCluster, maxMutCluster, minPvalue, maxPvalue, objects$Adjusted_Complex_number[names(genomeIndata)[jj]], noOfMutations, noOfMutComplexes )
        bufferTemp              <- ""
        bufferTemp              <- (paste(unlist(buffer),collapse="\t"))
        lapply(bufferTemp, write, summaryFile, ncolumns = length(buffer), append=TRUE, sep="")
        totalClusters           <- totalClusters + numOfClusters

        totalClusterLen         <- c(totalClusterLen,minClusterLen,maxClusterLen)
        totalNumComplexes       <- c(totalNumComplexes,minMutComplex,maxMutComplex )
        totalNumClusters        <- c(totalNumClusters,minMutCluster, maxMutCluster )
        totalminpValue          <- c(totalminpValue, as.double(minPvalue))
        totalmaxpValue          <- c(totalmaxpValue, maxPvalue)
        totalComplexAdj         <- totalComplexAdj + as.numeric(objects$Adjusted_Complex_number[names(genomeIndata)[jj]])
        totalnoMutations        <- totalnoMutations + noOfMutations
        totalnoMutCom           <- totalnoMutCom + noOfMutComplexes
    }
	
	##If any genomes were not processed due to filtering criteria, add them in summary file with "filtered out" value 
    if(length(objects$genomesRest) > 0 ) {
      for (i in 1:length(objects$genomesRest)) {
        buffer                  <- c()
        buffer                  <- c(objects$genomesRest[i], complex.threshold, minClusterSize, dist.threshold, pValue, approach, preFilteredMul,filterInfo, clusterFile)
        recordNotFound          <- rep("filtered_out",21 )
        buffer                  <- c(buffer, recordNotFound)
        bufferTemp              <- ""
        bufferTemp              <- (paste(unlist(buffer),collapse="\t"))
        lapply(bufferTemp, write, summaryFile, ncolumns = length(buffer), append=TRUE, sep="")
      }
    }
    ##write the total summary
	minTotalClusterLen        <- min(totalClusterLen, na.rm = TRUE)
    if(is.infinite(minTotalClusterLen)){
      minTotalClusterLen      <- 0
    }
    maxTotalClusterLen        <- max(totalClusterLen, na.rm = TRUE)
    if(is.infinite(maxTotalClusterLen)){
      maxTotalClusterLen      <- 0
    }
    minTotalNumComplexes        <- min(totalNumComplexes, na.rm = TRUE)
    if(is.infinite(minTotalNumComplexes)){
      minTotalNumComplexes      <- 0
    }
    maxTotalNumComplexes        <- max(totalNumComplexes, na.rm = TRUE)
    if(is.infinite(maxTotalNumComplexes)){
      maxTotalNumComplexes      <- 0
    }
    minTotalNumClusters        <- min(totalNumClusters, na.rm = TRUE)
    if(is.infinite(minTotalNumClusters)){
      minTotalNumClusters      <- 0
    }
    maxTotalNumClusters        <- max(totalNumClusters, na.rm = TRUE)
    if(is.infinite(maxTotalNumClusters)){
      maxTotalNumClusters      <- 0
    }
    minTotalminpValue     <- min(totalminpValue, na.rm = TRUE)
    if(is.infinite(minTotalminpValue)){
      minTotalminpValue      <- 0
    }
    maxTotalminpValue        <- max(totalminpValue, na.rm = TRUE)
    if(is.infinite(maxTotalminpValue)){
      maxTotalminpValue      <- 0
    }
    bufferTemp                  <- ""
    bufferTemp                  <- c("All Genomes",  complex.threshold, minClusterSize, dist.threshold, pValue, approach, preFilteredMul,"", "", totalClusters, totalA, totalA_non, totalC, totalC_non, totalG, totalG_non, totalT, totalT_non, totalN_non, minTotalClusterLen, maxTotalClusterLen)
    bufferTemp                  <- c(bufferTemp, minTotalNumComplexes, maxTotalNumComplexes, minTotalNumClusters, maxTotalNumClusters, minTotalminpValue, maxTotalminpValue, totalComplexAdj, totalnoMutations, totalnoMutCom )
	buffer                      <- ""
    buffer                      <- (paste(unlist(bufferTemp),collapse="\t"))
    lapply(buffer, write, summaryFile, ncolumns = length(buffer), append=TRUE, sep="")
	Sys.chmod(summaryFile, (file.info(summaryFile)$mode | "775"))
}

##Reads values of the parameters from the settings files
##and calls the Core-Script
Process.Series <- function(settingFile) {
	someenv <- new.env()
	startTime = proc.time()
	print(paste("Reading settings file ", settingFile))
	readFile <- read.table(settingFile, header = FALSE, row.names =1,sep="\t",quote="",comment.char = "",check.names = FALSE)
	##Reading the file and building hash key
	mutFileDir            = as.character(readFile["mutFileDir",][1])
	chromFile             = as.character(readFile["chromFile",][1])
	complex.threshold 	  = as.integer(as.character(readFile["complex.threshold",][1]))
	dist.threshold        = as.integer(as.character(readFile["dist.threshold",][1]))
	genome.pos            = as.character(readFile["genome.pos",][1])
	chrom.pos             = as.character(readFile["chrom.pos",][1])
	ref.pos               = as.character(readFile["ref.pos",][1])
	var.Type              = as.character(readFile["variant.pos",][1])
	ref.Allele            = as.character(readFile["ref.Allele",][1])
	clusterFile           = as.character(readFile["clusterResultFile_suff",][1])
	clusterFileLoc        = as.character(readFile["clusterResultFileDir",][1])
	minClusterSize        = as.integer(as.character(readFile["minClusterSize",][1]))
	pValue                = as.double(as.character(readFile["maxPvalue",][1]))
	ChromName.pos         = as.character(readFile["ChromName.pos",][1])
	ChromSize.pos         = as.character(readFile["ChromSize.pos",][1])
	filterInfo            = ""
	preFilteredMul        = 1.0
	filtering             = as.logical(readFile["filtering",][1])
	conservativeApproach  = as.logical(readFile["ConservativeApproach",][1])
	if(conservativeApproach) {
		preFilteredMul    = as.double(as.character(readFile["mutMultiplier",][1]))
	}
	if(filtering) {
		filterInfo        = as.character(readFile["filterInfo",][1])
	}
	else {
		filterInfo        = "No Filter Criteria Selected"
	}

	headerRowIdxMut		  = as.integer(as.character(readFile["headerRowIdxMut",][1]))
	headerRowIdxChr		  = as.integer(as.character(readFile["headerRowIdxChr",][1]))
	
	##Checking if FILESDIR exists, and if it does, to overwrite the value of clusterResultFileDir 
	##with its content and the value of mutFileDir with the content of DATADIR (Suggestion by Les)
	print(exists("FILESDIR"))
	if (exists("FILESDIR")) {
		clusterFileLoc <- FILESDIR
		mutFileDir <- DATADIR
	}
	if (exists("MUTMULTIPLIER")) {
		preFilteredMul <- MUTMULTIPLIER
	}

	##Processing multiple files, one at a time
	mutFileArray <- list.files(mutFileDir,pattern="*.txt$",full.names=TRUE)
	
    if(length(mutFileArray) == 0) { stop("        No MAF file is available in directory.
        Please check the MAF files with *.txt is present in specified directory.")
        warning("No file available")
    }

	print("----------")
	for (i in 1:length(mutFileArray)) {
		print("----------")
		print(paste("Processing file ", mutFileArray[i]))
		objects<-Load.data(mutFile=mutFileArray[i], headerRowIdxMut=headerRowIdxMut, chromFile=chromFile,headerRowIdxChr=headerRowIdxChr, genome.pos=genome.pos, var.type=var.Type, chrom.pos=chrom.pos, ref.pos=ref.pos, complex.threshold=complex.threshold, ChromName.pos=ChromName.pos, ChromSize.pos=ChromSize.pos,filtering = filtering, filterInfo=filterInfo, conservativeApproach);
		ComplexAdjustment <- objects$Adjusted_Complex_number
		if(filtering) {
			objectsOrig<-Load.data(mutFile=mutFileArray[i], headerRowIdxMut=headerRowIdxMut, chromFile=chromFile,headerRowIdxChr=headerRowIdxChr, genome.pos=genome.pos, var.type=var.Type, chrom.pos=chrom.pos, ref.pos=ref.pos, complex.threshold=complex.threshold, ChromName.pos=ChromName.pos, ChromSize.pos=ChromSize.pos,filtering = FALSE, filterInfo=filterInfo, conservativeApproach);
			ComplexAdjustment <- objectsOrig$Adjusted_Complex_number
			rm(objectsOrig)
		}
		objects<-Create.Clusters(objects,dist.threshold=dist.threshold, ref.Allele=ref.Allele, var.Type=var.Type, pValue=pValue, minClusterSize=minClusterSize, conservativeApproach,preFilteredMul ,clusterFile=file.path(clusterFileLoc, paste(gsub(".txt$","",basename(mutFileArray[i])),clusterFile,sep="_")), fileIndex=i, ComplexAdjustment=ComplexAdjustment)
		create.SummaryFile(objects, summaryFile = file.path(clusterFileLoc, paste(gsub(".txt$","",basename(mutFileArray[i])), "_sorted_clusterSum.txt",sep="")),clusterFile= paste(gsub(".txt$","",basename(mutFileArray[i])),clusterFile,sep="_"), complex.threshold=complex.threshold, minClusterSize=minClusterSize, dist.threshold=dist.threshold, pValue=pValue,conservativeApproach= conservativeApproach,preFilteredMul=preFilteredMul, filterInfo = filterInfo  )
	}
	finalTime = proc.time() - startTime
	print(paste("Total elapsed time = ", finalTime["elapsed"]))
	remove(list=ls());
}

##Executing the script.
##ClusterSettings.txt is the setting file, tab-delimited, with all the parameter information
##It should be kept in the same directory as of R code
settingFile <- "ClusterSettings.txt"
Process.Series(settingFile)
