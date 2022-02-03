The code released in this package has been developed under CentOS Linux and successfully tested/executed up to CentOS version 7.8 with R version up to 3.3.2 and Python version up to 2.7.13. Execution under other OS and language versions may require code modifications.

Configurable parameters
scriptDir/SCRIPTDIR - directory location of scripts and additional config files. The config files contain some additional parameters that can be modified but should be static in a production pipeline
dataDir/DATADIR - directory location of input files for the first script of pipeline (ClusterFinder)
filesDir/FILESDIR - directory location of output files for ClusterFinder and sequential inputs/outputs for the remaining scripts
refFasta/REF_FA - file location of genome reference sequence in FASTA format, its index file REF_FA.fai (created by 'samtools faidx') must be in the same directory
mutMultiplier/MUTMULTIPLIER - fraction of genome sequenced: 0.01 for whole-exome MAF and 1.0 for whole-genome MAF (used only to identify mutation clusters)
The mixed case versions are the names of the arguments passed to the main() function of control.R used to start the pipeline, resulting in setting the corresponding uppercase global variables. Their default values in control.R should be edited to match the installation setup.

Pipeline workflow
The pipeline's first script, ClusterFinder, starts in the data directory DATADIR and creates its outputs in the work directory FILESDIR where the remaining scripts operate. ClusterFinder has no restrictions on the names of the data files it accepts. The remaining scripts take only the files with the proper extension created by the preceding step in the pipeline.
The execution of the entire pipeline is coordinated by the script control.R in the "code" subdirectory. The preferred and recommended way of running the pipeline is to execute this script by sourcing it into an interactive R interpreter, setting the locations of the package directory and reference FASTA file, then calling the main() function:
$ R
> PARENTDIR = "/home/username/pmacd/"
> REFPATH = "/home/username/reference_genome/hg19.fa"
> source(paste0(PARENTDIR, "code/control.R"))
> main()
If needed, control.R can also start the pipeline by being executed as a standalone script using the Rscript or R CMD BATCH mechanism, which would require code modifications that are beyond the scope of this README (see https://stackoverflow.com/questions/33400312/use-main-function-in-r ).

To accelerate the execution of the pipeline for a given MAF file with multiple motifs, the retrieval of the sequence contexts around the mutations from the genome reference as well as adding the cluster annotations can be performed only once, since the outputs of these two steps are not dependent on a particular signature motif. The retrieved sequences are saved in a *seq output file when the pipeline for the first motif is run. This file and the *anz2.txt output can be copied or linked to the destination output directories for the subsequent motifs prior to launching their respective pipelines - with the source("ClusterFinder_standAlone.R") step commented out, and they will be detected and used bypassing the slow steps of adding cluster information and retrieving the same sequences from the genome again.
The standard version of the pipeline performs analysis of the tCw motif. In addition, this release includes configuration files for custom motifs tCa, rtCa, ytCa, yCn, nCg, rCg, yCg, nTt, rTt, and yTt that can be activated by manually editing the CUSTOMMOTIFS and CUSTOMTYPE switches in control.R. Any other custom motifs can be configured and analyzed by creating appropriate analogous customization files. The required syntax for representing the mutation motifs is to use an uppercase letter for the mutated base and lowercase letters for the surrounding bases.
The signature analysis in mutation clusters will evaluate signature-specific mutagenesis in all clusters as well in clusters containing only C- and/or G- mutated bases, and therefore the complete analysis is performed only for those custom motifs in which the mutated base is either C or G (for instance yCn).  For the custom motifs in which the mutated base is either A or T (for instance nTt),  the cluster analysis still evaluates motif-specific mutagenesis in all clusters, but produces zero values for clusters containing only mutations in C and/or G.

Input Formatting Requirements.
Mutation Annotation File (MAF): tab-delimited table of tumor-specific mutation calls per row in The Cancer Genome Atlas (TCGA) format, based on either whole exome sequencing (WES) or whole genome sequencing (WGS). Non-ASCII binary characters (>127 decimal) are not allowed in any column. The file is not required to be a TCGA MAF (see https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/ ), however the following column names and the value syntax and capitalization in these columns should be fixed even if analysis is performed for a non-tumor and/or non-human mutation calls:
Sample ID - column name: "Tumor_Sample_Barcode"; value syntax: not defined, each sample ID must be unique.
Chromosome - column name: "Chromosome"; value syntax: numbers and letters only, preferably prefixed with "chr" (human example: chr1-chr22, chrX, chrY, chrM or 1-22, X, Y, M)
(Note that P-MACD filters out mutation calls from chromosome assemblies whose names are not present in the index of the genome reference.
This syntax extends to genomes of other species by default as long as they follow the above naming conventions of the human genome.
It is crucial that exactly the same naming convention be used in three places: the input MAF, the config file *.chrom.sizes for ClusterFinder, and the genome reference.
Other naming conventions can be accommodated but have not been tested and may require code modifications.).
Position - column name: "Start_position"; value syntax: integer numbers
Reference allele - column name: "Reference_Allele"; value syntax; A, T, C, G (nucleotides in capitals.
Tumor allele - column name: "Tumor_Seq_Allele2"; value syntax: A, T, C, G (nucleotides in capitals).
Variant type - column name: "Variant_Type"; value syntax: SNP (single base substitutions), other - not defined.

Reference Genome: nucleotide sequence of the reference genome with one entry per chromosome in FASTA format.

Fraction of Genome Sequenced: 0.01 for whole-exome MAF and 1.0 for whole-genome MAF. (Used only to identify mutation clusters).

Output Data files
The pipeline files are described in detail in three APOBEC-specific Readme_*.txt files or in three Readme_ CUSTOMMOTIF_*.txt files , addressing the file content, and the columns of the extended MAF formats (one mutation per row) and the summary formats (one genome sample per row).

Each row in a summary file provides values calculated for a given sample. The last row provides totals across all samples. The following values (all present in "*_sorted_sum_all_fisher_Pcorr.txt" file) are the most useful for analyses of oligonucleotide-motif centered mutation signatures:

Fold-enrichment for a signature - "*_enrich" column;

Absolute number of signature mutations in a sample - "*+revcomp" column;

Fraction of signature mutations in a sample - [*+revcomp]_per_mut" column;

Benjamini-Hochberg-corrected p-values - "BH_Fisher_p-value_*" column.

Minimum estimate of the number of signature mutations caused by a signature-specific mutation process in a sample - "*_MutLoad_MinEstimate" column

The subdirectory "data" contains a small test MAF WES file that has been processed by the standard tCw version of the pipeline and by the CUSTOMMOTIF yCn pipeline. The results were saved to the "work_done" subdirectory.
Due to R's restrictions, binary characters cannot be present in data files and must be removed prior to launching the pipeline.

Background
The package was originally designed for the analysis of Pattern of Mutagenesis by APOBEC Cytidine Deaminases (P-MACD) and described at Roberts, S.A. et al., An APOBEC cytidine deaminase mutagenesis pattern is widespread in human cancers, Nature Genetics 45:970-976 (2013) for the analysis of canonical APOBEC mutation motif tCw to tTw or to tGw.  It was then extended to include various oligonucleotide-centered mutational motifs (CUSTOMMOTIFS) as one of configurable parameters (Chan, K et al. An APOBEC3A hypermutation signature is distinguishable from the signature of background mutagenesis by APOBEC3B in human cancers, Nature Genetics 47:1067-1072(2015) and Saini et al. The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts, PLOS Genet 12(10): e1006385. doi:10.1371/journal.pgen.1006385.  The outputs for the canonical APOBEC mutational motif tCw contain a number of APOBEC-specific analyses and are described in a set of three Readme_*.txt files.  The outputs for the user-defined configurable oligonucleotide-centered mutational motifs (CUSTOMMOTIFS) are very similar except that they replace the canonic APOBEC motif and its relevant text parts in file and in column names with the oligonucleotide sequences generated based on the corresponding CUSTOMMOTIF config.


