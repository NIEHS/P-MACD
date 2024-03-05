# P-MACD
P-MACD is a package originally designed for the analysis of patterns of mutagenesis by APOBEC Cytidine Deaminases that has been extended to support a variety of oligonucleotide-centered mutational motifs.

## Requirements
P-MACD was developed under CentOS Linux and successfully tested/executed in CentOS version 7.8 with R 3.3.2 and Python 2.7.13.

These legacy requirements can be provided - as tested in February 2024 under Rocky Linux release 8.8 - using the Miniconda package manager following these steps:
1. Install the Miniconda system following the instructions at: https://docs.anaconda.com/free/miniconda/miniconda-install/
2. Activate the base environment: <br>
`
eval "$(/home/username/miniconda3/bin/conda shell.bash hook)"
`
3. Create and activate an empty environment: <br>
`
conda create --name legacyR
`<br>
`
conda activate legacyR
`
4. Install R in the new environment:<br>
`
conda install r=3.3.2 -c conda-forge
`
5. Install Python:<br>
`
conda install python=2.7 -c conda-forge
`
6. Install the packages required by R (data.table) and Python (numpy)

Analogous solutions can be developed for the Docker or Apptainer platforms.


## Installation
No installation is required. Simply copy the contents of this repository to a convenient user-accessible location.

## Usage
P-MACD is most easily executed from within a interactive R session. Begin by defining `PARENTDIR` and `REFPATH` variables that store the locations of the P-MACD software and reference genome FASTA file, respectively:
```
> PARENTDIR = "/home/username/pmacd/"
> REFPATH = "/home/username/reference_genome/hg19.fa"
```
Launch the pipeline by importing the `control.R` script found in the `code` directory, and executing the `main` function:
```
> source(paste0(PARENTDIR, "code/control.R"))
> main()
```
See [README_1st.txt](README_1st.txt) for additional details.

## Author
P-MACD was authored by Leszek J. Klimczak, Ph.D.

## Citing P-MACD
On publication of any studies whose findings are in part derived from the output of the P-MACD software package, please include the following citation:
* An APOBEC cytidine deaminase mutagenesis pattern is widespread in human cancers. Roberts SA, Lawrence MS, Klimczak LJ, Grimm SA, Fargo D, Stojanov P, Kiezun A, Kryukov GV, Carter SL, Saksena G, Harris S, Shah RR, Resnick MA, Getz G, Gordenin DA.
Nat Genet. 2013 Sep;45(9):970-6. doi: 10.1038/ng.2702.

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
