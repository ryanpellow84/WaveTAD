# WaveTAD


## Ryan Pellow

## Installation

```bash
conda create -n WaveTAD python=3.8
conda activate WaveTAD
conda install -c r r-essentials
R
install.packages("remotes")
remotes::install_github("HenrikBengtsson/TopDom", ref="master")
install.packages("wavelets")
if(!requireNamespace("BiocManager",quietly=TRUE))
	install.packages("BiocManager")
BiocManager::install("IRanges")
install.packages("iotools")
conda install bwa -c bioconda
conda install samtools -c bioconda
conda install bedtools -c bioconda
conda install hicexplorer -c bioconda -c conda-forge
```

