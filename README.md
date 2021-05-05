# WaveTAD


## Citation
WaveTAD Ryan Pellow and Josep Comeron

## Installation

Create conda environment
```bash
conda create -n WaveTAD python=3.8
conda activate WaveTAD
```

Install packages
```bash
conda install -c r r-essentials
conda install bwa -c bioconda
conda install samtools -c bioconda
conda install bedtools -c bioconda
conda install hicexplorer -c bioconda -c conda-forge
```

Install R packages
```r
R
install.packages("remotes")
remotes::install_github("HenrikBengtsson/TopDom", ref="master")
install.packages("wavelets")
if(!requireNamespace("BiocManager", quietly=TRUE))
	install.packages("BiocManager")
BiocManager::install("IRanges")
install.packages("iotools")
quit()
```

## Quick Walkthrough
Download FASTQ Files:
[FASTQ1](https://www.dropbox.com/s/5laclm7m8gnr5hw/vaquerizas_3-4hpf_rep1_3R_25Mb_31Mb_1.fastq.gz?dl=0)
[FASTQ2](https://www.dropbox.com/s/sg96jevst7g1ko7/vaquerizas_3-4hpf_rep1_3R_25Mb_31Mb_2.fastq.gz?dl=0)

Download Reference File: 
[Reference](https://www.dropbox.com/s/9w2tnfa650ebh99/dmel-all-chromosome-r6.26.main.chr.10-5-19.fasta?dl=0)

Download Chromosome Sizes:
[Chromosome Sizes](https://www.dropbox.com/s/8dpdtikv35s85ze/chr3R_only_dm6.chrom.sizes?dl=0)

Restriction Enzyme Sequence/Dangling Sequence:
Mbol = GATC/GATC

Index Reference:
```bash
bwa index [Reference]
```

Execute WaveTAD:
```bash
sh WaveTAD.sh [FASTQ1] [FASTQ2] [Reference] [Chromosome Sizes] [Output Directory] [Restriction Enzyme Sequence] [Dangling Sequence]
```





