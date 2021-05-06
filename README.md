# WaveTAD


## Citation
WaveTAD Ryan Pellow and Josep Comeron

## Installation

Create conda environment
```bash
$ conda create -n WaveTAD python=3.8
$ conda activate WaveTAD
```

Install packages
```bash
$ conda install -c r r-essentials
$ conda install bwa -c bioconda
$ conda install samtools -c bioconda
$ conda install bedtools -c bioconda
$ conda install sra-tools -c bioconda
$ conda install hicexplorer -c bioconda -c conda-forge
```

Install R packages
```r
$ R
> install.packages("remotes")
> remotes::install_github("HenrikBengtsson/TopDom", ref="master")
> install.packages("wavelets")
> if(!requireNamespace("BiocManager", quietly=TRUE))
	install.packages("BiocManager")
> BiocManager::install("IRanges")
> install.packages("iotools")
> quit()
```

## Quick Walkthrough (eta 75mins) (ERR153317[0,1,2,3])
Download FASTQ Files:
[FASTQ1](https://www.dropbox.com/s/5laclm7m8gnr5hw/vaquerizas_3-4hpf_rep1_3R_25Mb_31Mb_1.fastq.gz?dl=0)
[FASTQ2](https://www.dropbox.com/s/sg96jevst7g1ko7/vaquerizas_3-4hpf_rep1_3R_25Mb_31Mb_2.fastq.gz?dl=0)

Download Reference File: 
[Reference](https://www.dropbox.com/s/9w2tnfa650ebh99/dmel-all-chromosome-r6.26.main.chr.10-5-19.fasta?dl=0)

Download Chromosome Sizes:
[Chromosome Size 3R only](https://www.dropbox.com/s/8dpdtikv35s85ze/chr3R_only_dm6.chrom.sizes?dl=0)

Restriction Enzyme Sequence/Dangling Sequence:
Mbol = GATC/GATC

Index Reference:
```bash
$ bwa index [Reference]
```

Execute WaveTAD:
```bash
$ sh WaveTAD.sh [FASTQ1] [FASTQ2] [Reference] [Chromosome Size 3R only] [Output Directory] [Restriction Enzyme Sequence] [Dangling Sequence]
```

Example Output:
| 1 | Chromosome | Start | End | Boundary_5_Prime_Pval | Boundary_3_Prime_Pval | Loop_Pval | Pval |
|---| --- | --- | --- | --- | --- | --- | --- |
| 2 | 3R | 29724331 | 29774976 | 0.0411328806941302 | 0.0481697032295724 | 0.389211597353144 | 0.000771167767436555 |
| 3 | 3R | 28902780 | 28967982 | 0.0468460209849528 | 0.0387259162186036 | 0.0568387748967424 | 0.000103114352438062 |
| 4 | 3R | 29056915 | 29122096 | 0.0498701054774847 | 0.0356704778510274 | 0.000409603200621821 | 0.000000728639239432421 |
| 5 | 3R | 29583384 | 29624658 | 0.034892758869705 | 0.0437696693996931 | 0.000891964347067345 | 0.00000136224766124579 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 578 | 3R | 26060150 | 28348436 | 0.0499273033452901 | 0.0469873372726448 | 0.00158143171430386 | 0.00000370996137707221 |

## Complete Walkthrough (eta ) (SRR1658528)
SRA Number: 1658528

Download Reference File: 
[Reference](https://www.dropbox.com/s/9w2tnfa650ebh99/dmel-all-chromosome-r6.26.main.chr.10-5-19.fasta?dl=0)

Download Chromosome Sizes:
[Chromosome Sizes All](https://www.dropbox.com/s/ach1722yxim6ly9/dm6.chrom.sizes?dl=0)

Restriction Enzyme Sequence/Dangling Sequence:
DpnII = GATC/GATC

Download FASTQ Files:
```bash
$ fastq-dump --split-3 --gzip -A [SRA Number] -O [Output Directory]
```

Index Reference:
```bash
$ bwa index [Reference]
```

Map Reads:
```bash
$ bwa mem -t 56 -E 50 -L 0 [Reference] [FASTQ1] | samtools view --threads 56 -bS - -o [BAM1]
$ bwa mem -t 56 -E 50 -L 0 [Reference] [FASTQ2] | samtools view --threads 56 -bS - -o [BAM2]
```

Calulate Coverage:
```bash
$ sh WaveCov.sh [BAM1] [BAM2] [Chromosome Sizes All]
```

Build Matrices:
```bash
$ hicFindRestSite -f [Reference] -p [Restriction Enzyme Sequence] -o [Restriction Enzyme Sites]
$ hicBuildMatrix -s [BAM1] [BAM2] -bs 1000 -o [1kb Matrix] --skipDuplicationCheck --QCfolder [Output Directory/qcfolder1] --threads 8 -rs [Restriction Enzyme Sites] -seq [Restriction Enzyme Sequence] --danglingSequence [Dangling Sequence]
$ hicBuildMatrix -s [BAM1] [BAM2] -bs 5000 -o [5kb Matrix] --skipDuplicationCheck --QCfolder [Output Directory/qcfolder2] --threads 8 -rs [Restriction Enzyme Sites] -seq [Restriction Enzyme Sequence] --danglingSequence [Dangling Sequence]
$ hicBuildMatrix -s [BAM1] [BAM2] -bs 10000 -o [10kb Matrix] --skipDuplicationCheck --QCfolder [Output Directory/qcfolder3] --threads 8 -rs [Restriction Enzyme Sites] -seq [Restriction Enzyme Sequence] --danglingSequence [Dangling Sequence]
$ hicBuildMatrix -s [BAM1] [BAM2] -bs 25000 -o [25kb Matrix] --skipDuplicationCheck --QCfolder [Output Directory/qcfolder4] --threads 8 -rs [Restriction Enzyme Sites] -seq [Restriction Enzyme Sequence] --danglingSequence [Dangling Sequence]
$ hicBuildMatrix -s [BAM1] [BAM2] -bs 50000 -o [50kb Matrix] --skipDuplicationCheck --QCfolder [Output Directory/qcfolder5] --threads 8 -rs [Restriction Enzyme Sites] -seq [Restriction Enzyme Sequence] --danglingSequence [Dangling Sequence]
```

Detect Loops:
```bash
$ hicDetectLoops -m [1kb Matrix] -o [1kb Loops BedGraph] --chromosomes [Chromosomes] -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
$ hicDetectLoops -m [5kb Matrix] -o [5kb Loops BedGraph] --chromosomes [Chromosomes] -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
$ hicDetectLoops -m [10kb Matrix] -o [10kb Loops BedGraph] --chromosomes [Chromosomes] -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
$ hicDetectLoops -m [25kb Matrix] -o [25kb Loops BedGraph] --chromosomes [Chromosomes] -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
$ hicDetectLoops -m [50kb Matrix] -o [50kb Loops BedGraph] --chromosomes [Chromosomes] -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
```

Correct Matrices:
```bash
$ hicCorrectMatrix correct --matrix [1kb Matrix] --chromosomes [Chromosomes] -o [1kb Corrected Matrix]
$ hicCorrectMatrix correct --matrix [5kb Matrix] --chromosomes [Chromosomes] -o [5kb Corrected Matrix]
$ hicCorrectMatrix correct --matrix [10kb Matrix] --chromosomes [Chromosomes] -o [10kb Corrected Matrix]
$ hicCorrectMatrix correct --matrix [25kb Matrix] --chromosomes [Chromosomes] -o [25kb Corrected Matrix]
$ hicCorrectMatrix correct --matrix [50kb Matrix] --chromosomes [Chromosomes] -o [50kb Corrected Matrix]
```

Normalize Matrices:
```bash
$ hicNormalize -m [1kb Corrected Matrix] --normalize norm_range -o [1kb Normalized Matrix]
$ hicNormalize -m [5kb Corrected Matrix] --normalize norm_range -o [5kb Normalized Matrix]
$ hicNormalize -m [10kb Corrected Matrix] --normalize norm_range -o [10kb Normalized Matrix]
$ hicNormalize -m [25kb Corrected Matrix] --normalize norm_range -o [25kb Normalized Matrix]
$ hicNormalize -m [50kb Corrected Matrix] --normalize norm_range -o [50kb Normalized Matrix]
```

Convert Matrices to Tables:
```bash
$ cooler dump -t pixels --header --join [1kb Normalized Matrix] -o [1kb Contacts Table]
$ cooler dump -t pixels --header --join [5kb Normalized Matrix] -o [5kb Contacts Table]
$ cooler dump -t pixels --header --join [10kb Normalized Matrix] -o [10kb Contacts Table]
$ cooler dump -t pixels --header --join [25kb Normalized Matrix] -o [25kb Contacts Table]
$ cooler dump -t pixels --header --join [50kb Normalized Matrix] -o [50kb Contacts Table]
```

Generate Diamond Area Scores:
```bash
$ R CMD BATCH --no-save --no-restore "--args [1kb Contacts Table] [1kb TopDom Matrix] [1kb TopDom Scores] 1000 [Chromsome Sizes All]" topdom.R
$ R CMD BATCH --no-save --no-restore "--args [5kb Contacts Table] [5kb TopDom Matrix] [5kb TopDom Scores] 5000 [Chromsome Sizes All]" topdom.R
$ R CMD BATCH --no-save --no-restore "--args [10kb Contacts Table] [10kb TopDom Matrix] [10kb TopDom Scores] 10000 [Chromsome Sizes All]" topdom.R
$ R CMD BATCH --no-save --no-restore "--args [25kb Contacts Table] [25kb TopDom Matrix] [25kb TopDom Scores] 25000 [Chromsome Sizes All]" topdom.R
$ R CMD BATCH --no-save --no-restore "--args [50kb Contacts Table] [50kb TopDom Matrix] [50kb TopDom Scores] 50000 [Chromsome Sizes All]" topdom.R
```

Use WaveTAD to Call TADs for each Chromsome Independently (small genomes and for this example):
```bash
$ awk -v myvar="[Chromosome]" '$1==myvar' [Left Coverage File] > [Left Coverage Chromosome]
$ awk -v myvar="[Chromosome]" '$1==myvar' [Right Coverage File] > [Right Coverage Chromosome]
$ R CMD BATCH --no-save --no-restore "--args [Left Coverage Chromosome] [Right Coverage Chromosome] [WaveTAD Results] [1kb TopDom Scores] [5kb TopDom Scores] [10kb TopDom Scores] [25kb TopDom Scores] [1kb Loops BedGraph] [5kb Loops BedGraph] [10kb Loops BedGraph] [25kb Loops BedGraph] [Chromosome] [Size of Biggest Chromsome]" WaveTAD.R  WaveTAD_[Chromsome].Rout
```

Use WaveTAD to Call TADs for each Chromsome Independently (big genomes):
```bash
$ awk -v myvar="[Chromosome]" '$1==myvar' [Left Coverage File] > [Left Coverage Chromosome]
$ awk -v myvar="[Chromosome]" '$1==myvar' [Right Coverage File] > [Right Coverage Chromosome]
$ R CMD BATCH --no-save --no-restore "--args [Left Coverage Chromosome] [Right Coverage Chromosome] [WaveTAD Results] [5kb TopDom Scores] [10kb TopDom Scores] [25kb TopDom Scores] [50kb TopDom Scores] [5kb Loops BedGraph] [10kb Loops BedGraph] [25kb Loops BedGraph] [50kb Loops BedGraph] [Chromosome] [Size of Biggest Chromsome]" WaveTAD.R  WaveTAD_[Chromsome].Rout
```


