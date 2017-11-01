# gonospy (in progress)

The gonospy (gonos = offspring) package contains a set a functions for filtering, transforming and summarizing (i.e. calculate diversity indices) high-throughput genomic data. The current version focuses on
data generated via sequencing of reduced representation libraries (RRLs, e.g. RAD-seq) using Illumina technology. Future versions will include long-read sequencing technologies as well.

## Prerequisites
The following python packages should also be already installed:
pandas, numpy

## Examples

Download the gonospy package and define it's directory
```python
import sys 
sys.path.append('pathTogonospyFolder/')

import gonospy
```

### Filter raw sequences / reads
The first step in any analysis that relies on Illumina data is to check the quality of the raw reads and filter them according to some criteria. The gonospy package follows the
recommendations of Minoche et al. (2011) Genome Biology. The user can define the following parameters:

**indir**: the directory of the zipped fastq.gz raw files  
**qThresh**: The quality threshold below which a base is considered of low quality (30)  
**pQ**: The percentage of each read to be checked for quality score (0.5, i.e. the first half of the read will be checked)  
**percQual**: The minimum percentage of total bases above the **qThresh** and within the **pQ** of each read for a read to pass the filtering (2/3)  
**trmBases**: Number of bases to trim from the end of the read (0)  
**repThrs**: The minimum number of identical copies of a read in order to be defined as highly repetitive and excluded from the filtered file (1000)  

The following code will walk through all fastq.gz in the 'test' directory, filter them by trimming the last 10 bases per read and save the filtered files in the 'test/cleanFiles' directory together with
a log file containing information about the filtering process. This information will also be printed on the screen. The user should also put a file named 'illumina_Q_scores_phred33.txt' (available from this repository) into the 'test' directory.

```python
gonospy.filterReads(indir = 'test/', trmBases = 10)
```

### Filter VCF files
The following function filters VCF files (https://samtools.github.io/hts-specs/VCFv4.2.pdf) based on several criteria defined by the user.
Alleles below the thresholds described below will be excluded and the genotype will be called as homozygous for the alternative allele. The following parameters can be defined:

**indir**: the directory of the .VCF files  
**refDepthThrs**: read coverage threshold for reference allele (5)  
**varDepthThrs**: read coverage threshold for the alternative allele (5)  
**maxAlleles**: maximum number of alternative alleles (1) (only binary states are allowed for the moment)  
**missSampThrs**: percentage of missing samples allowed (0.2)  
**totRefFreqThrs** minimum frequency of the reference allele across all genotypes (0.05)  
**totVarFreqThrs** minimum frequency of the variant allele across all genotypes (0.05)  
**totCovThrs** minimum frequency of the reference allele (coverage) in each genotype (0.05, accounts for extreme differences in read coverage between the reference and the alternative allele)  
**varCovThrs** minimum frequency of the variant allele (coverage) in each genotype (0.05, same as above)  
**locus**. if false (default) the chromosome number will be printed, else the loci number (for VCF files produced by STACKS; http://catchenlab.life.illinois.edu/stacks/)

```python
gonospy.filterVCF(indir = 'test/', RefDepthThrs = 10, MissSampThrs = 0.4, VarDepthThrs = 10, TotRefFreqThrs = 0.1, TotCovThrs = 1)
```
 
 ### Convert VCF
After filtering the VCF file, we may want to use our SNP data either to identify potential genetic structure across our samples, or reconstruct the phylogenetic relationships
among the samples. The following function converts the VCF to either of the three formats: 'phylip', 'structure' or 'fastStructure'. The user can modify the following parameters:

**infile**: input VCF file  
**outFormat**: 'phylip', 'structure' or 'fastStructure'  
**missChr**: the character to use for missing data. The default is ? but it should be changed to -9 for fastStructure  
**oneSNP**: export only one random SNP per locus (True) or all SNPs per locus (False)  
**nOfRep**: number of repetitions when exporting only one random SNP per locus (1)   
**excludeInv**: exclude invariant sites for downstream analyses using RAxML (False)  
**randomAllele**: for heterozygous genotypes in phylip format, the user can either use the IUPAC codes, or randomly choose one of the two states (False)  
**locusID**: if True (default), the program will try to get the locus information from the SND ID column (3rd column in VCF file). Otherwise it will use the chromosome ID (1st column)  
**populationFiles**: in structure format, the second column denotes the population ID where the sample belongs. If None (default value), the program will attempt to
retrieve the population ID from the sample name (e.g. sampleID_popID). Otherwise the user must supply a txt (tab delimited file) with one column for the sample names
and a second for the population IDs  

```python
#VCF to phylip
gonospy.converVCF(infile=infile, outFormat='phylip', excludeInv = False, oneSNP=True)
```
