# gonospy (in progress)

The gonospy (gonos = offspring) package contains a set a functions for filtering, transforming and summarizing high-throughput genomic data. The current version focuses on
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
The first step in any analysis that relies on Illumina data is to check the quality of the raw reads and filter them according to some attributes. The gonospy package follows the
recommendations of Minoche et al. (2011) Genome Biology. The user can define the following parameters:

**indir**: the directory of the zipped fastq.gz raw files  
**qThresh**: The quality threshold below which a base is considered of low quality. Default value 30  
**pQ**: The percentage of each read to be checked for quality score. Default value 0.5 (i.e. the first half of the read will be checked)  
**percQual**: The minimum percentage of total bases above the **qThresh** and within the **pQ** of each read for a read to pass the filtering. Default value 2/3  
**trmBases**: Number of bases to trim from the end of the read. Default value 0  
**repThrs**: The minimum number of identical copies of a read in order to be defined as highly repetitive and excluded from the filtered file  

The following code will walk through all fastq.gz in the 'test' directory, filter them by trimming the last 10 bases per read and save the filtered files in the 'test/cleanFiles' directory together with
a log file containing information about the filtering process. This information will also be printed on the screen. The user should also put a file named 'llumina_Q_scores_phred33.txt' (available from this repository) into the 'test' directory.

```python
gonospy.filterReads(indir = 'test/', trmBases = 10)
```

