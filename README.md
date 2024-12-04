# PhaseFinderDC

## Overview
PhaseFinderDC is a workflow adapted from the original [PhaseFinder](https://github.com/XiaofangJ/PhaseFinder) algorithm. Like PhaseFinder, it is designed to detect DNA inversion mediated phase variation in bacterial genomes using metagenomic sequencing data. PhaseFinderDC includes several modifications that optimize its performance for metagenomic sequencing data from mixtures of defined communities (DC). A defined community is one in which all bacterial strains present are known, and (ideally) have a high quality genome sequence. 

Like original PhaseFinder, PhaseFinderDC works by identifying regions flanked by inverted repeats, mimicking their inversion in silico, and identifying regions where sequencing reads support both orientations. New modifications in PhaseFinderDC include: (1) use of bowtie2 instead of bowtie for read alignment, (2) inclusion of full genome sequence in bowtie2 index instead of only sequences of inverted repeat regions, and (3) MAPQ filtering of bowtie2 read alignments to ignore cases of potential cross-mapping. These modifications allow PhaseFinderDC to identify invertons with high sensitivity and specificity in defined communities, even when closely related strains are present.

## Prerequisites
+ [Biopython](https://biopython.org/)
+ [pandas](https://pandas.pydata.org)
+ [samtools](http://samtools.sourceforge.net/) (>=1.4)
+ [bowtie2](https://github.com/BenLangmead/bowtie2)
+ [einverted](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [bedops](https://bedops.readthedocs.io/en/latest/)
+ [bedtools](https://bedtools.readthedocs.io/en/latest/)
+ [click](https://click.palletsprojects.com/en/8.1.x/)

To install PhaseFinderDC
```
git clone https://github.com/xiaofanjin/PhaseFinderDC.git
cd PhaseFinderDC
conda env create --file environment.yml
conda activate PhaseFinderDC
```

## Quick Start
As with the original PhaseFinder, all you need to get started with PhaseFinderDC is a genome (in fasta format) you would like to search for invertible DNA regions and genomic sequencing data (preferrably Illumina in fastq format) from the same organism, or metagenomic sequencing data from a sample containing the organism (preferrably Illumina in fastq format). For communities with multiple strains, concatenate the fasta files for each strain's genome into a single input genome fasta file.

To test PhaseFinderDC, you can use the example files supplied with the original PhaseFinder (genome: test.fa, genomic data: p1.fq, p2.fq) 

Example:
```
# Identify regions flanked by inverted repeats 
python PhaseFinderDC.py locate -f ./data/test.fa -t ./data/test.einverted.tab -g 15 85 -p 

# Mimic inversion
python PhaseFinderDC.py create -f ./data/test.fa -t ./data/test.einverted.tab -s 1000 -i ./data/test.ID.fasta

# Identify regions where sequencing reads support both orientations 
python PhaseFinderDC.py ratio -i ./data/test.ID.fasta -1 ./data/p1.fq -2 ./data/p2.fq -o ./data/out
```

If successful, the output will be in data/out.ratio.txt

In this example, there is one real example of an invertible DNA region "am_0171_0068_d5_0006:81079-81105-81368-81394" because only this region has reads supporting both the F and R orientation. 

---

## Tutorial
### 1. Generate a position table of regions flanked by inverted repeats 
Users can identify inverted repeats using the "PhaseFinderDC.py locate" command, or generate their own table. Note "PhaseFinderDC.py locate" is identical to the original "PhaseFinder.py locate"

#### 1.1. Generate the position table with the PhaseFinder script
```
Usage: PhaseFinderDC.py locate [OPTIONS]

  Locate putative inverted regions

Options:
  -f, --fasta PATH        Input genome sequence file in fasta format
                          [required]
  -t, --tab PATH          Output table with inverted repeats coordinates
                          [required]
  -e, --einv TEXT         Einverted parameters, if unspecified run with
                          PhaseFinder default pipeline
  -m, --mismatch INTEGER  Max number of mismatches allowed between IR pairs,
                          used with -einv (default:3)
  -r, --IRsize INTEGER    Max size of the inverted repeats, used with -einv
                          (default:50)
  -g, --gcRatio MIN MAX   The minimum and maximum value of GC ratio
  -p, --polymer           Remove homopolymer inverted repeats
  --help                  Show this message and exit.
```

##### Input: A fasta file containing the genome sequence
##### Output: A table file containing the postion information of invereted repeats in the genome

##### Examples:
* Run the default PhaseFinderDC locate parameters
```
python PhaseFinderDC.py locate -f ./data/test.fa -t ./data/test.einverted.tab 
```
* Run the default PhaseFinder locate parameters and remove inverted repeats with GC content lower than 15% and higher than 85% or with homopolymers
```
python PhaseFinderDC.py locate -f ./data/test.fa -t ./data/test.einverted.tab -g 15 85 -p 
```
* Run with the specified einverted parameters "-maxrepeat 750 -gap 100 -threshold 51 -match 5 -mismatch -9" 
```
python PhaseFinderDC.py locate -f ./data/test.fa -t ./data/test.einverted.tab -e "-maxrepeat 750 -gap 100 -threshold 51 -match 5 -mismatch -9" 
```


#### 1.2. Generate the position table with other tools
You can identify regions flanked by inverted repeats directly with tools such as [einverted](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html) and [palindrome](http://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html). 

![inverted repeats](https://github.com/XiaofangJ/PhaseFinder/blob/master/IR.png)

Prepare the output into the following format:

A table file with five columns (tab delimited):

 Column name | Explanation                                                   |
-------------|---------------------------------------------------------------|
 Scaffold    | The scaffold or contig name where the inverted repeat is detected
 pos A       | The start coordinate of the first inverted repeat (0-based)
 pos B       | The end coordinate of the first inverted repeat (1-based)
 pos C       | The start coordinate of the second inverted repeat (0-based)
 pos D       | The end coordinate of the second inverted repeat (1-based)

---
### 2. Mimic inversion in silico to create a database of inverted sequences
```
Usage: PhaseFinderDC.py create [OPTIONS]

  Create inverted fasta file

Options:
  -f, --fasta PATH         Input genome sequence file in fasta format
                           [required]
  -t, --tab PATH           Table with inverted repeat coordinates  [required]
  -s, --flanksize INTEGER  Base pairs of flanking DNA on both sides of the
                           identified inverted repeats  [required]
  -i, --inv PATH           Output path of the inverted fasta file  [required]
  -p, --threads INTEGER    Number of threads to use for building bowtie2 index [default 1]
  --help                   Show this message and exit.
```

#### Input
* The position table from step 1

#### Output
* A fasta file containing inverted (R) and non-inverted (F) putative invertible DNA regions flanked by sequences of specified length (bowtie2 indexed). Note that unlike the original PhaseFinder, PhaseFinderDC also includes in this fasta file sequences of intervening regions between putative invertible DNA regions, thus the final bowtie2 index covers the full sequence of original input genomes.
* A table file (with suffix ".info.tab") describing the location of inverted repeats in the above fasta file

---
### 3. Align sequence reads to inverted sequence database and calculate the ratio of reads aligning to the F or R orienation. 
```
Usage: PhaseFinderDC.py ratio [OPTIONS]

  Align reads to inverted fasta file

Options:
  -i, --inv PATH         Input path of the inverted fasta file  [required]
  -1, --fastq1 PATH      First pair in fastq  [required]
  -2, --fastq2 PATH      Second pair in fastq  [required]
  -p, --threads INTEGER  Number of threads
  -q, --minmapq INTEGER  bowtie2 mapQ threshold to filter read alignments
  -a, --bt2args TEXT     bowtie2 arguments
  -o, --output TEXT      Output prefix  [required]
  -kb, --keepbam         Keep bam file output
  --help                 Show this message and exit.
```

#### Input
* Output from step 2
* fastq file of genomic or metagenomic sequence used to verify DNA inversion
* Number of threads used for bowtie2 alignment and samtools process. PhaseFinderDC uses bowtie2 instead of bowtie
* MAPQ filter threshold (new feature in PhaseFinderDC): bowtie2 alignments with MAPQ score less than this threshold are ignored when counting reads (default = 30)
* flexible bowtie2 align arguments: string specifying flags to pass on to bowtie2 aligner, use double quotes (default = "--very-sensitive")
#### Output
* A table file (with suffix ".ratio.txt") containing the reads that supporting either R or F orientation of invertible DNA - formatting follows that of the original PhaseFinder algorithm

 Column name | Explanation                                                                 |
-------------|-----------------------------------------------------------------------------|
Sequence     | Putative invertible regions(Format:Scaffold:posA-posB-posC-posD)
Pe_F         | The number of reads supprting the F orientation with paired-end information
Pe_R         | The number of reads supprting the R orientation with paired-end information
Pe_ratio     | Pe_R/(Pe_F + Pe_R). The percent of reads supporting the R orientation with the paired-end method
Span_F       | The number of reads supporting the F orientation spanning the inverted repeat by at least 10 bp on either side
Span_R       | The number of reads supporting the R orientation spanning the inverted repeat by at least 10 bp on either side
Span_ratio   | Span_R/(Span_F + Span_R). The percent of reads supporting the R orientation with the spanning method. 

As with the original PhaseFinder algorithm, true invertible regions have reads supporting both the F and R orientation. We recommend combining the information from both the paired-end (Pe) and spanning (Span) methods to find valid invertible DNA regions. Our default is similar to that recommended by the original PhaseFinder, to classify a region as invertible if Pe_R >= 5 and Span_R >= 5 and Pe_F >= 5 and Span_F >= 5. 

## Citation
Jin X, et al. Comprehensive profiling of genomic invertons in defined gut microbial community reveals associations with intestinal colonization and surface adhesion, *bioRxiv* (2024) 

## Citation for original PhaseFinder
Jiang X, Hall AB, et al. Invertible promoters mediate bacterial phase variation, antibiotic resistance, and host adaptation in the gut, *Science* (2019) [DOI: 10.1126/science.aau5238](http://science.sciencemag.org/content/363/6423/181)
