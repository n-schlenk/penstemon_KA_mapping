## Project Overview

**Title:** The _Penstemon kunthii_ draft genome: Integrating a genetic map with assembled sequence data  
**Purpose:** To verify the existing _Penstemon kunthii_ genome assembly using a linkage map derived from a cross between _P. kunthii_ and _P. amphorellae_. This was my master's thesis project at the University of Kansas.  
**Summary:** We generated MSG sequence reads from an F2 population, which were then used to make a genetic map for each of the 8 _Penstemon_ chromosomes. Then, I investigated other genomic elements including GC content and repeat element distribution. This repository holds the code used to complete this project, as well as my notes/instructions.
  
**Required software:** anaconda (msg_ipyrad, NumPy), bwa, samtools, BCFtools, java, python3, Lep-MAP3  
**Required files:** picard.jar (attached), sequences/barcodes (multiplexed), pedigree

## Step 1: Demultiplex and Filter
|||
|-----|-----|
|Starts with:|FASTQ (seq output), TXT (barcodes)|
|Ends with:|FASTQ (concatenated by ID), TXT (stats)|
##### Workflow:
First, you must load the msg_ipyrad (bioconda) environment. Use this to create parameters file. Then run/submit.
```
ipyrad -n <title>
ipyrad -p <parameter file> -s 12 -c 1 -f
```
##### Flags and parameters:
```
ipyrad -n <title>                       # creates new parameter file, which must be manually edited
ipyrad -s 12                            # only complete steps 1 and 2 of ipyrad workflow (demultiplex, filter)
ipyrad -c 1                             # print results (if you sbatch, these will end up in a log file)
ipyrad -f                               # forces ipyrad to run steps, even if they have already been completed (applies when some data is missing)                         
```

## Step 2: Align to Reference Genome
|||
|-----|-----|
|Starts with:|FASTQ (one per ID), FASTA (ref genome)|
|Ends with:|BAM (one per ID)|
##### Workflow:
```
bwa index | bwa mem | java AddOrReplaceReadGroups | samtools view | samtools sort > <final BAM file>
```
##### Flags and parameters:
```
samtools view -b                        # translates the SAM to BAM (binary)
```
## Step 3: Call SNPs
|||
|-----|-----|
|Starts with:|BAM (one per ID)|
|Ends with:|VCF (single)|
##### Workflow:
```
bcftools mpileup | bcftools call > <VCF>
```
##### Flags and parameters:
```
bcftools mpileup -Ou                     # output to uncompressed BCF
bcftools mpileup -I                      # exclude indels (only include SNPs)
bcftools mpileup -a FORMAT/AD            # outputs allelic depth in addition to the defaults (genotype and phred-scaled genotypic likelihood) 
bcftools mpileup --max-depth=100         # include a maximum of 100 reads per sample per SNP
bcftools call -v                         # output variant sites only 
bcftools call -m                         # use alterate model for multiallelic and rare-variant calling (recommended)
bcftools -O v                            # output to uncompressed VCF
```

## Step 4: Filter, ParentCall, OrderMarkers
|||
|-----|-----|
|Starts with:|VCF, pedigree|
|Ends with:|better VCF|

I used several different filters on our dataset before generating the linkage maps. These filters and recommended parameters are included in the above code. Before ordering markers on the filtered VCF, you must apply ParentCall2. I ordered the markers **with** and **without** the genome as a reference. Below, I show the code used **with** the genome. The instructions for ordering **without** the genome are described in **project_steps_explained.md** above.

_java -cp <filepath to Lep-MAP3 bin> is required before calling Lep-MAP3 modules (ParentCall2 and OrderMarkers2)_
##### Workflow:
```
ParentCall2 data=pedigree.txt vcfFile=filtered.vcf > pc.vcf
OrderMarkers2 evaluateOrder=map.txt data=pc.vcf > order.txt
```
##### Flags and Parameters:
```
ParentCall2 removeNonInfomative=1                       # removes markers that are monomorphic or homozygous for both parents
OrderMarkers2 sexAveraged=1                             cM values for markers are averaged between sexes
OrderMarkers2 improveOrder=0                            calculate cM values with the orders in which SNPs are aligned
```
