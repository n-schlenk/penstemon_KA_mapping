# KA_mapping_project
Scripts used to generate recombination map from MSG sequence data.  
Steps include demultiplexing/filtering, aligning to reference genome, variant calling, and generating map with Lep-MAP3.
##### Required software:
- bioconda (to install msg_ipyrad environment)
- anaconda
- bwa
- samtools
- bcftools
- java

##### Required files
- picard.jar (attached)
- Lep-MAP3 (download from SourceForge)

## Demultiplex and Filter
|||
|-----|-----|
|Starts with:|FASTQ (seq output), TXT (barcodes)|
|Ends with:|FASTQ (concatenated by ID), TXT (stats)|
##### Workflow:
First, you must install and activate msg_ipyrad environment and use ipyrad to create parameters file
```
ipyrad -n [title]
```
Then, submit the job (recommended sbatch to cluster)
```
conda activate msg_ipyrad
ipyrad -p [parameter file] -s 12 -c 1 -f
```
##### Flags and parameters:
```
ipyrad -n [title]                       creates new parameter file, which must be manually edited
ipyrad -s 12                            only complete steps 1 and 2 of ipyrad workflow (demultiplex, filter)
ipyrad -c 1                             print results (if you sbatch, these will end up in the serial_test file)
ipyrad -f                               forces ipyrad to run steps, even if they have already been completed (applies when some data is missing)                         
```

## Align to Reference Genome
|||
|-----|-----|
|Starts with:|FASTQ (one per ID), FASTA (ref genome)|
|Ends with:|BAM (one per ID)|
##### Workflow:
```
bwa index | bwa mem | java AddOrReplaceReadGroups | samtools view -b | samtools sort > [final BAM]
```
##### Flags and parameters:
```
samtools view -b                        translates the SAM to BAM (binary)
```
## Call Variants
|||
|-----|-----|
|Starts with:|BAM (one per ID)|
|Ends with:|VCF (single)|
##### Workflow:
```
bcftools mpileup | bcftools call > [VCF]
```
##### Flags and parameters:
```
bcftools mpileup -Ou                     output to uncompressed BCF
bcftools mpileup -I                      exclude indels (only include SNPs)
bcftools mpileup -a FORMAT/AD            outputs allelic depth in addition to the defaults (genotype and phred-scaled genotypic likelihood) 
bcftools mpileup --max-depth=100         include a maximum of 100 reads per sample per SNP
```
```
bcftools call -v                         output variant sites only 
bcftools call -m                         use alterate model for multiallelic and rare-variant calling (recommended)
bcftools -O v                            output to uncompressed VCF
```

## Lep-MAP3
|||
|-----|-----|
|Starts with:|VCF, TXT (pedigree)|
|Ends with:|TXT (genetic map coordinates)|
##### Workflow:
```
java -cp [filepath to Lep-MAP3 bin] is required before each step so that the server knows where to pull Lep-MAP3 information
```
```
ParentCall2 data=[pedigree TXT] vcfFile=[VCF] > [parentcalled VCF]
Filtering2 data=[parentcalled VCF] > [filtered VCF]
SeparateChromosomes2 data=[filtered VCF] > [map TXT]
JoinSingles2All map=[map TXT] data=[filtered VCF] > [joined map TXT]
OrderMarkers2 map=[joined map TXT] data=[filtered VCF] > [final map TXT]
```
##### Flags and parameters:
```
Filtering2 removeNonInfomative=1            removes markers that are monomorphic or homozygous for both parents
Filtering2 dataTolerance=0.0000001          removes distorted markers at give p-value threshold (applies HWE)
SeparateChromosomes2 lodLimit=5             separates linkage groups based on LOD limit (markers with LOD higher than the provided limit are grouped together)
JoinSingles2All lodLimit=4                  distributes singles into existing linkage groups using a less strict LOD limit
OrderMarkers2 chromosome=1                  orders markers within linkage groups (sometimes it is easier to do this one linkage group at a time, thus the chromosome flag)
OrderMarkers2 sexAveraged=1                 cM values for markers are averaged between sexes
OrderMarkers2 grandparentPhase=1            phases data with consideration to the presence of grandparents in the dataset
OrderMarkers2 evaluateOrder=[previously ordered map] should only be used on a pre-existing ordered map. It is recommended that this step is run on your compeleted ordered map to produce a more confident ordering. I typically do this 2-5 times.
```
