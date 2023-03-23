# KA_mapping_project
Scripts used to generate recombination map from MSG sequence data.
- Parent 1 = Penstemon kunthii, not inbred (same species as genome assembly)
- Parent 2 = Penstemon amphorellae, not inbred
- F1 (not included) = Cross between kunthii and amphorellae
- F2 (298) = Selfed F1s 

##### Required software:
- bioconda (to install msg_ipyrad environment)
- anaconda
- bwa
- samtools
- bcftools
- java
- python3
- NumPy
- SciPy

##### Required files:
- picard.jar (attached)
- Lep-MAP3 (download from SourceForge)

## Demultiplex and Filter
|||
|-----|-----|
|Starts with:|FASTQ (seq output), TXT (barcodes)|
|Ends with:|FASTQ (concatenated by ID), TXT (stats)|
##### Workflow:
First, you must install and load the msg_ipyrad module (I did this in a conda environment). Use this command to create parameters file.
```
ipyrad -n [title]
```
Then, submit the job (recommended sbatch to cluster)
```
ipyrad -p [parameter file] -s 12 -c 1 -f
```
##### Flags and parameters:
```
ipyrad -n [title]                       creates new parameter file, which must be manually edited
ipyrad -s 12                            only complete steps 1 and 2 of ipyrad workflow (demultiplex, filter)
ipyrad -c 1                             print results (if you sbatch, these will end up in a log file)
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

## Additional Filtering
|||
|-----|-----|
|Starts with:|VCF, pedigree|
|Ends with:|better VCF|

java -cp [filepath to Lep-MAP3 bin] is required before ParentCall2 and Filtering2 because they require Lep-MAP3 information
##### Filters (complete in order):
```
python filter_mq_aldepth_parents.py                     removes markers with low MQ score, low call counts, more than one alt allele, and those which display noninformative/unexpected parental genotypes (parental genotypes must be the last 2 columns)
python filter_refparent.py                              removes markers where the parental genotype of the genome species is not homozygous for the ref allele
python filter_bestsnp.py | filter_extractsnp.py         filters VCF for best SNPs per radtag (assumes that you have already filtered for quality)
head -n37 unfiltered.vcf > pre_pc.vcf                   add original headers back to VCF before using ParentCall2
cat filtered.vcf >> pre_pc.vcf
java -cp <LM filepath> ParentCall2 data=pedigree.txt vcfFile=filtered.vcf > pc.vcf

```
##### Parameters:
`````
filter_mq_aldepth_parents.py                            minIndiv = 50 (minimum number of called genotypes at the SNP)      
filter_mq_aldepth_parents.py                            minMQ = 30 (minimum MQ score for the SNP)
filter_bestsnp.py                                       MinMinor = 1 (number of individuals that will have the minor allele)
ParentCall2 removeNonInfomative=1                       removes markers that are monomorphic or homozygous for both parents
`````

## Map
|||
|-----|-----|
|Starts with:|filtered and parentcalled VCF|
|Ends with:|genetic map coordinates for each scaffold|
##### Workflow:
```
sh separate_scaffold filtered_parentcalled.vcf          generates VCF and map file for each scaffold
sh ordermarkers.sh                                      applies Lep-MAP3 OrderMarkers2 on each scaffold
python parentcall_reformat.py                           create new reference VCF from parentcalled VCF file
python get_coords.py                                    generates coordinates for output of OrderMarkers2 (per scaffold)
```
##### Additional Arguments:
```
OrderMarkers2 sexAveraged=1                             cM values for markers are averaged between sexes
OrderMarkers2 improveOrder=0                            calculate cM values with the orders in which SNPs are aligned
```

### Pre-Graphing
|||
|-----|-----|
|Starts with:|map coordinates|
|Ends with:|data for plotting|
##### Workflow:
```
python find_gc_per_chr.py                                 input fastas, output txt files with GC content per <bin> 
python find_repeats_per_chr.py                            input stockholm file (RepeatModeler) and target repeat element family, output repeat element count per <bin>

```
