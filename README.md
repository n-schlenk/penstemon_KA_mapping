# KA_mapping_project
Scripts used to generate recombination map from MSG sequence data.  
1. Demultiplex and filter    
2. Align to reference genome
3. Call variants  
4. Filter
5. LepMap  

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

