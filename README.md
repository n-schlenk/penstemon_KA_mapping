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

```
testing testing
```

## Align to Reference Genome
|||
|-----|-----|
|Starts with:|FASTQ (one per ID), FASTA (ref genome)|
|Ends with:|BAM (one per ID)|

## Call Variants
|||
|-----|-----|
|Starts with:|BAM (one per ID)|
|Ends with:|VCF (single)|
