#!/bin/bash
# call variants on BAM files (must be in 'bam_files' directory), output one VCF
# created by N. Schlenk

module load bcftools
ls ./bam_files > bam_list.txt
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth=100 -f PGA_assembly_trimmed.fasta -b bam_list.txt | bcftools call -vmO v > called.vcf
