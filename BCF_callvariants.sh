# call variants on BAM files, output one VCF
# created by N. Schlenk


#!/bin/bash

module load bcftools
ls ./SRGMD > bam_list.txt
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth=100 -f PGA_assembly_trimmed.fasta -b bam_list.txt | bcftools call -vmO v > called.vcf
