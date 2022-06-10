#! /bin/bash
# created by N. Schlenk
# extracts SNPs by scaffold alignment in VCF before Lep-MAP3
# RUN AS 'sh seperate_vcf_by_scaffold.sh [INPUT VCF]'
# assumes scaffold titles contain 'scaffold[#]'

first_scaffold=0                                                      # define first scaffold
last_scaffold=7                                                       # define last scaffold

header_linecount=$(grep -n -m1 '#CHROM' $1 | grep -o '^..')           # identify line numbers for VCF header section

for i in $(eval echo {$first_scaffold..$last_scaffold})               # for every scaffold...
do
    SNP_count=$(grep scaffold$i $1 | wc -l)                           # count number of scaffold in VCF (must subtract 1 header)
    SNP_count=$(($SNP_count - 1))
    j=$(($i + 1))                                                     # optional adjustment to base 1
    head -n $header_linecount $1 > LG_$j'.vcf'                        # create LG file with VCF headers
    grep scaffold$i $1 | tail -n $SNP_count >> LG_$j'.vcf'            # append SNP lines to new file
done
