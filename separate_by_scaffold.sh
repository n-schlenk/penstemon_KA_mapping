#! /bin/bash
# created by N. Schlenk
# separates VCF into scaffolds, creates map for each scaffold using order in VCF
# RUN AS 'sh separate_by_scaffold.sh [INPUT VCF]'

first_scaffold=0
last_scaffold=7

for i in $(eval echo {$first_scaffold..$last_scaffold})
do
    grep '##' $1 > scaffold_$i'.vcf'
    SNP_count1=$(grep scaffold$i $1 | wc -l)
    SNP_count2=$(($SNP_count1 - 1))
    grep scaffold$i $1 | tail -n $SNP_count2 >> scaffold_$i'.vcf'
    for n in $(eval echo {1..$SNP_count1})
    do
        echo $n >> map_scaffold_$i'.txt'
    done
done

