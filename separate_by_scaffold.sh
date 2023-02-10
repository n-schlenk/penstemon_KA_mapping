#! /bin/bash
# created by N. Schlenk
# separates VCF into scaffolds, creates map for each scaffold using order in VCF
# RUN AS 'sh separate_by_scaffold.sh [INPUT VCF]'

first_scaffold=0
last_scaffold=7

for i in $(eval echo {$first_scaffold..$last_scaffold})
do
    rm map_scaffold$i'.txt'
    head -n7 $1 > scaffold_$i'.vcf'
    SNP_count=$(grep scaffold$i $1 | wc -l)
    grep scaffold$i $1 >> scaffold_$i'.vcf'
    for n in $(eval echo {1..$SNP_count})
    do
        echo $n >> map_scaffold_$i'.txt'
    done
done
