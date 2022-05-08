#! /usr/bin/env python3
import re
import os


# cut -f1,2,308,309 KA_called_nodup.vcf | grep -P '/\d' >  amph_kunth_genotypes
# grep -v '^#' KA_called_nodup.vcf > edited_KA_called_nodup.vcf

with open('amph_kunth_genotypes', 'r') as file:
    for line in file:
        geno = line.split()
        amph=(str(re.findall("^...", geno[2]))[2:5])        # amph genotypes
        kunth=(str(re.findall("^...", geno[3]))[2:5])       # kunth genotypes
        if amph != kunth:                                   # eliminate non-informative markers
            if kunth != '1/1':                              # eliminate markers where kunth does not have $
                with open('hand_filtered_parents', 'a') as file2:
                    file2.write(line)


#for i in {1..675841}
#do
#    chr=$(cat fil_pos | head -n $i | tail -n 1 | cut -f1)
#    pos=$(cat fil_pos | head -n $i | tail -n 1 | cut -f2)
#    for j in {1..1619187}
#    do
#        x=$(cat edited_KA_called_nodup.vcf | head -n $j | tail -n 1 | cut -f1)
#        y=$(cat edited_KA_called_nodup.vcf | head -n $j | tail -n 1 | cut -f2)
#        if [ $x == $chr ] && [ $y == $pos ]
#        then
#            cat edited_KA_called_nodup.vcf | head -n $j | tail -n 1 >> filtered_all.vcf
#        fi
#    done
#done
