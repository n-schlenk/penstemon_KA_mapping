#! bin/bash

java -cp [path] ParentCall2 data=pedigree.txt vcfFile=KA.vcf >> parentcalled.vcf
java -cp [path] Filtering2 data=parentcalled.vcf removeNonInformative=1 dataTolerance=0.0000001 >> filtered.vcf
# use separate_scaffold.sh, then run OrderMarkers2 on each scaffold individually
java -cp [path] OrderMarkers2 evaluateOrder=[scaffold map] data=[scaffold VCF] sexAveraged=1 improveOrder=0 > [linkage map for one scaffold]




#### most recent script ####
module load java

for i in {0..7}; do
    java -cp /panfs/pfs.local/home/n752s925/Lep-MAP3/bin OrderMarkers2 evaluateOrder = map_scaffold_$i'.txt' data = scaffold_$i'.vcf' improveOrder=0 sexAveraged=1 grandparentPhase=1 scale=2M/N > scaled_gp_$i'.txt'
done
