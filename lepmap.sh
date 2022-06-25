#! bin/bash

java -cp [path] ParentCall2 data=pedigree.txt vcfFile=KA.vcf >> parentcalled.vcf
java -cp [path] Filtering2 data=parentcalled.vcf removeNonInformative=1 dataTolerance=0.0000001 >> filtered.vcf
# use separate_scaffold.sh, then run OrderMarkers2 on each scaffold individually
java -cp [path] OrderMarkers2 evaluateOrder=[scaffold map] data=[scaffold VCF] sexAveraged=1 improveOrder=0 > [linkage map for one scaffold]

