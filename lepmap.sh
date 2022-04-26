#! bin/bash

java -cp /panfs/pfs.local/work/hileman/n752s925/bin ParentCall2 data=pedigree.txt vcfFile=KA_called.vcf >> parentcalled.vcf
java -cp /panfs/pfs.local/work/hileman/n752s925/bin Filtering2 data=parentcalled.vcf removeNonInformative=1 dataTolerance=0.0000001 >> filtered.vcf
java -cp /panfs/pfs.local/work/hileman/n752s925/bin SeparateChromosomes2 data=filtered.vcf lodLimit=5 > map.txt
java -cp /panfs/pfs.local/work/hileman/n752s925/bin JoinSingles2All map=map.txt data=filtered.vcf lodLimit=4 > joined_map.txt
java -cp /panfs/pfs.local/work/hileman/n752s925/bin OrderMarkers2 map=joined_map.txt data=filtered.vcf chromosome=PGA_scaffold0__171_contigs__length_68081433 > ordered1.txt

