# index reference genome, align each fastq to reference genome (results in SAM file)
# AddOrReplaceReadGroups, translate to BAM, sort, optional MarkDuplicates
# index resulting BAM file
# created by N. Schlenk

#!/bin/bash
module load bwa
module load java
module load picard
module load samtools

bwa index -p kunthii /home/n752s925/scratch/BWA/PGA_assembly_trimmed.fasta
for i in $(ls /home/n752s925/scratch/fin_edits/KA_edits)
do
    ID=$(echo $i | grep -Po '^.{5}')
    bwa mem kunthii /home/n752s925/scratch/fin_edits/KA_edits/$i > /home/n752s925/scratch/BWA/$ID'_bwa_output.sam'
    java -jar /home/n752s925/scratch/picard.jar AddOrReplaceReadGroups I=/home/n752s925/scratch/BWA/$ID'_bwa_output.sam' O=/home/n752s925/scratch/BWA/$ID'_withRG.sam'
    samtools view -b /home/n752s925/scratch/BWA/$ID'_withRG.sam' | samtools sort > /home/n752s925/scratch/BWA/$ID'.bam'
    # java -jar /home/n752s925/scratch/picard.jar MarkDuplicates /home/n752s925/scratch/BWA/$ID'.bam' > /home/n752s925/scratch/BWA/$ID'_sortedRGMD.bam'
    samtools index $ID'.bam'
done



