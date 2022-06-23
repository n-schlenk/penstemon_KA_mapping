#! /bin/bash
# index reference genome, align each fastq to reference genome (results in SAM file)
# AddOrReplaceReadGroups, translate to BAM, sort
# fastq files must be in 'KA_fastqs' directory
# created by N. Schlenk

module load bwa
module load java
module load picard
module load samtools

bwa index -p kunthii PGA_assembly_trimmed.fasta
for i in $(ls KA_fastqs/)
do
    ID=$(echo $i | grep -Po '^.{5}')
    bwa mem kunthii KA_fastqs/$i > $ID'_bwa_output.sam'
    java -jar picard.jar AddOrReplaceReadGroups \
        I=$ID'_bwa_output.sam'\
        O=$ID'_withRG.sam' \
        RGLB=$ID \
        RGPL=ILLUMINA \
        RGPU=barcode \
        RGSM=$ID
    samtools view -b $ID'_withRG.sam' | samtools sort > $ID'.bam'
done



