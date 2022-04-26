#!/bin/bash
#SBATCH --job-name=bwa_stuff    # Job name
#SBATCH --partition=eeb           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=n752s925@ku.edu      # Where to send mail
#SBATCH --ntasks=1                    # Run one task
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=2-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd $SLURM_SUBMIT_DIR

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
    java -jar /home/n752s925/scratch/picard.jar MarkDuplicates /home/n752s925/scratch/BWA/$ID'.bam' > /home/n752s925/scratch/BWA/$ID'_sortedRGMD.bam'
done



