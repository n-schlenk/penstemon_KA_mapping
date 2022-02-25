#!/bin/bash
#SBATCH --job-name=bwa_alignment    # Job name
#SBATCH --partition=eeb           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=n752s925@ku.edu      # Where to send mail
#SBATCH --ntasks=1                    # Run one task
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=0-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd $SLURM_SUBMIT_DIR

module load bwa                                                           # load module

bwa index -p kunthii ./PGA_assembly.fasta                                 # create index with prefix kunthii
for i in $(ls ../fin_edits/KA_edits)                                      # for every KA fastq file...
do
    ID=$(echo $i | grep -o ^[^.]*)                                        # grab ID
    bwa mem kunthii ../fin_edits/KA_edits/$i > $ID'_bwa_output.sam'       # align fastq to kunthii genome and output to ID_bwa_output.sam
done

