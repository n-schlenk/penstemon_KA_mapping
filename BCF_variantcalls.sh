#!/bin/bash
#SBATCH --job-name=bcf_var    # Job name
#SBATCH --partition=eeb           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=n752s925@ku.edu      # Where to send mail
#SBATCH --ntasks=1                    # Run one task
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=2-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd $SLURM_SUBMIT_DIR

module load bcftools

ls ./SRGMD > bam_list.txt
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth=100 -f PGA_assembly_trimmed.fasta -b bam_list.txt | bcftools call -vmO v > called.vcf
