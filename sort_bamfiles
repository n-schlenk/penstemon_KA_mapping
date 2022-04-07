module load samtools
for bam_file in $(ls ./RG_bamfiles)
do
    ID=$(echo $bam_file | grep -Po '^.{5}')
    samtools sort $bam_file > ./RG_bamfiles/$ID'_sorted_RG.bam'
done
