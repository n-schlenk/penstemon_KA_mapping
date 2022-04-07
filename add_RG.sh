module load java
module load picard

for bam_file in $(ls ./BAM)
do
    ID=$(echo $bam_file | grep -Po '^.{5}')
    java -jar picard.jar AddOrReplaceReadGroups I=./BAM/$bam_file O=$ID'_withRG$
done



