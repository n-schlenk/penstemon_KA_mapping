#! /bin/bash

# run as 'sh stats.sh'

# awk '{print $1}' [stats file]                           # look at every stats file individually and identify errors in sample format, write code for each file as needed
mkdir edit_stats
    ## aug2017_KA: kaF2 -> KA, parent ids
cat raw_stats/aug2017_KAL1_stats.txt | sed 's/kaF2-/KA/g' | sed 's/kunth[0-9]/kunth/g' | sed 's/amph[0-9]/amph/g' > edit_stats/aug2017_KAL1_edited.txt
cat raw_stats/aug2017_KAL2_stats.txt | sed 's/kaF2-/KA/g' | sed 's/kunth[0-9]/kunth/g' | sed 's/amph[0-9]/amph/g' > edit_stats/aug2017_KAL2_edited.txt
    ## aug2017_NB: KA -> NB, parent ids
cat raw_stats/aug2017_NBL1_stats.txt | sed 's/KA/NB/g' | sed 's/barb[0-9]/barb/g' | sed 's/neom[0-9]/neom/g' > edit_stats/aug2017_NBL1_edited.txt
cat raw_stats/aug2017_NBL2_stats.txt | sed 's/KA/NB/g' | sed 's/barb[0-9]/barb/g' | sed 's/neom[0-9]/neom/g' > edit_stats/aug2017_NBL2_edited.txt
    ## aug2019: remove '-'
sed 's/-//g' raw_stats/aug2019_P1_stats.txt > edit_stats/aug2019_P1_edited.txt
    ## dec2017_P1: F2 -> NB, parent ids
cat raw_stats/dec2017_P1L1_stats.txt | sed 's/F2/NB/g' | sed 's/barb[0-9]/barb/g' | sed 's/neom[0-9]/neom/g' > edit_stats/dec2017_P1L1_edited.txt
cat raw_stats/dec2017_P1L2_stats.txt | sed 's/F2/NB/g' | sed 's/barb[0-9]/barb/g' | sed 's/neom[0-9]/neom/g' > edit_stats/dec2017_P1L2_edited.txt
    ## dec2017_P2: kaF2 -> KA, parent ids
cat raw_stats/dec2017_P2L1_stats.txt | sed 's/ka/KA/g' | sed 's/F2//g' | sed 's/kunth[0-9]/kunth/g' | sed 's/amph[0-9]/amph/g' > edit_stats/dec2017_P2L1_edited.txt
cat raw_stats/dec2017_P2L2_stats.txt | sed 's/ka/KA/g' | sed 's/F2//g' | sed 's/kunth[0-9]/kunth/g' | sed 's/amph[0-9]/amph/g' > edit_stats/dec2017_P2L2_edited.txt
    ## march2019_P1: NB -> NB(0), parent ids
cat raw_stats/march2019_P1_stats.txt | sed 's/^NB\([0-9]\{2\}\W\)/NB0\1/g' | sed 's/kunthii../kunth/g' > edit_stats/march2019_P1_edited.txt
    ## march2019_P2: CAW -> CAW(.)
cat raw_stats/march2019_P2_stats.txt | sed 's/CAW\([0-9]\{2\}\)/CAW\1./g' | sed 's/\(CAW[0-9.]*\)[A|B]/\1/g' | sed 's/CAW49./CAW49.0/g' | sed 's/CAW41.PL5[A|B]/CAW41.5/g' > edit_stats/march2019_P2_edited.txt

cp raw_stats/april2021* edit_stats                         # move rest of files to new folder, even if they didn't need editing
cp raw_stats/septnov2021* edit_stats

for file in $(ls edit_stats)                               # this loop will merge all of the stats files into one big one with duplicate IDs
do
    sed 1d edit_stats/$file >> cat_stats.txt
    echo -e "\n" >> cat_stats.txt
    sort cat_stats.txt > sorted_cat_stats.txt
done

rm cat_stats.txt                                           # remove intermediate file, create list of all unique IDs
awk '{print $1}' sorted_cat_stats.txt | uniq > unique.txt

for i in $(cat unique.txt)                                 # this loop creates a file with a sum of raw reads and a sum of filtered reads for every ID
do
    cat sorted_cat_stats.txt | grep $i | awk '{raw+=$2}{filtered+=$7} END{print $1,raw,filtered}' >> ID_raw_fil.txt
done

awk '{if ($3 < 1000000) print $1, $3}' ID_raw_fil.txt > not_enough_reads.txt    # returns individuals without a million reads 
awk '{if ($3 > 1000000) print $1, $3}' ID_raw_fil.txt > yes_enough_reads.txt    # returns individuals with a million or more reads
