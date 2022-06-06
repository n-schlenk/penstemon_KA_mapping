# fixes ID errors in the original edited files (EXAMPLE), merges them into one big txt file with ID, raw read count, and filtered read count
# created by N. Schlenk and K. Brubeck

#! /bin/bash
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

-------------------------------------------------------------------------------------------------

# concatenates all fastq files with the same ID into one file under that ID
# assumes that titles are fixed like they were fixed in the above code, ID.txt contains each ID once
# created by N. Schlenk

#!/bin/bash
for i in $(cat ID.txt)
do
    A=$(find -name "$i.*.gz")
    cat $A >> $i'.cat.fastq.gz'
done

-------------------------------------------------------------------------------------------------

# created by N. Schlenk
# INPUT 2 arguments
  # file with IDs (no dups) in one columns and total read count in another (in this case, the 3rd col)
  # minimum read threshold
# OUTPUT 1 file
  # file with ID and count reads that did not meet threshold

#!/bin/bash
N=$(cat $1 | wc -l)                                     # calculate number of lines in file
for i in $(seq 1 $N)
do
    A=$(head -n $i $1 | tail -n 1 | awk '{print $1}')   # extract individual ID
    B=$(head -n $i $1 | tail -n 1 | awk '{print $3}')   # extract read count
    if [[ $B -lt $2 ]]
    then
        echo $A $B >> not_enough_reads.txt              # print ID and count to file if count is less than given threshold
    fi
done

-------------------------------------------------------------------------------------------------

# verify that your fastq files were concatenated properly by counting lines in each file
# pre-concatenated fastq files (multiple per ID) in cat_edits, concatenated fastqs (one per ID) in fin_edits
# created by N. Schlenk

for i in $(cat ID.txt)                                  # for every ID...
do
    A=$(find -name "$i.*.gz" | grep "cat_edits")        # find filepath(s) in cat_edits with ID in title
    for j in $(echo A)                                  # for every filepath... 
    do
        wc -l $j >> $i'counts.txt'                      # count number of lines and append to <ID>counts.txt
    done
    B=$(awk '{sum+=1} END {print sum}' $i'.counts.txt') # calculate line count sum for each <ID>counts.txt
    echo $i $B >> cat_edits_linecount.txt               # print ID and total linecount to output file
    rm $i'counts.txt'                                   # delete <ID>counts.txt file
done

# after running this script, determine line count for each ID in fin_edits
for i in $(cat ID.txt)                                  # for every ID...
do
    A=$(find -name "$i.*.gz" | grep "fin_edits")        # find filepath in fin_edits with ID in title
    B=$(wc -l $A)                                       # count number of lines for filepath
    echo $i $B >> fin_edits_linecount.txt               # print ID and linecount to output file
done

# lastly, you can compare the files
comm cat_edits_linecount.txt fin_edits_linecount.txt    # if cat went as planned, every ID should be in the 3rd column
