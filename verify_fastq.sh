# this script takes every filtered and demultiplexed fastq file for each unique ID and combines them, each ID in one cat file
# ID.txt contains each ID once
-------------------------------------------------------------------------------------------------
#!/bin/bash
for i in $(cat ID.txt)
do
    A=$(find -name "$i.*.gz")
    cat $A >> $i'.cat.fastq.gz'
done


# INPUT 2 arguments
  # file with IDs (no dups) in one columns and total read count in the other
  # minimum read threshold
# OUTPUT 1 file
  # file with ID and count reads that did not meet threshold
-------------------------------------------------------------------------------------------------
#!/bin/bash
N=$(cat $1 | wc -l)                                     # calculate number of lines in file
for i in $(seq 1 $N)
do
    A=$(head -n $i $1 | tail -n 1 | awk '{print $1}')   # extract individual ID
    B=$(head -n $i $1 | tail -n 1 | awk '{print $2}')   # extract read count
    if [[ $B -lt $2 ]]
    then
        echo $A $B >> not_enough_reads.txt              # print ID and count to file if count is less than given threshold
    fi
done


# the purpose of this workflow is to verify that your fastq files were concatenated properly by counting lines in each file
# pre-concatenated fastq files (multiple per ID) in cat_edits, concatenated fastqs (one per ID) in fin_edits
-------------------------------------------------------------------------------------------------
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

# lastly, you need to compare the files
comm cat_edits_linecount.txt fin_edits_linecount.txt    # if cat went as planned, every ID should be in the 3rd column
