#!/bin/bash
mkdir cat_edits                                                             # make directory for concatenated fastq files
A=$(awk '{print $1}' $1 | uniq )                                            # assemble file with one ID per line, no duplicates
N=$(awk '{print $1}' $1 | uniq | sed '/^$/d' | wc -l)                       # calculate number of unique IDs
for i in $(seq 1 $N)                                                        # start loop from 1 to last unique ID in list
do
    B=$(echo $A | cut -d ' ' -f $i)                                         # define ID within loop
    C=$(grep $B $1 | awk '{print $3}')                                      # extract location(s) for ID
    D=$(grep -c $B $1)                                                      # count number of duplicates
    if [[ $D == 1 ]]                                                        # if there are no duplicates
    then                                                                    # then copy file to new directory
        cp ./$C'_edits/'$B'.trimmed_R1_.fastq.gz' ./cat_edits
    elif [[ $D == 2 ]]                                                      # if there are 2 duplicates
    then                                                                    # then cat the, and move to new directory
        E=$(echo $C | cut -d ' ' -f 1)
        F=$(echo $C | cut -d ' ' -f 2)
        cat ./$E'_edits/'$B'.trimmed_R1_.fastq.gz' ./$F'_edits/'$B'.trimmed_R1_.fastq.gz' > ./cat_edits/$B'_cat_edts_trimmed_R1_.fastq.gz'
    elif [[ $D == 3 ]]                                                      # if there are 3 duplicates
    then                                                                    # then cat them and move to new directory
        E=$(echo $C | cut -d ' ' -f 1)
        F=$(echo $C | cut -d ' ' -f 2)
        G=$(echo $C | cut -d ' ' -f 3)
        cat ./$E'_edits/'$B'.trimmed_R1_.fastq.gz' ./$F'_edits/'$B'.trimmed_R1_.fastq.qz' ./$G'_edits/'$B'.trimmed_R1_.fastq.gz' > ./cat_edits/$B'_cat_edits_trimmed_R1_.fastq.gz'
    else
        echo $B has more than 3 duplicates                                  # this could go on forever, change if needed
    fi
done
