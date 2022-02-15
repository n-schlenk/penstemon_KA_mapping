# INPUT one variable
  # file with ID (duplicates allowed) and counts
# OUTPUT
  # file with ID (no duplicates) and concatenated read count


#!/bin/bash
N=$(cat $1 | wc -l)																	                                  # determine number of lines
for i in $(seq 1 $N)
do
    A=$(awk '{print $1}' $1 | head -n $i | tail -n 1)                                 # extract one ID
    B=$(grep -c $A $1)                                                                # determine number of duplicates for ID
    if [[ $B == 1 ]]                                                                  # if it is not a duplicate...
    then
        D=$(awk '{print $1, $2}' $1 | head -n $i | tail -n 1)                         # save D as ID and count as-is
    else                                                                              # if it is a duplicate (2 or more)
        C=$(grep $A $1 | awk '{print $2}' | awk '{total += $1} END {print total}')    # combine number of counts
        D=$(echo $A $C)                                                               # save D as ID and combined counts for that ID
    fi
    echo $D > cat_counts.txt                                                          # print out D (ID and counts) to new file
done
