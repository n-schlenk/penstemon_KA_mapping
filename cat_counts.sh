#!/bin/bash
for i in {1..1193}
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
    echo $D                                                                           # print out D (ID and counts)
done
