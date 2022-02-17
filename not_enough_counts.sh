# INPUT 2 arguments
  # file with ID (no dups) and count reads (cat)
  # minimum read threshold
# OUTPUT 1 file
  # file with ID and count reads that did not meet threshold


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
