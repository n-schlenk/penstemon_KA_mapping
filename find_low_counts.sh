#!/bin/bash
for i in {1..689}
do
    A=$(head -n $i $1 | tail -n 1 | awk '{print $1}')   # extract individual ID
    B=$(head -n $i $1 | tail -n 1 | awk '{print $2}')   # extract read count
    if [[ $B -lt 1000000 ]]
    then
        echo $A $B                                      # print out ID and count if count is less than one million
    fi
done
