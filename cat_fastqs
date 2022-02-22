# this script takes every fastq edit file for every unique ID and combines, each ID in one cat file
# ID.txt contains each ID once

#!/bin/bash
for i in $(cat ID.txt)
do
    A=$(find -name "$i.*.gz")
    cat $A >> $i'.cat.fastq.gz'
done
