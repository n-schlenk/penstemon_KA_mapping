#! /bin/bash

for file in $(ls ./BAM)
do
    new_name=$(echo $file | sed 's/.sam//g')
    mv ./BAM/$file ./BAM/$new_name
done
