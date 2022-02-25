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
-------------------------------------------------------------------------------------------------
# after running this script, determine line count for each ID in fin_edits
for i in $(cat ID.txt)                                  # for every ID...
do
    A=$(find -name "$i.*.gz" | grep "fin_edits")        # find filepath in fin_edits with ID in title
    B=$(wc -l $A)                                       # count number of lines for filepath
    echo $i $B >> fin_edits_linecount.txt               # print ID and linecount to output file
done
-------------------------------------------------------------------------------------------------
# lastly, you need to compare the files
comm cat_edits_linecount.txt fin_edits_linecount.txt    # if cat went as planned, every ID should be in the 3rd column
