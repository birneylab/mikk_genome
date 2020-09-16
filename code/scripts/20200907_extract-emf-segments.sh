#!/bin/bash

# get variables from command line
i=$(echo $1)
j=$(echo $2)
bname_short=$(echo $3)

#  get line that states which segment it corresponds to in HdrR
target_line=$(grep -m $j "^SEQ oryzias_latipes " $i | tail -n1 ) # get jth instance of the line with segment coords for HdrR
chr=$(echo $target_line | cut -f3 -d" " ) # pull out chr from that line
segment_start=$(echo $target_line | cut -f4 -d" " ) # pull out start coords from that line
segment_end=$(echo $target_line | cut -f5 -d" " ) # pull out end coords from that line
strand=$(echo $target_line | cut -f6 -d" " ) # pull out strand direction from that line
#echo -e "$chr\t$segment_start\t$segment_end\t$strand"
file_name_data=$(echo $chr\_$segment_start\_$segment_end\_$strand.data.txt) # incorporate into file name

# pull out data and copy to data file
chunk_start=$(grep -m $j -n "^DATA" $i | tail -n1 | cut -f1 -d":") # get jth instance of ^DATA (start of chunk)
chunk_end=$(grep -m $j -n "//" $i | tail -n1 | cut -f1 -d":") # get jth instance of // (end of chunk)
## pull out just the data, and insert tabs between each character
awk "NR==$(expr $(echo $chunk_start) + 1),NR==$(expr $(echo $chunk_end) - 1)" $i | \
  awk '$1=$1' FS="" OFS="\t" \
  > emfs/segmented/$bname_short/$file_name_data

# pull out metadata and copy to metadata files
file_name_seq=$(echo $chr\_$segment_start\_$segment_end\_$strand.seq.txt)
file_name_tree=$(echo $chr\_$segment_start\_$segment_end\_$strand.tree.txt)
blank_line_start=$(grep -m $j -n ^$ $i | tail -n1 | cut -f1 -d":") # get line number of jth instance of blank line (start of each segment)
counter_plus_1=$(expr $(echo $j) + 1) # get j + 1
blank_line_end=$(grep -m $counter_plus_1 -n ^$ $i | tail -n1 | cut -f1 -d":") # get line number of jth +1 instance of blank line (end of segment)
## send SEQ lines to file
awk "NR==$blank_line_start,NR==$blank_line_end" $i | grep "^SEQ" \
  > emfs/segmented/$bname_short/$file_name_seq
## send TREE lines to file
awk "NR==$blank_line_start,NR==$blank_line_end" $i | grep "^TREE" \
  > emfs/segmented/$bname_short/$file_name_tree 
