#!/bin/bash

# get variables from command line
i=$(echo $1)
j=$(echo $2)
new_path=$(echo $3)

# create j + 1
j_plus=$(expr $(echo $j) + 1)

#  get line numbers for start and end of relevant chunk
chunk_start=$(grep -m $j -n "^$" $i | tail -n1 | cut -f1 -d":") # get jth instance of blank line (start of chunk)
chunk_end=$(grep -m $j_plus -n "^$" $i | tail -n1 | cut -f1 -d":") # get j+1th instance of blank line (end of chunk)
if (( $chunk_start == $chunk_end )); then # if it's the last segment, then set $chunk_end to the final line number
    chunk_end=$(cat $i | wc -l )
fi

echo "Chunk start: $chunk_start"
echo "Chunk end: $chunk_end"

# get line with the segment's chr, start, end and strand for file names
hit_count=$(awk "NR==$chunk_start,NR==$chunk_end" $i | grep "^SEQ oryzias_latipes " | wc -l) # get number of matches for HdrR within chunk

echo "Hit count: $hit_count"
echo "Matching lines:"
awk "NR==$chunk_start,NR==$chunk_end" $i | grep "^SEQ oryzias_latipes "

# Extract target line, sequence, tree and data for each HdrR segment match

for k in $( seq 1 $hit_count ) ; do

  # Get target line
  target_line=$(awk "NR==$chunk_start,NR==$chunk_end" $i | grep -m $k "^SEQ oryzias_latipes " | tail -1)

  echo "HdrR metadata line $k: $target_line"

  # create name for data file
  chr=$(echo $target_line | cut -f3 -d" " ) # pull out chr from that line
  segment_start=$(echo $target_line | cut -f4 -d" " ) # pull out start coords from that line
  segment_end=$(echo $target_line | cut -f5 -d" " ) # pull out end coords from that line
  strand=$(echo $target_line | cut -f6 -d" " ) # pull out strand direction from that line

  file_name_data=$(echo $chr\_$segment_start\_$segment_end\_$strand.data.txt) # incorporate into file name

  echo "File name for data $k: $file_name_data"

  # pull out data and copy to data file
  awk "NR==$chunk_start,NR==$chunk_end" $i | egrep -v "^SEQ|TREE|ID|DATA|\/\/|^$" | \
    awk '$1=$1' FS="" OFS="\t" \
    > $new_path/$file_name_data

  # pull out metadata and copy to metadata files
  file_name_seq=$(echo $chr\_$segment_start\_$segment_end\_$strand.seq.txt)
  file_name_tree=$(echo $chr\_$segment_start\_$segment_end\_$strand.tree.txt)
  ## send SEQ lines to file
  awk "NR==$chunk_start,NR==$chunk_end" $i | grep "^SEQ" \
    > $new_path/$file_name_seq
  ## send TREE lines to file
  awk "NR==$chunk_start,NR==$chunk_end" $i | grep "^TREE" \
    > $new_path/$file_name_tree ;

done
