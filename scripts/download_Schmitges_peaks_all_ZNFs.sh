#! /bin/bash

# list of URLs for the MACS 1e-5 hg19 peaks
mkdir data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/
# Loop through each line in the file
while IFS= read -r line; do
    # Extract the specified column
    col=$(echo "$line" | awk '{print $2}')
    
    # Download the file using wget
    wget -q -P data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/ "https://chip-atlas.dbcls.jp/data/hg19/eachData/bed05/${col}.05.bed"
done < data/schmitges_all_znfs_chips/schmitges_all_znfs.txt
