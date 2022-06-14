#!/bin/bash
for filename in /home/leonid/PycharmProjects/python/minor_2022/project/*.bed; do
    sort -k1,1 -k2,2n $filename > sorted_$filename
    bedtools merge -i sorted_$filename  -c 5 -o max > merged_$filename.bed
done

