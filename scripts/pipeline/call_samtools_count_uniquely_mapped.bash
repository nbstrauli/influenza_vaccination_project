#!/bin/bash

# $1 = the input filepath for the BAM file that I want to count the number of mapped reads in it
# $2 = the output filepath to write the number of mapped reads to

#This gives the number of reads that mapped to the ref genome, but don't
#count those who map more then once (well, only count them once that is).
#Don't use this for large BAM files though, cause the sorting step is
#memory intensive.

samtools view $1 | cut -f1 | sort | uniq | wc -l > $2

exit
