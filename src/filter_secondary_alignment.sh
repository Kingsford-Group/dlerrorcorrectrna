#!/bin/bash
#
# Usage: filter_secondary_alignment.sh <sam_file> 
#

sam_file=$1

# get the header
grep "^@" $sam_file > tmp1

# filter out secondary alignments
samtools view -F 256 $sam_file > tmp2

cat tmp1 tmp2 > tmp2.1

# extract supplementary alignments to identify chimeric alignment reads
samtools view -f 2048 tmp2.1 > suppl_algn
awk '{print $1}' suppl_algn > chimeric_readnames
sort chimeric_readnames | uniq > unique_chimeric_readnames

# filter out supplementary alignments
samtools view -F 2048 tmp2.1 > tmp3

cat tmp1 tmp3 > primary_$sam_file

rm tmp1 tmp2 tmp3 tmp2.1 suppl_algn

