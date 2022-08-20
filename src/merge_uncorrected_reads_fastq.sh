#!/bin/bash
#
# Laura Tung
#
# Usage: merge_uncorrected_reads_fastq.sh <source_dir> 

source_dir=$1

touch merged_uncorrected_reads.fastq
for i in {54339..332299..1}
do
    cat ${source_dir}/${i}.fastq >> merged_uncorrected_reads.fastq
done

