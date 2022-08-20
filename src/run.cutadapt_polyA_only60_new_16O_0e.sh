#!/bin/bash
#
# Usage: run.cutadapt_polyA_only60_new_16O_0e.sh <longreads_data_filename> 
#
# <longreads_data_filename>: fastq file with adapters already trimmed (by porechop or pychopper).
# 
# This is to further trim the poly-A sequences only.

longreads_data=$1


bin_dir=/home/ltung/.local/bin


# Run cutadapt
echo "Running cutadapt..."
/usr/bin/time -p $bin_dir/cutadapt -a "A{60}" -O 16 -m 20 -e 0 -j 80 -o $(basename $longreads_data .fastq)_poly60_16O_0e.fastq $longreads_data &> run.cutadapt60_new_16O_0e.log
echo "Done cutadapt"

