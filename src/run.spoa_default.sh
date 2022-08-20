#!/bin/bash
#
# Usage: run.spoa_default.sh <in_fastq> <outdir> <out_fasta>
#
# <in_fastq>: input spoa set fastq. e.g. *.fixed.fastq
# <outdir>: output directory. e.g. *.spoa_output
# <out_fasta>: output alignment matrix fasta. e.g. *.spoa_sma_output.fa

in_fastq=$1
outdir=$2
out_fasta=$3


bin_dir=/home/ltung/.local/bin

if [ ! -d $outdir ]
then
    mkdir $outdir
fi

# Run spoa
{ time $bin_dir/spoa -l 0 -r 2 $in_fastq > $outdir/$out_fasta; } 2> $outdir/time.log
