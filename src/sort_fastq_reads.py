# sort_fastq_reads.py
#
# Name: Laura Tung
#
# Usage: python sort_fastq_reads.py --fastq <fastq_file> --outfolder <outfolder> --outfile <outfile> --t <n_cores>
#


#import pdb; pdb.set_trace() # Uncomment to debug code using pdb (like gdb)

import sys
from modules import get_sorted_fastq_for_cluster
import argparse


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--outfile', type=str,  default=False, help='Output filename.')
    parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq file.')
    parser.add_argument('--outfolder', type=str,  default=None, help='Output folder.')
    parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores to use.')
    
    args = parser.parse_args()
    
    args.k = 13
    args.quality_threshold = 7.0
    
    get_sorted_fastq_for_cluster.main(args)
