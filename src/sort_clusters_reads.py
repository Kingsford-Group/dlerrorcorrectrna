# sort_clusters_reads.py
#
# Name: Laura Tung
#
# Usage: python sort_clusters_reads.py <final_clusters.tsv> <fastq_files_dir> <n_cores>
#
# <final_clusters.tsv>: cluster file outputted by isONclust: final_clusters.tsv
# <fastq_files_dir>: directory containing the fastq files for all clusters: fastq_files

#import pdb; pdb.set_trace() # Uncomment to debug code using pdb (like gdb)

import sys
import numpy as np
import pandas as pd

import subprocess


def sort_reads_for_clusters(cluster_df, fastq_files_dir, n_cores, bin_dir):
    
    min_cluster_size = 3
    max_spoa_size = 800
    
    top_lines = max_spoa_size*4
    
    num_clusters = cluster_df.iloc[-1,0] + 1
    print("num_clusters:", num_clusters)
    
    cluster_size_list = []
    for cluster_id in range(num_clusters):
        # find the cluster size
        single_cluster_df = cluster_df[cluster_df.iloc[:,0] == cluster_id]
        cluster_size = single_cluster_df.shape[0]
        
        if cluster_size >= min_cluster_size:
            cluster_size_list.append([cluster_id, cluster_size])
            
            cluster_fastq = fastq_files_dir + '/' + str(cluster_id) + ".fastq"
            
            if cluster_size > max_spoa_size:
                # sort fastq reads by length and quality
                sorted_file = fastq_files_dir + '/' + str(cluster_id) + ".sorted.fastq"
                print("Sorting fastq for cluster", cluster_id)
                subprocess.call('python3.6 ' + bin_dir + 'sort_fastq_reads.py ' + '--fastq ' + cluster_fastq + ' --outfolder ' + fastq_files_dir + ' --outfile ' + sorted_file + ' --t ' + n_cores, shell=True)
                print("Done sorting fastq for cluster", cluster_id)
                
                # fix read names in the sorted fastq
                fixed_file = fastq_files_dir + '/' + str(cluster_id) + ".sorted.fixed.fastq"
                subprocess.call("sed 's/_runid=/ runid=/' " + sorted_file + " > " + fixed_file, shell=True)
                
                # take top max_spoa_size reads
                top_file = fastq_files_dir + '/' + str(cluster_id) + ".top.fixed.fastq"
                subprocess.call('head -' + str(top_lines) + ' ' + fixed_file + ' > ' + top_file, shell=True)
                
            else:
                # fix read names in fastq
                fixed_file = fastq_files_dir + '/' + str(cluster_id) + ".fixed.fastq"
                subprocess.call("sed 's/_runid=/ runid=/' " + cluster_fastq + " > " + fixed_file, shell=True)
                       
    return np.array(cluster_size_list)

    
if __name__ == "__main__":
    
    cluster_file = sys.argv[1]
    fastq_files_dir = sys.argv[2]
    n_cores = sys.argv[3]
    
    bin_dir = "/mnt/disk39/user/ltung/ML_error_correction/bin/"    
    
    cluster_df = pd.read_csv(cluster_file, sep='\t', engine='c', header=None, low_memory=False)
    
    cluster_size_array = sort_reads_for_clusters(cluster_df, fastq_files_dir, n_cores, bin_dir)
    
    np.savetxt('final_clusters_sizes', cluster_size_array, fmt='%i %i')

