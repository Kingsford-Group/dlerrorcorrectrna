# process_msa_spoa.py
#
# Name: Laura Tung
#
# Usage: python process_msa_spoa.py <msa_spoa_file> <reads_fastq_file> <read_id>
#
# <read_id>: "set": create features for all reads in <msa_spoa_file>, OR <read_id>: create features for <read_id> only.


#import pdb; pdb.set_trace() # Uncomment to debug code using pdb (like gdb)

import sys
import numpy as np
import pandas as pd
import gzip, time


start = time.time()


def construct_features_spoa_set(msa_df, reads_quality_df, new_reads_quality_df, consensus_seq, consensus_occur_freq, consensus_avg_qual, strands_df):
    
    data_list = []
    
    num_spoa_reads = len(msa_df)
    
    for key in msa_df:
        ref_read_seq = msa_df[key]
        ref_read_qual = new_reads_quality_df[key]
        read_length = len(reads_quality_df[key])
        strand = strands_df[key]
        
        new_ref_read_seq, new_ref_read_qual, new_consensus_seq, new_consensus_occur_freq, new_consensus_avg_qual, GC_content = construct_features_per_read(ref_read_seq, consensus_seq, ref_read_qual, consensus_occur_freq, consensus_avg_qual, read_length)
        
        data_list.append([key, new_ref_read_seq, new_ref_read_qual, new_consensus_seq, new_consensus_occur_freq, new_consensus_avg_qual, read_length, GC_content, num_spoa_reads, strand])
    
    data_df = pd.DataFrame(data_list, columns = ["read_name", "read_sequence", "read_quality", "consensus_sequence", "consensus_occur_frequency", "consensus_avg_quality", "read_length", "GC_content", "number_spoa_reads", "strand"])
    
    return data_df


def construct_features_single_read(msa_df, reads_quality_df, new_reads_quality_df, consensus_seq, consensus_occur_freq, consensus_avg_qual, strands_df, read_id):
    
    num_spoa_reads = len(msa_df)
    
    key = read_id
    ref_read_seq = msa_df[key]
    ref_read_qual = new_reads_quality_df[key]
    read_length = len(reads_quality_df[key])
    strand = strands_df[key]
    
    new_ref_read_seq, new_ref_read_qual, new_consensus_seq, new_consensus_occur_freq, new_consensus_avg_qual, GC_content = construct_features_per_read(ref_read_seq, consensus_seq, ref_read_qual, consensus_occur_freq, consensus_avg_qual, read_length)
    
    data_list = [[key, new_ref_read_seq, new_ref_read_qual, new_consensus_seq, new_consensus_occur_freq, new_consensus_avg_qual, read_length, GC_content, num_spoa_reads, strand]]
    
    data_df = pd.DataFrame(data_list, columns = ["read_name", "read_sequence", "read_quality", "consensus_sequence", "consensus_occur_frequency", "consensus_avg_quality", "read_length", "GC_content", "number_spoa_reads", "strand"])
    
    return data_df
        

def construct_features_per_read(ref_read_seq, consensus_seq, ref_read_qual, consensus_occur_freq, consensus_avg_qual, read_length):
    
    # ignore leading and ending "-" in ref_read_seq
    start_idx, end_idx = find_start_end_index(ref_read_seq)
    
    new_ref_read_seq = []
    new_consensus_seq = []    
    new_ref_read_qual = []
    
    new_consensus_occur_freq = []
    new_consensus_avg_qual = []
        
    GC_count = 0
    
    # skip common '-' in both ref_read_seq and consensus_seq
    for i in range(start_idx, end_idx):
        if (ref_read_seq[i] == '-') and (consensus_seq[i] == '-'):
            continue
        else:
            new_ref_read_seq.append(ref_read_seq[i])
            new_consensus_seq.append(consensus_seq[i])
            new_ref_read_qual.append(ref_read_qual[i])
            
            new_consensus_occur_freq.append(consensus_occur_freq[i])
            new_consensus_avg_qual.append(consensus_avg_qual[i])
            
            if (ref_read_seq[i] == 'G') or (ref_read_seq[i] == 'C'):
                GC_count += 1
                
    GC_content = GC_count/read_length
             
    return new_ref_read_seq, new_ref_read_qual, new_consensus_seq, new_consensus_occur_freq, new_consensus_avg_qual, GC_content


def compute_consensus_freq_qual(msa_seq_array, msa_quality_array, consensus_seq):
    
    consensus_occur_freq = []
    consensus_avg_qual = []
    
    spoa_set_size = msa_seq_array.shape[0]
    avg_qual_sum = 0
    letter_col_count = 0
    
    for i in range(len(consensus_seq)):            
        occurences = [n for n in range(spoa_set_size) if msa_seq_array[n,i] == consensus_seq[i]]
        freq = len(occurences)/spoa_set_size
        consensus_occur_freq.append(freq)
        
        if consensus_seq[i] != '-':
            occur_qual = msa_quality_array[occurences,i]
            avg_qual = np.mean(occur_qual)
            avg_qual_sum += avg_qual
            letter_col_count += 1
        else:
            avg_qual = 0
        consensus_avg_qual.append(avg_qual)
    
    # for '-' in consensus, use the mean of average quality of letter columns.        
    mean_avg_qual = avg_qual_sum/letter_col_count
    modified_consensus_avg_qual = [mean_avg_qual if x == 0 else x for x in consensus_avg_qual]
    
    return consensus_occur_freq, modified_consensus_avg_qual


def find_start_end_index(ref_read):
    
    # ignore leading and ending "-" in the ref_read
    cut_leading = ref_read.lstrip('-')
    start_idx = len(ref_read)-len(cut_leading)
    
    cut_ending = ref_read.rstrip('-')
    end_idx = len(cut_ending)    # the end index is not at the last letter; it is "the index of last letter"+1.
    
    return start_idx, end_idx


def correspond_quality_with_sequence(msa_df, reads_quality_df):
    
    msa_seq_list = []
    msa_quality_list = []
    new_reads_quality_df = {}
    
    for key in msa_df:
        read_msa_seq = msa_df[key]
        read_qual_scores = reads_quality_df[key]
        
        # insert 0 quality score for '-' bases
        new_read_qual_scores = np.zeros(len(read_msa_seq), dtype=int)
        j = 0
        for i in range(len(read_msa_seq)):
            if read_msa_seq[i] != '-':
                new_read_qual_scores[i] = read_qual_scores[j]
                j += 1
                
        msa_seq_list.append(list(read_msa_seq))
        msa_quality_list.append(new_read_qual_scores)
        
        new_reads_quality_df[key] = new_read_qual_scores
        
    return np.array(msa_seq_list), np.array(msa_quality_list), new_reads_quality_df
                

def read_msa_spoa_matrix(msa_spoa_file):     # reading MSA SPOA output fasta file (MSA SPOA matrix)
    
    msa_file = open(msa_spoa_file, 'r')
    
    msa_df, num, current = {}, -1, ''
    read_count = 0
    
    for line in msa_file:
        line = line.strip()
        num = (num + 1) % 2
        if num == 0:
            current = line[1:]
            read_count += 1          
        else:
            if current == "Consensus":
                consensus_seq = line
            else:
                msa_df[current] = line
    
    msa_file.close()
    
    return msa_df, consensus_seq


def get_reads_fastq_quality(reads_fastq_file):  # reading fastq file, creating a dictionary of mapping read names to quality scores

    echo('Reading reads fastq file...')
    if reads_fastq_file[-3:] == '.gz':
        reads_file = gzip.open(reads_fastq_file, 'r')
    else:
        reads_file = open(reads_fastq_file, 'r')
        
    reads_df, num, current = {}, -1, ''
    read_count = 0
    strands_df = {}
    
    for line in reads_file:
        if reads_fastq_file[-3:] == '.gz':
            line = line.decode('utf8').strip()
        else:
            line = line.strip()
        num = (num + 1) % 4
        if num == 0:
            current = line[1:].split(' ')[0]
            read_count += 1
            
            strand = line.split('strand=')[-1][0]
            strands_df[current] = strand

        elif num == 3:
            #scores = [int((ord(ch)-33)*2.75) for ch in line]
            scores = [(ord(ch)-33) for ch in line]     # Phred+33: 0 to 93, using ASCII 33 to 126.
            #scores = [ord(ch) for ch in line]
            reads_df[current] = scores

    reads_file.close()
    echo('Finished reading reads fastq file.')

    return reads_df, strands_df


def echo(msg, prepend='', postpend=''):
    global start
    seconds = time.time() - start
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    hms = "%02d:%02d:%02d" % (h, m, s)
    print(prepend + '[' + hms + '] ' + msg + postpend)
    

if __name__ == "__main__":

    msa_spoa_file = sys.argv[1]
    reads_fastq_file = sys.argv[2]
    read_id = sys.argv[3]
    
    reads_quality_df, strands_df = get_reads_fastq_quality(reads_fastq_file)
    msa_df, consensus_seq = read_msa_spoa_matrix(msa_spoa_file)
    
    msa_seq_array, msa_quality_array, new_reads_quality_df = correspond_quality_with_sequence(msa_df, reads_quality_df)
    
    print("msa_seq_array shape:", msa_seq_array.shape)
    print("msa_quality_array shape:", msa_quality_array.shape)
    print("new_reads_quality_df length:", len(new_reads_quality_df))
     
    consensus_occur_freq, consensus_avg_qual = compute_consensus_freq_qual(msa_seq_array, msa_quality_array, consensus_seq)

    print("consensus_occur_freq length:", len(consensus_occur_freq))
    print("consensus_avg_qual length:", len(consensus_avg_qual))
    
    if read_id == "set":
        data_df = construct_features_spoa_set(msa_df, reads_quality_df, new_reads_quality_df, consensus_seq, consensus_occur_freq, consensus_avg_qual, strands_df)
    else:
        data_df = construct_features_single_read(msa_df, reads_quality_df, new_reads_quality_df, consensus_seq, consensus_occur_freq, consensus_avg_qual, strands_df, read_id)

    #print("data_df:")
    #print(data_df.info())
    
    pd.DataFrame.to_csv(data_df, path_or_buf=msa_spoa_file[:-3] + "_features.csv", index=False)
    
    
