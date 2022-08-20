# align_create_labels.py
#
# Name: Laura Tung
#
# Usage: python align_create_labels.py <alignment_file_sam> <unique_chimeric_readnames>
#
# <unique_chimeric_readnames>: a file containing chimeric reads names: "unique_chimeric_readnames"

#import pdb; pdb.set_trace() # Uncomment to debug code using pdb (like gdb)

import sys
import numpy as np
import pandas as pd

import re
import pysam
import pickle


def construct_aligned_sequences_labels(read_seq, cs_string, cigar_str):
    
    start_clips = 0
    end_clips = 0
    
    # cs tag does not include the beginning/ending soft clips, extract from CIGAR
    cigar_list = re.findall('(\d+)([MIDSHX=])', cigar_str)
    if cigar_list[0][1] == 'S':
        start_clips = int(cigar_list[0][0])
    if cigar_list[-1][1] == 'S':
        end_clips = int(cigar_list[-1][0])
        
    if start_clips > 0:
        align_read_seq = read_seq[0:start_clips]
        align_true_read_seq = "-"*start_clips
        true_labels = ["Insertion"]*start_clips
    else:
        align_read_seq = ""
        align_true_read_seq = ""
        true_labels = []
        
    # process cs tag
    i = start_clips         # i is the index of read_seq
    for item in re.findall('(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)', cs_string):
        op = item[0]
        if op == ':':       # Identical sequence (match)
            num_matches = int(item[1:])
            match_seg = read_seq[i:i+num_matches]
            
            align_read_seq += match_seg
            align_true_read_seq += match_seg
            true_labels += ["No_correction"]*num_matches
            
            i += num_matches
        
        elif op == '*':     # Substitution: ref to query (mismatch)
            ref_base = item[1].upper()
            query_base = item[2].upper()
            
            align_read_seq += query_base
            
            if ref_base == 'N':
                align_true_read_seq += query_base
                true_labels.append("No_correction")
            else:
                align_true_read_seq += ref_base
                true_labels.append("Substituted_" + ref_base)
            
            i += 1
            
        elif op == '+':     # Insertion to the reference
            insert_seg = item[1:].upper()
            num_insert = len(insert_seg)
            
            align_read_seq += insert_seg
            align_true_read_seq += "-"*num_insert
            true_labels += ["Insertion"]*num_insert
            
            i += num_insert
            
        elif op == '-':     # Deletion from the reference
            delete_seg = item[1:].upper()
            if 'N' in delete_seg:
                delete_seg = delete_seg.replace("N", "")
                
            num_delete = len(delete_seg)
            
            align_read_seq += "-"*num_delete
            align_true_read_seq += delete_seg
            for j in range(num_delete):
                true_labels.append("Deleted_" + delete_seg[j])
    
    # take care of the ending clips    
    if end_clips > 0:
        end_pos = -1*end_clips
        align_read_seq += read_seq[end_pos:]
        align_true_read_seq += "-"*end_clips
        true_labels += ["Insertion"]*end_clips
      
    return align_read_seq, align_true_read_seq, true_labels


def reverse_complement_seq(seq):
    
    mapping = str.maketrans('ATCG-', 'TAGC-')
    
    reverse_complemented_seq = seq.translate(mapping)[::-1]
    
    return reverse_complemented_seq


def reverse_complement_labels(labels):
    
    Label_ReverseComplement = {'Insertion': 'Insertion', 'Substituted_A': 'Substituted_T', \
                               'Substituted_C': 'Substituted_G', 'Substituted_G': 'Substituted_C', \
                               'Substituted_T': 'Substituted_A', 'Deleted_A': 'Deleted_T', \
                               'Deleted_C': 'Deleted_G', 'Deleted_G': 'Deleted_C', \
                               'Deleted_T': 'Deleted_A', 'No_correction': 'No_correction'}
    
    reverse_complemented_labels = [Label_ReverseComplement[nuc] for nuc in labels][::-1]
    
    return reverse_complemented_labels


if __name__ == "__main__":
    
    alignment_file_sam = sys.argv[1]
    unique_chimeric_readnames_file = sys.argv[2]
    
    unique_chimeric_readnames = np.loadtxt(unique_chimeric_readnames_file, dtype='str')
    
    align_labels_dict = {}
    
    in_sam_file = pysam.AlignmentFile(alignment_file_sam, 'r')
    
    for alnm in in_sam_file.fetch(until_eof=True):
        #read_name = alnm.query_name.split("-runid=")[0]
        read_name = alnm.query_name
        if read_name in unique_chimeric_readnames:
            continue
        
        try:
            cs_string = alnm.get_tag('cs')
        except KeyError:
            continue
        cigar_str = alnm.cigarstring
        read_seq = alnm.query_sequence

        align_read_seq, align_true_read_seq, true_labels = construct_aligned_sequences_labels(read_seq, cs_string, cigar_str)
        
        if alnm.is_reverse == True:
            # reverse complemented: convert back to the original sequence
            rev_comp_align_read_seq = reverse_complement_seq(align_read_seq)
            rev_comp_align_true_read_seq = reverse_complement_seq(align_true_read_seq)
            rev_comp_true_labels = reverse_complement_labels(true_labels)
            
            align_labels_dict[read_name] = (rev_comp_align_read_seq, rev_comp_align_true_read_seq, rev_comp_true_labels)
        else:
            align_labels_dict[read_name] = (align_read_seq, align_true_read_seq, true_labels)

    in_sam_file.close()
        
    print("align_labels_dict size:", len(align_labels_dict))
    
    # TESTING
    qname = list(align_labels_dict.keys())[0]
    print(qname)
    value = align_labels_dict[qname]
    print(value[0])
    print(value[1])
    print(value[2])

    qname = list(align_labels_dict.keys())[1]
    print(qname)
    value = align_labels_dict[qname]
    print(value[0])
    print(value[1])
    print(value[2])
    
    qname = list(align_labels_dict.keys())[4]
    print(qname)
    value = align_labels_dict[qname]
    print(value[0])
    print(value[1])
    print(value[2])    
    

    # save the dictionary
    pickle_out = open("align_labels_dict.pickle","wb")
    pickle.dump(align_labels_dict, pickle_out)
    pickle_out.close()
