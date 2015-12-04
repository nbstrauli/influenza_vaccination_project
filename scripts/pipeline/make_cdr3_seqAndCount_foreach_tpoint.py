#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

import sys
import os

def make_seq_and_count(pooled_cdr3_seqs_dirpath, seq_and_count_output_dirpath):
    """This script takes the pooled CDR3 sequeces from each time-point ('pooled_cdr3_seqs_dirpath') and formats them into the 'seq_and_count' format of storing CDR3 sequences ('seq_and_count_output_dirpath'). That is the format where each unique CDR3 seq has one line, and this line begins with its count, followed by its actual sequence. I wrote pretty much this same script for the Snyderome data, but it was a bit wonky, so I just rewrote it here for the new influenza vacc. dataset."""
    if pooled_cdr3_seqs_dirpath[-1] != '/':
        pooled_cdr3_seqs_dirpath += '/'
    if seq_and_count_output_dirpath[-1] != '/':
        seq_and_count_output_dirpath += '/'
    if not os.path.exists(seq_and_count_output_dirpath):
        os.makedirs(seq_and_count_output_dirpath)
    for i in os.listdir(pooled_cdr3_seqs_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        print i
        unique_cdr3s = {}
        filein = open(pooled_cdr3_seqs_dirpath + i, "r")
        for j in filein:
            cdr3_seq = j[:-1].split('\t')[1]
            if cdr3_seq == 'N/A' or cdr3_seq == 'Null' or cdr3_seq == 'rejected' or cdr3_seq == '':
                continue
            try:
                unique_cdr3s[cdr3_seq] += 1
            except KeyError:
                unique_cdr3s[cdr3_seq] = 1
        filein.close()
        unique_cdr3s_list = []
        for j in unique_cdr3s:
            unique_cdr3s_list.append([unique_cdr3s[j], j])
        unique_cdr3s_list = sorted(unique_cdr3s_list)
        fileout = open(seq_and_count_output_dirpath + i, "w")
        for j in unique_cdr3s_list:
            fileout.write(str(j[0]) + '\t' + j[1] + '\n')
        fileout.close()
    return

if __name__ == '__main__':
    make_seq_and_count(sys.argv[1], sys.argv[2])
