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
import random
import subprocess
from scipy.misc import comb
import os

def read_in_cdr3_seqs_and_counts(cdr3_seqs_and_counts_filepath):
    """This script simply reads in the input file and save it to memory, as both a list of cdr3_counts (list of integers) and cdr3 sequences (as a list of strings), where the elements of each list correspond to each other."""
    filein = open(cdr3_seqs_and_counts_filepath, "r")
    cdr3_counts = []
    cdr3_seqs = []
    for i in filein:
        if i[0] == '0':
            continue
        line = i[:-1].split('\t')
        cdr3_seq = line[1]
        if cdr3_seq == 'unmapped':
            break
        cdr3_count = int(line[0])
        cdr3_seqs.append(cdr3_seq)
        cdr3_counts.append(cdr3_count)
    filein.close()
    return cdr3_counts, cdr3_seqs

def call_needle_align(temp_cdr3_seq_1_filepath, temp_cdr3_seq_2_filepath, alignment_filepath):
    """This script simply calls a very short bash script that I wrote that runs the needle global alignment program from EMBOSS."""
    subprocess.call(['bash', 'call_needle_align.bash', temp_cdr3_seq_1_filepath, temp_cdr3_seq_2_filepath, alignment_filepath])
    return

def parse_needle_output_and_get_sub_sum(alignment_filepath, cdr3_counts_sub_set):
    """This script simultaneously parses the needle output in order to get the percent distance for each alignment, and scales that value by the count of cdr3 sequences that are of each corresponding type. Returns the sum of these scaled values."""
    filein = open(alignment_filepath, "r")
    sub_sum_of_distances = 0
    index = 0
    for i in filein:
        if i[:11] == '# Identity:':
            percent_dist = 100 - float(i[-7:-3])
            sub_sum_of_distances += cdr3_counts_sub_set[index] * percent_dist
            index += 1
    filein.close()
#    print '----------------------------------'
#    print 'percent_identity_values_found:', index
#    print 'number_that_should_of_been_found:', len(cdr3_counts_sub_set)
#    print '----------------------------------'
    sys.stdout.flush()
    return sub_sum_of_distances

def pairwise_comps(cdr3_counts, cdr3_seqs, temp_output_dirpath):
    """This script is the heart of the module. It runs pairwise comparisons of each sequence in 'cdr3_seqs' by using 'needle' from EMBOSS to align them. It then weights each of these comparisons by takings into account the counts (from 'cdr3_counts') for each sequence. It does this according to a standard equation for pi (a measure of diversity), which is basically the same as mean pairwise genetic distance. It returns 'mean_percent_distance'."""
    num_unique_cdr3_seqs = len(cdr3_counts)
    total_num_cdr3_seqs = sum(cdr3_counts)
    sum_of_distances = 0
    for i in xrange(num_unique_cdr3_seqs-1):
        sub_sum_of_distances = 0
        temp_file_random_suffix = str(random.random())
        temp_cdr3_seq_1_filepath = temp_output_dirpath + '1_' + temp_file_random_suffix + '.fasta'
        fileout = open(temp_cdr3_seq_1_filepath, "w")
        fileout.write('>' + str(i) + '\n' + cdr3_seqs[i] + '\n')
        fileout.close()
        temp_cdr3_seq_2_filepath = temp_output_dirpath + '2_' + temp_file_random_suffix + '.fasta'
        fileout = open(temp_cdr3_seq_2_filepath, "w")
        for j in xrange(i+1, num_unique_cdr3_seqs):
            #I think needle runs faster if the '-bsequence' is a fasta file with many sequences in it.
            #This is why I save the sequences to disk in such a way. It makes for kind of wonky code,
            #as I have to loop through each 'bsequence' to write them to disk, and then loop through 
            #each entry in the needle output to get the percent distance. But I think it ends up
            #being faster.
            fileout.write('>' + str(j) + '\n' + cdr3_seqs[j] + '\n')
        fileout.close()
        alignment_filepath = temp_output_dirpath + 'alignment_' + temp_file_random_suffix + '.needle'
        call_needle_align(temp_cdr3_seq_1_filepath, temp_cdr3_seq_2_filepath, alignment_filepath)
        sub_sum_of_distances = parse_needle_output_and_get_sub_sum(alignment_filepath, cdr3_counts[i+1:])
        subprocess.call(['rm', temp_cdr3_seq_1_filepath, temp_cdr3_seq_2_filepath, alignment_filepath])
        sum_of_distances += cdr3_counts[i] * sub_sum_of_distances
    mean_percent_distance = sum_of_distances / round(comb(total_num_cdr3_seqs, 2))
    return mean_percent_distance

def calc_diversity(cdr3_seqs_and_counts_filepath, temp_output_dirpath):
    """This module is similar to 'calc_pairwise_dist' except it uses a better way of calculating pi that takes advantage of having multiple cdr3 sequences that are identical (ie having multiple cdr3 seqs of the same class). Because of this it involves much less alignments to be made using 'needle' from EMBOSS, and it also uses the 'cdr3_seq_and_count' filetype. Which may be more reliable (the 'cdr3_seqs' filetype may have been corrupted when the script that made it was writing the output to disk). This script in particular just runs the pipeline, and returns the mean percent difference ('mean_percent_distance'). 'cdr3_seqs_and_counts_filepath' = input filepath to the filetype that has all the unique cdr3 sequences with their counts as well."""
    cdr3_counts, cdr3_seqs = read_in_cdr3_seqs_and_counts(cdr3_seqs_and_counts_filepath)
    if not os.path.exists(temp_output_dirpath):
        os.makedirs(temp_output_dirpath)
    mean_percent_distance = pairwise_comps(cdr3_counts, cdr3_seqs, temp_output_dirpath)
    subprocess.call(['rm', '-r', temp_output_dirpath])
    print 'For file:', cdr3_seqs_and_counts_filepath
    print 'mean percent distance is:', mean_percent_distance
    sys.stdout.flush()
    return mean_percent_distance

def run(cdr3_seqs_master_dirpath, output_dirpath):
    """This script's funciton is to call the module that calculates the mean pairwise distance when given a set of cdr3 sequences. Originally it was designed to call 'calc_pairwise_distance.py', but this script was not functioning very well and I found a major inefficieny with it as well. So, I then made 'calc_pairwise_dist_improved.py' to correct this inefficeincy (see script for details). Other than that, this script is pretty simple, it just formats the input and output filepaths, gets the right inputs and outputs corresponding the provided SGE_TASK_ID environment variable, and send it all to 'calc_pairwise_dist_improved'. 'cdr3_seqs_master_dirpath' = the input directory path that contiains each timepoint's cdr3 sequence and counts filetypes, in addition each timepoint has many trials of these cdr3 seq and count filetypes, as they were created by simulating the removal of reads, which requires multiple trials; output_dirpath = directory path of which to save the output for each timepoint, consisting of mean pairwise genetic distance for each trial, with the overall average over all trials."""
    if cdr3_seqs_master_dirpath[-1] != '/':
        cdr3_seqs_master_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    cdr3_seqs_tpoint_dirpaths = [cdr3_seqs_master_dirpath + i + '/' for i in os.listdir(cdr3_seqs_master_dirpath) if i[0] != '.']
    cdr3_seqs_tpoint_dirpath = cdr3_seqs_tpoint_dirpaths[sge_task_id - 1]
    output_filepath = output_dirpath + os.path.basename(cdr3_seqs_tpoint_dirpath[:-1])
    fileout = open(output_filepath, "w")
    fileout.write('trial\tmean_pairwise_distance\n')
    mean_percent_dist_sum = 0
    count = 0
    for i in os.listdir(cdr3_seqs_tpoint_dirpath):
        if i[0] == '.':
            continue
        cdr3_seqs_filepath = cdr3_seqs_tpoint_dirpath + i
        temp_output_dirpath = output_dirpath + 'temp_alignment_output_' + str(sge_task_id) + '_' + i + '/'
        mean_percent_dist = calc_diversity(cdr3_seqs_filepath, temp_output_dirpath)
        mean_percent_dist_sum += mean_percent_dist
        count += 1
        fileout.write(i + '\t' + str(mean_percent_dist) + '\n')
        print 'File:', cdr3_seqs_filepath
        print 'done'
        sys.stdout.flush()
    fileout.write('avrg\t' + str(mean_percent_dist_sum / count) + '\n')
    fileout.close()
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
