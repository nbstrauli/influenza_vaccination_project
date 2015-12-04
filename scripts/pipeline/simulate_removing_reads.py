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
import re
import random
import os

def count_reads_fastq(fastq_filepath):
    """Counts the number of reads in a fastq file and returns that number. Terribly simple"""
    filein = open(fastq_filepath, "r")
    line_count = 0
    for i in filein:
        line_count += 1
    read_count = line_count / 4.0
    filein.close()
    return read_count

def make_empirical_pdf(cdr3_seq_and_count_input_filepath, num_of_reads):
    """This script takes the cdr3 seqs and counts file, as well as the total number of reads for the corresponding time-point and turns that information into a discrete empirical distribution. Each element in the distribution corresponds to a unique cdr3 sequence and the value of that element is its frequency in the dataset. In addition, the last element in the pdf corresponds to the unmapped reads in the dataset (ie those reads that did not yeild a cdr3 sequence). Also this script returns the actual counts for each element in the pdf ('cdr3_counts'), and the sequences for each element in the pdf ('cdr3_seqs'), where the last sequence is simply 'unmapped'.
    UPDATE: We use this script in other modules, and I added a functionality to it where it will not add the number of unmapped reads to the end if 'num_or_reads' = False. This functionality is used in the script, 'calc_genetic_dist_between_tpoints.py' for example."""
    cdr3_seq_and_count_regex = '^([0-9]+)\t([^\t]*)\n$'
    filein = open(cdr3_seq_and_count_input_filepath, "r")
    total_num_of_cdr3s = 0
    cdr3_seqs = [] #this is a catalogue of the different unique cdr3 sequences in the data. The last entry in this list is 'unmapped' corresponding to the insidence of unmapped reads in the data
    cdr3_pdf = [] #this will be a vector of the frequencies for each cdr3 sequence in 'cdr3_seqs', where each element in 'cdr3_seqs' corresponds to that of 'cdr3_pdf'
    cdr3_counts = [] #this is the raw counts for each cdr3 seq
    for i in filein:
        if i[:-1].split('\t')[1] == 'unmapped':
            break
        match = re.search(cdr3_seq_and_count_regex, i)
        cdr3_count = float(match.group(1))
        cdr3_seq = match.group(2)
        total_num_of_cdr3s += cdr3_count
        cdr3_seqs.append(cdr3_seq)
        if num_of_reads != False:
            cdr3_pdf.append(cdr3_count / num_of_reads)
        cdr3_counts.append(int(cdr3_count))
    filein.close()
    #added this conditional later
    if num_of_reads == False:
        cdr3_pdf = [float(i) / total_num_of_cdr3s for i in cdr3_counts]
        return cdr3_pdf, cdr3_seqs, cdr3_counts
    num_of_unmapped_reads = num_of_reads - total_num_of_cdr3s
    cdr3_seqs.append('unmapped')
    cdr3_pdf.append(num_of_unmapped_reads / num_of_reads)
    cdr3_counts.append(num_of_unmapped_reads)
    return cdr3_pdf, cdr3_seqs, cdr3_counts

def convert_pdf2cdf(cdr3_pdf):
    """This script takes the pdf and simply converts it to a cdf. Pretty straight forward."""
    total = 0
    cdr3_cdf = []
    for i in cdr3_pdf:
        cdr3_cdf.append(total + i)
        total += i
    return cdr3_cdf

def inverse_transform_sampling(cdr3_cdf):
    """This is a cool general way to sample from any distribution (especially empirical ones). Generates a random # b/t 0 and 1, then finds the largest element in the cdf (x-axis) that has a value (y-axis) that is greater or equal to the random number. Returns the index of the element that it finds ('chosen_cdr3_index')."""
    random_float = random.random()
    for index, i in enumerate(cdr3_cdf):
        if i >= random_float:
            chosen_cdr3_index = index
            break
    return chosen_cdr3_index

def update_variables(chosen_cdr3_index, cdr3_pdf, cdr3_cdf, cdr3_counts, num_of_reads):
    """This script updates all of the variables after a read has been removed. It is running pretty slow, and I may have to rework it if I figure out how to make it faster. The problem is that each element in the distribution needs to be updated, which involves a division of some relatively large numbers. 'chosen_cdr3_index' = the index (for a list) that marks the sequence where 1 has to be removed; 'cdr3_pdf' = the pdf (list of probabilities); 'cdr3_cdf' = list of probs, but cumulative; 'cdr3_counts' = actual counts of each cdr3 and unmapped sequences; 'num_of_reads' = total number of reads currently present in the distributions."""
    num_of_reads -= 1
    total_prob = 0
    for i in xrange(len(cdr3_pdf)):
        if i == chosen_cdr3_index:
            cdr3_counts[i] -= 1
            cdr3_pdf[i] = (cdr3_counts[i] / num_of_reads)
            cdr3_cdf[i] = cdr3_pdf[i] + total_prob
            total_prob += cdr3_pdf[i]
        else:
            cdr3_pdf[i] = cdr3_counts[i] / num_of_reads
            cdr3_cdf[i] = cdr3_pdf[i] + total_prob
            total_prob += cdr3_pdf[i]
    return cdr3_pdf, cdr3_cdf, cdr3_counts, num_of_reads

def write_output(cdr3_seqs, cdr3_counts, cdr3_seq_and_count_output_dirpath, cdr3_seq_output_dirpath, sge_task_id):
    """This script simply writes the output to file, and in the correct format. I makes a cdr3 sequence and count format file, as well as just a cdr3 sequence file (for this type the unique identifier for each sequence is just a count over all sequences, so if there are a total of 10 cdr3 sequences then their ids will be 1-10). For each time-point there are 'trials' number of trials, so it saves output to the directories that pertain to a given time-point (given by 'cdr3_seq_and_count_output_dirpath', and 'cdr3_seq_and_count_output_dirpath') and saves the filename as the trial (given by 'sge_task_id')."""
    if cdr3_seq_and_count_output_dirpath[-1] != '/':
        cdr3_seq_and_count_output_dirpath += '/'
    if not os.path.exists(cdr3_seq_and_count_output_dirpath):
        os.makedirs(cdr3_seq_and_count_output_dirpath)
    if cdr3_seq_output_dirpath[-1] != '/':
        cdr3_seq_output_dirpath += '/'
    if not os.path.exists(cdr3_seq_output_dirpath):
        os.makedirs(cdr3_seq_output_dirpath)
    seq_and_count_fileout = open(cdr3_seq_and_count_output_dirpath + str(sge_task_id), "w")
    seq_fileout = open(cdr3_seq_output_dirpath + str(sge_task_id), "w")
    seq_id = 0
    for i in xrange(len(cdr3_seqs)):
        seq_and_count_fileout.write(str(cdr3_counts[i]) + '\t' + cdr3_seqs[i] + '\n')
        if cdr3_seqs[i] == 'unmapped':
            break
    seq_and_count_fileout.close()
    seq_fileout.close()
    return

def run(cdr3_seq_and_count_input_filepath, timepoint_fastq_filepath, cdr3_seq_and_count_output_dirpath, cdr3_seq_output_dirpath, num_of_reads_to_reduce_to):
    """This script runs the pipeline. It starts by getting the total number of sequencing reads from the timepoint that you are removing reads from, then it makes the pdf, and cdf of the cdr3 seqs and unmapped reads, then call the subscripts to remove reads one by one until it reached the 'num_of_reads_to_reduce_to' value (this value should be 55,832,787 which is the number of reads from the 1st time-point). 'cdr3_seq_and_count_input_filepath' = filepath to the timepoint's cdr3 seqs with counts data; 'timepoint_fastq_dirpath' = directory path to the timepoint's fastq data; 'output_filepath' = filepath to write a cdr3 seqs and counts file to after removing the reads."""
    num_of_reads = count_reads_fastq(timepoint_fastq_filepath)
    sys.stdout.flush()
    cdr3_pdf, cdr3_seqs, cdr3_counts = make_empirical_pdf(cdr3_seq_and_count_input_filepath, num_of_reads)
    cdr3_cdf = convert_pdf2cdf(cdr3_pdf)
    num_of_reads_to_remove = int(num_of_reads - num_of_reads_to_reduce_to)
    for i in xrange(num_of_reads_to_remove):
        chosen_cdr3_index = inverse_transform_sampling(cdr3_cdf)
        cdr3_pdf, cdr3_cdf, cdr3_counts, num_of_reads = update_variables(chosen_cdr3_index, cdr3_pdf, cdr3_cdf, cdr3_counts, num_of_reads)
        if str(i)[-6:] == '000000':
            print 'removed', i, 'reads'
            sys.stdout.flush()
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    write_output(cdr3_seqs, cdr3_counts, cdr3_seq_and_count_output_dirpath, cdr3_seq_output_dirpath, sge_task_id)
    print num_of_reads
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
