#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

import sys
import os

def breakup_fastq(fastq_input_filepath, output_dir_filepath, partition_size):
    filein_name = os.path.basename(fastq_input_filepath)[:-6]
    filein = open(fastq_input_filepath, "r")
    total_seq_count = 0
    seq_count = 0
    line_count = 0
    file_count = 1
    fileout = open(output_dir_filepath + '/' + filein_name + '_' + str(file_count) + '.fastq', "w")
    for i in filein:
        if seq_count == partition_size:
            seq_count = 0
            fileout.close()
            file_count += 1
            fileout = open(output_dir_filepath + '/' + filein_name + '_' + str(file_count) + '.fastq', "w")
        line_count += 1
        if line_count == 4:
            seq_count += 1
            total_seq_count += 1
            line_count = 0
        fileout.write(i)
    fileout.close()
    filein.close()
    return total_seq_count

def run_foreach_file_in_directory(fastq_input_dirpath, output_dir_filepath, partition_size):
    if fastq_input_dirpath[-1] != '/':
        fastq_input_dirpath += '/'
    if output_dir_filepath[-1] != '/':
        output_dir_filepath += '/'
    if not os.path.exists(output_dir_filepath):
        os.makedirs(output_dir_filepath)
    read_counts = []
    for i in os.listdir(fastq_input_dirpath):
        if i[0] == '.' or i == 'README' or i == os.path.basename(output_dir_filepath[:-1]):
            continue
        fastq_input_filepath = fastq_input_dirpath + i
        total_seq_count = breakup_fastq(fastq_input_filepath, output_dir_filepath, partition_size)
        read_counts.append(total_seq_count)
    lowest_read_count = min(read_counts)
    print lowest_read_count
    return

if __name__ == '__main__':
    #breakup_fastq(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    run_foreach_file_in_directory(sys.argv[1], sys.argv[2], int(sys.argv[3]))
