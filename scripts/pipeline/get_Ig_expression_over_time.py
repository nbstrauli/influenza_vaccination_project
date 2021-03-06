#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=2G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00


import sys
import os
import itertools
import subprocess

def count_Ig_mapped_reads(evalues_filepath, evalue_cutoff):
    """This script counts the number of reads that mapped to an Ig locus by way of counting the number of evalues (generated from IgBLAST alignmnets) that are equal to or below 'evalue_cutoff'. It does this by going through each of the evalues that are contained in the file called 'unmapped_igblast_output_evalues' in the map_reads.bash output directories."""
    filein = open(evalues_filepath, "r")
    num_Ig_mapped_reads = 0
    for i in filein:
        evalue = i[:-1].split('\t')[1]
        if evalue == 'Null':
            continue
        if float(evalue) <= evalue_cutoff:
            num_Ig_mapped_reads += 1
    filein.close()
    return num_Ig_mapped_reads

def call_samtools(mapped_reads_BAM_filepath):
    """This script determines the number of mapped reads in the provided BAM file. This BAM file is the 'accepted_hits.bam' file that was generated by tophat for each of the output directories for each of the time-points ('mapped_reads_BAM_filepath'). It does this by usig Samtools."""
    numMappedReads_output_filepath = mapped_reads_BAM_filepath[:-4] + '_numUniquelyMappedReads'
    subprocess.call(['bash', 'call_samtools_count_uniquely_mapped.bash', mapped_reads_BAM_filepath, numMappedReads_output_filepath])
    filein = open(numMappedReads_output_filepath, "r")
    num_mapped_reads = int(filein.readline())
    filein.close()
    return num_mapped_reads

def write_output(mapped_reads_overtime, Ig_mapped_reads_overtime, tpoints, Ig_expression_output_filepath):
    """This script simply parses the output from 'get_num_of_Ig_mapped_reads' and 'get_total_number_of_mapped_reads' to get the expression level of Ig for each timepoint. It then writes all of this data to 'Ig_expression_output_filepath'."""
    fileout = open(Ig_expression_output_filepath, "w")
    fileout.write('timepoint_id\ttotal_number_of_mapped_reads\tnumber_of_Ig_mapped_reads\tratio_of_total_mapped_reads_to_Ig_mapped_reads\n')
    for i, j, k in itertools.izip(mapped_reads_overtime, Ig_mapped_reads_overtime, tpoints):
        fileout.write(str(k) + '\t' + str(i) + '\t' + str(j) + '\t' + str(float(j) / float(i)) + '\n')
    fileout.close()
    return

def get_counts(map_reads_output_master_dirpath, Ig_expression_output_filepath, evalue_cutoff):
    """This script simply calls 'get_total_number_of_mapped_reads.py' and 'get_num_of_Ig_mapped_reads' and retrievs their output. The purpose of this small module is to combine these two outputs to get the expression level of Ig over time. 'map_reads_output_master_dirpath' = the output directory (for a given patient and sample type) of map_reads.bash. 'Ig_expression_output_filepath' = path for the output data to be written to. 'evalue_cutoff' = the cutoff for evalues that is passed to 'get_num_of_Ig_mapped_reads'."""
    """UPDATE: This script was modified so that it instructs 'get_total_number_of_mapped_reads' and 'get_num_of_Ig_mapped_reads' to actually write the output of those scripts to an output file. Because of this, the input for this module had to include the output filepaths of those two outputs ('mapped_reads_output_filepath', and 'Ig_mapped_reads_output_filepath' respectively). This was changed in 'run_pipeline.bash' as well. If one wants to make this script not write output for total mapped reads and Ig mapped reads, then simply make 'mapped_reads_output_filepath' and 'Ig_mapped_reads_output_filepath' equal to 'None'"""
    if map_reads_output_master_dirpath[-1] != '/':
        map_reads_output_master_dirpath += '/'

    #this is a dictionary where each entry is a time point, and
    #each definition is a list of the map_reads output directories
    #that correspond to that timepoint
    tpoint_dirpaths = {}
    for i in os.listdir(map_reads_output_master_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        tpoint = int(i.split('_')[-2])
        try:
            tpoint_dirpaths[tpoint].append(map_reads_output_master_dirpath + i + '/')
        except KeyError:
            tpoint_dirpaths[tpoint] = [map_reads_output_master_dirpath + i + '/']

    mapped_reads_overtime = []
    Ig_mapped_reads_overtime = []
    tpoints = []
    for i in sorted(tpoint_dirpaths):
        num_mapped_reads = 0
        num_Ig_mapped_reads = 0
        for j in tpoint_dirpaths[i]:
            mapped_reads_BAM_filepath = j + 'accepted_hits.bam'
            num_mapped_reads += call_samtools(mapped_reads_BAM_filepath)
            evalues_filepath = j + '/unmapped_igblast_output_evalues'
            num_Ig_mapped_reads += count_Ig_mapped_reads(evalues_filepath, evalue_cutoff)
        mapped_reads_overtime.append(num_mapped_reads)
        Ig_mapped_reads_overtime.append(num_Ig_mapped_reads)
        tpoints.append(i)

    write_output(mapped_reads_overtime, Ig_mapped_reads_overtime, tpoints, Ig_expression_output_filepath)
    return

if __name__ == '__main__':
    get_counts(sys.argv[1], sys.argv[2], float(sys.argv[3]))
