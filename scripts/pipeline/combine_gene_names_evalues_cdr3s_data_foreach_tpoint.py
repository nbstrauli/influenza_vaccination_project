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

def pool_cdr3s(map_reads_output_dirpaths, combined_cdr3_seqs_output_filepath):
    """This script """
    fileout = open(combined_cdr3_seqs_output_filepath, "w")
    for i in map_reads_output_dirpaths:
        cdr3_seqs_filepath = i + '/unmapped_igblast_output_cdr3_seq_filtered'
        filein = open(cdr3_seqs_filepath, "r")
        for j in filein:
            if j[:-1].split('\t')[1] == 'Null':
                continue
            else:
                fileout.write(j)
        filein.close()
    fileout.close()
    return

def add_to_file(combined_fileout, partitioned_filein, count):
    """THis script does the actual concatenation files. The count variable is recording the line number for the total number of evalues/gene_names that have been recorded. 'combined_fileout' = the file that is the combined data for all the files (the output file). 'partitioned_filein' = the input files that are partitions of the data from a single time-point."""
    for i in partitioned_filein:
        count += 1
        combined_fileout.write(str(count) + '\t' + i.split('\t')[1])
    partitioned_filein.close()
    return count

def combine(map_reads_output_dirpaths, combined_gene_names_output_dirpath, combined_evalues_output_dirpath):
    """This script is actually super simple but long. It is long because each of the different gene classes needs its own variable assigned for each of the steps of this script. This script cycles through the output of the map_reads.bash script to find the particular files that have the gene names and evalues of each of the unmapped fastq reads. It takes these files and basically combines them together by a simple concatenation. 'map_reads_output_dirpath' = input directory. 'combined_gene_names_output_dirpath' = output directory for each of the gene classes for the gene names data. 'combined_evalues_output_dirpath' = output directory for each of the gene classes for the evalues data."""
    vgene_heavy_name_fileout = open(combined_gene_names_output_dirpath + 'vgene_heavy_name', "w")
    dgene_heavy_name_fileout = open(combined_gene_names_output_dirpath + 'dgene_heavy_name', "w")
    jgene_heavy_name_fileout = open(combined_gene_names_output_dirpath + 'jgene_heavy_name', "w")
    vgene_lambda_name_fileout = open(combined_gene_names_output_dirpath + 'vgene_lambda_name', "w")
    jgene_lambda_name_fileout = open(combined_gene_names_output_dirpath + 'jgene_lambda_name', "w")
    vgene_kappa_name_fileout = open(combined_gene_names_output_dirpath + 'vgene_kappa_name', "w")
    jgene_kappa_name_fileout = open(combined_gene_names_output_dirpath + 'jgene_kappa_name', "w")
    vgene_heavy_evalue_fileout = open(combined_evalues_output_dirpath + 'vgene_heavy_evalue', "w")
    dgene_heavy_evalue_fileout = open(combined_evalues_output_dirpath + 'dgene_heavy_evalue', "w")
    jgene_heavy_evalue_fileout = open(combined_evalues_output_dirpath + 'jgene_heavy_evalue', "w")
    vgene_lambda_evalue_fileout = open(combined_evalues_output_dirpath + 'vgene_lambda_evalue', "w")
    jgene_lambda_evalue_fileout = open(combined_evalues_output_dirpath + 'jgene_lambda_evalue', "w")
    vgene_kappa_evalue_fileout = open(combined_evalues_output_dirpath + 'vgene_kappa_evalue', "w")
    jgene_kappa_evalue_fileout = open(combined_evalues_output_dirpath + 'jgene_kappa_evalue', "w")
    count_vgene_heavy_name = 0
    count_dgene_heavy_name = 0
    count_jgene_heavy_name = 0
    count_vgene_lambda_name = 0
    count_jgene_lambda_name = 0
    count_vgene_kappa_name = 0
    count_jgene_kappa_name = 0
    count_vgene_heavy_evalue = 0
    count_dgene_heavy_evalue = 0
    count_jgene_heavy_evalue = 0
    count_vgene_lambda_evalue = 0
    count_jgene_lambda_evalue = 0
    count_vgene_kappa_evalue = 0
    count_jgene_kappa_evalue = 0
    for i in map_reads_output_dirpaths:
        vgene_heavy_name_filein = open(i + '/unmapped_igblast_output_vgene_heavy_name', "r")
        count_vgene_heavy_name = add_to_file(vgene_heavy_name_fileout, vgene_heavy_name_filein, count_vgene_heavy_name)
        dgene_heavy_name_filein = open(i + '/unmapped_igblast_output_dgene_heavy_name', "r")
        count_dgene_heavy_name = add_to_file(dgene_heavy_name_fileout, dgene_heavy_name_filein, count_dgene_heavy_name)
        jgene_heavy_name_filein = open(i + '/unmapped_igblast_output_jgene_heavy_name', "r")
        count_jgene_heavy_name = add_to_file(jgene_heavy_name_fileout, jgene_heavy_name_filein, count_jgene_heavy_name)
        vgene_lambda_name_filein = open(i + '/unmapped_igblast_output_vgene_lambda_name', "r")
        count_vgene_lambda_name = add_to_file(vgene_lambda_name_fileout, vgene_lambda_name_filein, count_vgene_lambda_name)
        jgene_lambda_name_filein = open(i + '/unmapped_igblast_output_jgene_lambda_name', "r")
        count_jgene_lambda_name = add_to_file(jgene_lambda_name_fileout, jgene_lambda_name_filein, count_jgene_lambda_name)
        vgene_kappa_name_filein = open(i + '/unmapped_igblast_output_vgene_kappa_name', "r")
        count_vgene_kappa_name = add_to_file(vgene_kappa_name_fileout, vgene_kappa_name_filein, count_vgene_kappa_name)
        jgene_kappa_name_filein = open(i + '/unmapped_igblast_output_jgene_kappa_name', "r")
        count_jgene_kappa_name = add_to_file(jgene_kappa_name_fileout, jgene_kappa_name_filein, count_jgene_kappa_name)
        vgene_heavy_evalue_filein = open(i + '/unmapped_igblast_output_vgene_heavy_evalue', "r")
        count_vgene_heavy_evalue = add_to_file(vgene_heavy_evalue_fileout, vgene_heavy_evalue_filein, count_vgene_heavy_evalue)
        dgene_heavy_evalue_filein = open(i + '/unmapped_igblast_output_dgene_heavy_evalue', "r")
        count_dgene_heavy_evalue = add_to_file(dgene_heavy_evalue_fileout, dgene_heavy_evalue_filein, count_dgene_heavy_evalue)
        jgene_heavy_evalue_filein = open(i + '/unmapped_igblast_output_jgene_heavy_evalue', "r")
        count_jgene_heavy_evalue = add_to_file(jgene_heavy_evalue_fileout, jgene_heavy_evalue_filein, count_jgene_heavy_evalue)
        vgene_lambda_evalue_filein= open(i + '/unmapped_igblast_output_vgene_lambda_evalue', "r")
        count_vgene_lambda_evalue = add_to_file(vgene_lambda_evalue_fileout, vgene_lambda_evalue_filein, count_vgene_lambda_evalue)
        jgene_lambda_evalue_filein= open(i + '/unmapped_igblast_output_jgene_lambda_evalue', "r")
        count_jgene_lambda_evalue = add_to_file(jgene_lambda_evalue_fileout, jgene_lambda_evalue_filein, count_jgene_lambda_evalue)
        vgene_kappa_evalue_filein = open(i + '/unmapped_igblast_output_vgene_kappa_evalue', "r")
        count_vgene_kappa_evalue = add_to_file(vgene_kappa_evalue_fileout, vgene_kappa_evalue_filein, count_vgene_kappa_evalue)
        jgene_kappa_evalue_filein = open(i + '/unmapped_igblast_output_jgene_kappa_evalue', "r")
        count_jgene_kappa_evalue = add_to_file(jgene_kappa_evalue_fileout, jgene_kappa_evalue_filein, count_jgene_kappa_evalue)
    vgene_heavy_name_fileout.close()
    dgene_heavy_name_fileout.close()
    jgene_heavy_name_fileout.close()
    vgene_lambda_name_fileout.close()
    jgene_lambda_name_fileout.close()
    vgene_kappa_name_fileout.close()
    jgene_kappa_name_fileout.close()
    vgene_heavy_evalue_fileout.close()
    dgene_heavy_evalue_fileout.close()
    jgene_heavy_evalue_fileout.close()
    vgene_lambda_evalue_fileout.close()
    jgene_lambda_evalue_fileout.close()
    vgene_kappa_evalue_fileout.close()
    jgene_kappa_evalue_fileout.close()
    print count_vgene_heavy_name
    print count_dgene_heavy_name
    print count_jgene_heavy_name
    print count_vgene_lambda_name
    print count_jgene_lambda_name
    print count_vgene_kappa_name
    print count_jgene_kappa_name
    print count_vgene_heavy_evalue
    print count_dgene_heavy_evalue
    print count_jgene_heavy_evalue
    print count_vgene_lambda_evalue
    print count_jgene_lambda_evalue
    print count_vgene_kappa_evalue
    print count_jgene_kappa_evalue
    return

def run_foreach_tpoint(map_reads_output_master_dirpath, combined_gene_names_output_master_dirpath, combined_evalues_output_master_dirpath, combined_cdr3_seqs_output_dirpath):
    """This module takes the output from the 'map_reads.bash' script that consist of the evalues and gene names of each of the unmapped fastq reads for the same influenza vaccination time-point, and combines them into one file. This script does this for each of the time-points. 'map_reads_output_master_dirpath' = the directory that contains the output from the 'map_reads.bash' script. 'combined_gene_names_output_master_dirpath' = the master directory that the combined gene names data for each of the time-points will be written to. 'combined_evalues_output_master_dirpath' = This is the master directory that the combined evalues data will be written to for each of the time-points."""
    if map_reads_output_master_dirpath[-1] != '/':
        map_reads_output_master_dirpath += '/'
    if combined_gene_names_output_master_dirpath[-1] != '/':
        combined_gene_names_output_master_dirpath += '/'
    if not os.path.exists(combined_gene_names_output_master_dirpath):
        os.makedirs(combined_gene_names_output_master_dirpath)
    if combined_evalues_output_master_dirpath[-1] != '/':
        combined_evalues_output_master_dirpath += '/'
    if not os.path.exists(combined_evalues_output_master_dirpath):
        os.makedirs(combined_evalues_output_master_dirpath)
    if combined_cdr3_seqs_output_dirpath[-1] != '/':
        combined_cdr3_seqs_output_dirpath += '/'
    if not os.path.exists(combined_cdr3_seqs_output_dirpath):
        os.makedirs(combined_cdr3_seqs_output_dirpath)

    #this is a dictionary where each entry is a time point, and
    #each definition is a list of the map_reads output directories
    #that correspond to that timepoint
    tpoint_dirpaths = {}
    for i in os.listdir(map_reads_output_master_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        tpoint = i.split('_')[-2]
        try:
            tpoint_dirpaths[tpoint].append(map_reads_output_master_dirpath + i + '/')
        except KeyError:
            tpoint_dirpaths[tpoint] = [map_reads_output_master_dirpath + i + '/']

    for i in tpoint_dirpaths:
        combined_gene_names_output_dirpath = combined_gene_names_output_master_dirpath + i + '/'
        if not os.path.exists(combined_gene_names_output_dirpath):
            os.makedirs(combined_gene_names_output_dirpath)
        combined_evalues_output_dirpath = combined_evalues_output_master_dirpath + i + '/'
        if not os.path.exists(combined_evalues_output_dirpath):
            os.makedirs(combined_evalues_output_dirpath)
        combine(tpoint_dirpaths[i], combined_gene_names_output_dirpath, combined_evalues_output_dirpath)
        combined_cdr3_seqs_output_filepath = combined_cdr3_seqs_output_dirpath + i
        pool_cdr3s(tpoint_dirpaths[i], combined_cdr3_seqs_output_filepath)
    return

if __name__ == '__main__':
    run_foreach_tpoint(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
