#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=2G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=180G
#$ -l h_rt=336:00:00

import sys
import os
import itertools

def make_dic_of_all_gene_names(gene_name_input_filepaths, evalue_input_filepaths, evalue_cutoff, drop_allele_info):
    """This script goes through all of the gene input files (provided by 'gene_name_input_filepaths') as well as the evalue input files (provided by 'evalue_input_filepaths') for all time-points. It does this to create a dictionary of all the gene names that exist (for this gene class), which have a corresponding evalue that is lower than 'evalue_cutoff'. Returns this dictionary, where each index is a gene name, and its definition is the number of counts that this gene name was found (and satisfied the evalue cutoff)."""
    """evalues can only be a float or 'N/A'"""
    gene_names_dic = {}
    for i, j, in itertools.izip(gene_name_input_filepaths, evalue_input_filepaths):
        gene_name_filein = open(i, "r")
        evalue_filein = open(j, "r")
        for k, l in itertools.izip(gene_name_filein, evalue_filein):
            evalue = l.split('\t')[1][:-1]
            if evalue == 'N/A':
                continue
            elif float(evalue) > evalue_cutoff:
                continue
            elif float(evalue) <= evalue_cutoff:
                gene_name = k.split('\t')[1][:-1]
                if gene_name == 'N/A':
                    print 'error: found a defined evalue with and undefined gene name'
                    sys.stdout.flush()
                    return
                if drop_allele_info:
                    #remove the allele information from the gene name
                    gene_name = gene_name.split('*')[0]
                try:
                    gene_names_dic[gene_name] += 1
                except KeyError:
                    gene_names_dic[gene_name] = 1
        gene_name_filein.close()
        evalue_filein.close()
        print '\tgot gene names for:', i
        sys.stdout.flush()
    return gene_names_dic

def make_dic_of_all_cdrs(cdr3_seqs_input_filepaths):
    cdr3_seqs_dic = {}
    for i in cdr3_seqs_input_filepaths:
        filein = open(i, "r")
        for j in filein:
            cdr3_seq = j[:-1].split('\t')[1]
            try:
                cdr3_seqs_dic[cdr3_seq] += 1
            except KeyError:
                cdr3_seqs_dic[cdr3_seq] = 1
        filein.close()
    return cdr3_seqs_dic

def get_gene_usage_foreach_tpoint(gene_name_input_filepaths, evalue_input_filepaths, gene_names_dic, evalue_cutoff, drop_allele_info):
    """This script is very similar to 'make_dic_of_all_gene_names', except instead of getting a gene name dictionary for all time-points, its make an individual dictionary for each time-point. It uses the dictionary for all time-points ('gene_names_dic') to prime the dics for each time-point (which are stored in memory as a list of dics (in chronological order)). It then loops through the gene name and evalue files (provided by 'gene_name_input_filepaths' and 'evalue_input_filepaths', respectively) to get all the gene names from alignments that satisfy the evalue cutoff, for each time-point. Returns a list of dictionaries, where each dic corresponds to a time-point, and the last element in the dictionary-list is all time-points combined."""
    all_tpoints_gene_names_dics = [] #list of dicitonaries for each time-point
    for i, j in itertools.izip(gene_name_input_filepaths, evalue_input_filepaths):
        gene_name_filein = open(i, "r")
        evalue_filein = open(j, "r")
        tpoint_gene_names_dic = {}
        #prime the gene name dictionay for this time-point
        for l in gene_names_dic:
            tpoint_gene_names_dic[l] = 0
        for l, m in itertools.izip(gene_name_filein, evalue_filein):
            evalue = m.split('\t')[1][:-1]
            if evalue == 'N/A':
                continue
            elif float(evalue) > evalue_cutoff:
                continue
            elif float(evalue) <= evalue_cutoff:
                gene_name = l.split('\t')[1][:-1]
                if drop_allele_info:
                    #remove allele information
                    gene_name = gene_name.split('*')[0]
                tpoint_gene_names_dic[gene_name] += 1
        gene_name_filein.close()
        evalue_filein.close()
        all_tpoints_gene_names_dics.append(tpoint_gene_names_dic)
    return all_tpoints_gene_names_dics

def get_cdr3_usage_foreach_tpoint(cdr3_seqs_input_filepaths, cdr3_seqs_dic):
    all_tpoints_cdr3_seqs_dics = [] #list of dicitonaries for each time-point
    for i in cdr3_seqs_input_filepaths:
        filein = open(i, "r")
        tpoint_cdr3_seqs_dic = {}
        #prime the dictionay for this time-point
        for j in cdr3_seqs_dic:
            tpoint_cdr3_seqs_dic[j] = 0
        for j in filein:
            cdr3_seq = j[:-1].split('\t')[1]
            tpoint_cdr3_seqs_dic[cdr3_seq] += 1
        filein.close()
        all_tpoints_cdr3_seqs_dics.append(tpoint_cdr3_seqs_dic)
    return all_tpoints_cdr3_seqs_dics

def get_gene_lengths(ref_seq_dirpath, gene_class, drop_allele_info):
    gene_class_dic = {'vgene_heavy':'IGHV', 'dgene_heavy':'IGHD', 'jgene_heavy':'IGHJ', 'vgene_lambda':'IGLV', 'jgene_lambda':'IGLJ', 'vgene_kappa':'IGKV', 'jgene_kappa':'IGKJ'}
    ref_seq_filepath = '%s%s.fasta' % (ref_seq_dirpath, gene_class_dic[gene_class])
    filein = open(ref_seq_filepath, "r")
    gene_lens_dic = {}
    for i in filein:
        if i[0] == '>':
            gene_name = i[1:-1]
            if drop_allele_info:
                gene_name = gene_name.split('*')[0]
        else:
            length = float(len(i[:-1]))
            try:
                gene_lens_dic[gene_name].append(length)
            except KeyError:
                gene_lens_dic[gene_name] = [length]
    new_gene_lens_dic = {}
    for i in gene_lens_dic:
        if drop_allele_info:
            new_gene_lens_dic[i] = round(sum(gene_lens_dic[i]) / len(gene_lens_dic[i]), 0)
        else:
            new_gene_lens_dic[i] = gene_lens_dic[i][0]
    return new_gene_lens_dic

def get_total_mapped_reads(ig_expression_filepath):
    filein = open(ig_expression_filepath, "r")
    filein.readline()
    total_mapped_reads = []
    for i in filein:
        total_mapped_reads.append(float(i.split('\t')[1]))
    filein.close()
    return total_mapped_reads

def normalize_counts(all_tpoints_gene_names_dics, gene_lens_dic, total_mapped_reads, scaling_factor):
    #for each time-point
    for i in xrange(len(total_mapped_reads)):
        mapped_reads = total_mapped_reads[i]
        #for each gene found in the data
        for j in all_tpoints_gene_names_dics[i]:
            #if this is CDR3 seq data
            if gene_lens_dic == 'cdr3':
                expression_level = (all_tpoints_gene_names_dics[i][j] / mapped_reads) * scaling_factor
            else:
                gene_length = gene_lens_dic[j]
                expression_level = (all_tpoints_gene_names_dics[i][j] / (mapped_reads * gene_length)) * scaling_factor
            all_tpoints_gene_names_dics[i][j] = expression_level
    return all_tpoints_gene_names_dics

def write_output(all_tpoints_gene_names_dics, tpoint_ids, gene_usage_output_dirpath, gene_class):
    #this will be a list if lists, where each element corresponds to a gene,
    #within a gene entry, the 1st element is the range of that gene's
    #expression trajectory, followed by the gene name, then the actual
    #expression trajectory
    gene_expr_trajs = []
    for i in all_tpoints_gene_names_dics[0]:
        gene_name = i
        expr_traj = []
        for j in all_tpoints_gene_names_dics:
            expr_traj.append(j[gene_name])
        range = max(expr_traj) - min(expr_traj)
        gene_expr_trajs.append([range, gene_name, expr_traj])
    #sort according to the range of the expression trajectories
    gene_expr_trajs = sorted(gene_expr_trajs)
    #if this is CDR3 seq data
    if gene_class == 'cdr3':
        fileout = open(gene_usage_output_dirpath, "w")
    else:
        output_filepath = gene_usage_output_dirpath + gene_class
        fileout = open(output_filepath, "w")
    fileout.write('\t' + '\t'.join([str(i) for i in tpoint_ids]) + '\n')
    for i in gene_expr_trajs:
        fileout.write(i[1] + '\t' + '\t'.join([str(j) for j in i[2]]) + '\n')
    fileout.close()
    return

def run(gene_names_master_dirpath, evalues_master_dirpath, gene_usage_output_dirpath, gene_class, evalue_cutoff, drop_allele_info, ref_seq_dirpath, ig_expression_filepath, scaling_factor):
    """This script runs the pipeline. The program gets the frequency of each reference gene (of a given gene class) for each snyderome time-point. It writes these results to 'gene_usage_output_dirpath'."""
    if gene_names_master_dirpath[-1] != '/':
        gene_names_master_dirpath += '/'
    if evalues_master_dirpath[-1] != '/':
        evalues_master_dirpath += '/'
    if gene_usage_output_dirpath[-1] != '/':
        gene_usage_output_dirpath += '/'
    if not os.path.exists(gene_usage_output_dirpath):
        os.makedirs(gene_usage_output_dirpath)
    if ref_seq_dirpath[-1] != '/':
        ref_seq_dirpath += '/'

    #get input gene name and evalue filepaths, and make sure
    #they are in numerical order according to the time points
    tpoint_ids = []
    gene_name_input_filepaths = []
    for i in os.listdir(gene_names_master_dirpath):
        if i[0] == '.':
            continue
        tpoint_ids.append(int(i))
        gene_name_input_filepaths.append([int(i), gene_names_master_dirpath + i + '/' + gene_class + '_name'])
    gene_name_input_filepaths = [i[1] for i in sorted(gene_name_input_filepaths)]
    tpoint_ids = sorted(tpoint_ids)
    evalue_input_filepaths = []
    for i in os.listdir(evalues_master_dirpath):
        if i[0]== '.':
            continue
        evalue_input_filepaths.append([int(i), evalues_master_dirpath + i + '/' + gene_class + '_evalue'])
    evalue_input_filepaths = [i[1] for i in sorted(evalue_input_filepaths)]

    print 'getting dic of all gene names'
    sys.stdout.flush()
    gene_names_dic = make_dic_of_all_gene_names(gene_name_input_filepaths, evalue_input_filepaths, evalue_cutoff, drop_allele_info)
    print 'getting dic for individual time-points'
    sys.stdout.flush()
    all_tpoints_gene_names_dics = get_gene_usage_foreach_tpoint(gene_name_input_filepaths, evalue_input_filepaths, gene_names_dic, evalue_cutoff, drop_allele_info)
    print "normalizing counts"
    gene_lens_dic = get_gene_lengths(ref_seq_dirpath, gene_class, drop_allele_info)
    total_mapped_reads = get_total_mapped_reads(ig_expression_filepath)
    all_tpoints_gene_names_dics = normalize_counts(all_tpoints_gene_names_dics, gene_lens_dic, total_mapped_reads, scaling_factor)
    print "writing output"
    write_output(all_tpoints_gene_names_dics, tpoint_ids, gene_usage_output_dirpath, gene_class)
    return

def run_array(gene_names_master_dirpath, evalues_master_dirpath, gene_usage_output_dirpath, evalue_cutoff, drop_allele_info, ref_seq_dirpath, ig_expression_filepath, scaling_factor, gene_classes):
    """NOTE: This module should be updated to write the order of gene names in the gene_usage_output_files as a seperate output file in its own directory. This was done for 'get_gene_usage_time_progression_no_normalize.py', but didn't get a chance to do it for this script. Mostly because we would have to go back to the influenza pipeline and update that as well. Should do this update if we ever need this module again."""
    """UPDATE: We added the feature to write the gene name order files, but have not updated the pipeline. Must do this soon!"""
    """UPDATE: We added it to the pipeline."""
    if drop_allele_info == 'True':
        drop_allele_info = True
    elif drop_allele_info == 'False':
        drop_allele_info = False
    else:
        print 'The "drop_allele_info" variable must be either True or False'
        print 'It is currently set to:', drop_allele_info
        return
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    gene_class = gene_classes[sge_task_id - 1]
    run(gene_names_master_dirpath, evalues_master_dirpath, gene_usage_output_dirpath, gene_class, evalue_cutoff, drop_allele_info, ref_seq_dirpath, ig_expression_filepath, scaling_factor)
    return

def run_cdr3(cdr3_seqs_dirpath, cdr3_usage_output_filepath, ig_expression_filepath, scaling_factor):
    if cdr3_seqs_dirpath[-1] != '/':
        cdr3_seqs_dirpath += '/'

    tpoint_ids = []
    cdr3_seqs_input_filepaths = []
    for i in os.listdir(cdr3_seqs_dirpath):
        if i[0] == '.':
            continue
        tpoint_ids.append(int(i))
        cdr3_seqs_input_filepaths.append([int(i), cdr3_seqs_dirpath + i])
    cdr3_seqs_input_filepaths = [i[1] for i in sorted(cdr3_seqs_input_filepaths)]
    tpoint_ids = sorted(tpoint_ids)
    
    print 'getting dic of all cdr3 seqs'
    sys.stdout.flush()
    cdr3_seqs_dic = make_dic_of_all_cdrs(cdr3_seqs_input_filepaths)
    print 'getting dic for individual time-points'
    all_tpoints_cdr3_seqs_dics = get_cdr3_usage_foreach_tpoint(cdr3_seqs_input_filepaths, cdr3_seqs_dic)
    print "normalizing counts"
    total_mapped_reads = get_total_mapped_reads(ig_expression_filepath)
    all_tpoints_cdr3_seqs_dics = normalize_counts(all_tpoints_cdr3_seqs_dics, 'cdr3', total_mapped_reads, scaling_factor)
    print "writing output"
    write_output(all_tpoints_cdr3_seqs_dics, tpoint_ids, cdr3_usage_output_filepath, 'cdr3')
    return

if __name__ == '__main__':
    if len(sys.argv[1:]) > 9:
        run_array(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7], float(sys.argv[8]), sys.argv[9:])
    elif len(sys.argv[1:]) == 9:
        run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), sys.argv[6], sys.argv[7], sys.argv[8], float(sys.argv[9]))
    elif len(sys.argv[1:]) == 4:
        run_cdr3(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))
