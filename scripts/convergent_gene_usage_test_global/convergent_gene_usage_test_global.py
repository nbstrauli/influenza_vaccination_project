import sys
import os
import random
import itertools
import math
from scipy import stats

def get_day0_pdf(freq_info):
    """This script simply converts the frequecies of the genes at day0 into a pdf. To do this it gets the first column from the frequency progession file (which is given by the filepath in 'freq_info'), and then normalizes each element of this resulting list of floats by the total of that list. This essentially is an empirical pdf. Returns a list of floats, which is the empirical probability of each gene."""
    """Added something here, where if 'freq_info' is actually a list, then it will treat that as a direct input of frequencies, and just turn that into a pdf. Added this so that the pdfs can be updated when one of the elements are removed from them below."""
    if type(freq_info) is str:
        pdf = []
        gene_names = []
        total = 0.0
        filein = open(freq_info, "r")
        filein.readline()
        for i in filein:
            line = i.split('\t')
            day0_freq = float(line[1])
            pdf.append(day0_freq)
            gene_name = line[0]
            gene_names.append(gene_name)
            total += day0_freq
        filein.close()
        if total == 0:
            return []
        else:
            pdf = [i/total for i in pdf]
            return pdf, gene_names
    else:
        pdf = freq_info
        total = sum(pdf)
        if total == 0:
            return []
        else:
            pdf = [i/total for i in pdf]
            return pdf

def convert_pdf2cdf(pdf):
    """This script takes the pdf and simply converts it to a cdf. Pretty straight forward."""
    total = 0
    cdf = []
    for i in pdf:
        cdf.append(total + i)
        total += i
    return cdf

def get_num_diff_expr(diff_expr_genes_filepath, FDR):
    filein = open(diff_expr_genes_filepath, "r")
    filein.readline()
    num_diff_expr = 0
    for i in filein:
        pval = float(i.split('\t')[1])
        if pval <= FDR:
            num_diff_expr += 1
    filein.close()
    return num_diff_expr

def inverse_transform_sampling(cdf):
    """This is a cool general way to sample from any discrete distribution (especially empirical ones). Generates a random # b/t 0 and 1, then finds the largest element in the cdf (x-axis) that has a value (y-axis) that is greater or equal to the random number. Returns the index of the element that it finds ('chosen_cdr3_index')."""
    random_float = random.random()
    for index, i in enumerate(cdf):
        if i >= random_float:
            chosen_index = index
            break
    return chosen_index

def sample_from_cdf(pdf, cdf, num_diff_expr, gene_names):
    chosen_genes = []
    for i in xrange(num_diff_expr):
        chosen_index = inverse_transform_sampling(cdf)
        del pdf[chosen_index]
        pdf = get_day0_pdf(pdf)
        cdf = convert_pdf2cdf(pdf)
        chosen_gene = gene_names[chosen_index]
        del gene_names[chosen_index]
        chosen_genes.append(chosen_gene)
    return chosen_genes

def get_SGS(chosen_genes_list, gene_name_lists):
    num_patients = len(chosen_genes_list)
    #first prime genename_dic
    genename_dic = {}
    for i in gene_name_lists:
        for j in i:
            try:
                genename_dic[j]
            except KeyError:
                genename_dic[j] = [0 for k in xrange(num_patients)]
    for patient_num, i in enumerate(chosen_genes_list):
        for j in i:
            genename_dic[j][patient_num] = 1
    return genename_dic

def update_row_sums(genename_dic, row_sums):
    """This script uses the information within 'genename_dic' (which is the sig_or_not info for 1 trial) and updates the counts for each possible row sum from that interation."""
    for i in genename_dic:
        row_sum = sum([int(j) for j in genename_dic[i]])
        row_sums[row_sum] += 1
    return row_sums

def get_obs_row_sum_counts(obs_sgs_filepath):
    """This script gets the counts for each possible row sum. It first checks what format the file is in. If it is in the 'is_each_gene_sig_or_not' format, this script will calculate the counts for each possible row sum from this file, then return that as a list of the counts, where the 1st element of the list is the counts for how many rows summed to 0, and the 2nd element is how rows summed to 1, etc. This is encoded in 'row_sum_counts', which is retured. If the file is in the format of already having the row_sum_counts calculated, then it simply reads this data into memory, and returns it."""
    row_sum_counts = []
    filein = open(obs_sgs_filepath, "r")
    #prime row_sum_counts list
    row_sum_counts = [0 for i in filein.readline().split('\t')[1:]]
    for i in filein:
        row_sum = int(i[:-1].split('\t')[-1])
        row_sum_counts[row_sum] += 1
    return row_sum_counts

def multinom_likelihood_ratio_test(obs_row_sum_counts, null_row_sum_counts):
    """This script runs the multinomial log likelihood ratio test, and returns the resulting p value. Also applies the Williams correction to the G statistic. This procedure is described on pg. 699 of the Biometry (Sokal) book in Box 17.1."""
    null_total_counts = float(sum(null_row_sum_counts))
    null_probs = [i/null_total_counts for i in null_row_sum_counts]
    obs_total_counts = float(sum(obs_row_sum_counts))
    log_ratio_sum = 0
    for obs_count, null_prob in itertools.izip(obs_row_sum_counts, null_probs):
        if obs_count == 0:
            continue
        log_ratio_sum += obs_count * math.log( null_prob/(obs_count/obs_total_counts) )
    G = -2 * log_ratio_sum

    #apply William's correction:
    a = len(obs_row_sum_counts)
    q = 1 + ( (a + 1)/(6 * obs_total_counts) )
    G_adj = G / q

    df = len(null_probs) - 1
    p_value = 1 - stats.chi2.cdf(G_adj, df)
    return p_value

def write_output(null_row_sum_counts, obs_row_sum_counts, p_value, output_filepath):
    """This script writes the null row sum counts and observed row sum counts as well as the multinomial likelihood ratio test p value to the output filepath."""
    row_sum_values = [str(i) for i in xrange(len(null_row_sum_counts))]
    fileout = open(output_filepath, "w")
    fileout.write('SGS_values\tobserved_counts\tnull_expectation\tmultinom_likelihood_ratio_test_p_value\n')
    count = 0
    start = 'yup'
    for i, j in itertools.izip(null_row_sum_counts, obs_row_sum_counts):
        if start:
            start = None
            fileout.write(str(count) + '\t' + str(i) + '\t' + str(j) + '\t' + str(p_value) + '\n')
        else:
            fileout.write(str(count) + '\t' + str(i) + '\t' + str(j) + '\n')
        count += 1
    fileout.close()
    return

def run(freq_prog_dirpath, diff_expr_genes_dirpath, obs_sgs_filepath, output_filepath, trials, FDR=0.05):
    if freq_prog_dirpath[-1] != '/':
        freq_prog_dirpath += '/'
    if diff_expr_genes_dirpath != '/':
        diff_expr_genes_dirpath += '/'

    pdfs = []
    cdfs = []
    gene_name_lists = []
    patient_ids = []
    num_diff_expr_list = []
    for i in sorted(os.listdir(freq_prog_dirpath)):
        if i[0] == '.' or i == 'README':
            continue
        freq_prog_filepath = freq_prog_dirpath + i
        pdf, gene_names = get_day0_pdf(freq_prog_filepath)
        pdfs.append(pdf)
        gene_name_lists.append(gene_names)
        cdf = convert_pdf2cdf(pdf)
        cdfs.append(cdf)
        patient_ids.append(i)
        diff_expr_genes_filepath = diff_expr_genes_dirpath + i
        num_diff_expr = get_num_diff_expr(diff_expr_genes_filepath, FDR)
        num_diff_expr_list.append(num_diff_expr)
    
    #Now it's time to sample from the cdfs (without replacement) to
    #create lists of genes sampled from day 0.
    #The index number in the below list corresponds to each possible row
    #sum, and the value corresponds to the count for that row sum
    row_sums = [0 for j in xrange(len(patient_ids) + 1)]
    #for each trial
    for i in xrange(trials):
        #Cycle through patients to sample from each cdf
        chosen_genes_list = []
        for j in xrange(len(cdfs)):
            #now sample genes from day0
             pdf = pdfs[j][:]
             cdf = cdfs[j][:]
             gene_names = gene_name_lists[j][:]
             #This is what samples the genes from the cdf
             chosen_genes = sample_from_cdf(pdf, cdf, num_diff_expr, gene_names)
             chosen_genes_list.append(chosen_genes)
        genename_dic = get_SGS(chosen_genes_list, gene_name_lists)
        #update row sum counts
        update_row_sums(genename_dic, row_sums)
    
    #now compare the observed SGS distribution to the simulated
    #null distribution
    obs_row_sum_counts = get_obs_row_sum_counts(obs_sgs_filepath)
    p_value = multinom_likelihood_ratio_test(obs_row_sum_counts, row_sums)
    write_output(obs_row_sum_counts, row_sums, p_value, output_filepath)
    return

if __name__ == '__main__':
    if len(sys.argv[1:]) == 6:
        run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]))
    else:
        run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
