import sys
import os
import itertools
import math
from scipy import stats

def get_row_sum_counts(row_sums_filepath):
    """This script gets the counts for each possible row sum. It first checks what format the file is in. If it is in the 'is_each_gene_sig_or_not' format, this script will calculate the counts for each possible row sum from this file, then return that as a list of the counts, where the 1st element of the list is the counts for how many rows summed to 0, and the 2nd element is how rows summed to 1, etc. This is encoded in 'row_sum_counts', which is retured. If the file is in the format of already having the row_sum_counts calculated, then it simply reads this data into memory, and returns it."""
    row_sum_counts = []
    filein = open(row_sums_filepath, "r")
    header = filein.readline().split('\t')
    #check if input file is a list of row sum counts, or is in the 
    #'is_each_gene_sig_or_not' format
    if len(header) > 3:
        #prime row_sum_counts list
        row_sum_counts = [0 for i in header]
        for i in filein:
            sigs = [int(j) for j in i[:-1].split('\t')[1:]]
            row_sum = sum(sigs)
            row_sum_counts[row_sum] += 1
    else:
        for i in filein:
            row_sum_counts.append(int(i[:-1].split('\t')[1]))
    filein.close()
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
    fileout.write('row_sum_count_values\tnull_expectation\tobserved_counts\tmultinom_likelihood_ratio_test_p_value\n')
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

def run(null_row_sums_dirpath, obs_row_sums_dirpath, output_dirpath):
    """This scripts purpose is to compare row sum counts from an observed dataset to that of a null dataset. First, what are 'row_sum_counts'? Well, they are what you get when you sum the rows of the 'is_each_gene_sig_or_not' files and then count how many times you get each possible number from these row sums. For eaxample, if there are 5 patients, then each row's sum will be somewhere between 0 and 5. Thus, row_sum_counts are the counts for how many times one gets a 0 or a 1 or a 2 etc. for these row sums. This module compares the two distributions by running a multinomial likelihood ratio test.
    'null_row_sums_dirpath' = this is the directory that contains the information for the row sums for the null data set. This file can be in the format of already being processed as row_sum_counts, or it can be in the form of 'is_each_gene_sig_or_not' and the row_sum_counts will be calculated from that.
    'obs_row_sums_dirpath' = same as the 'null_row_sums_dirpath', except for the observed dataset.
    'output_dirpath' = directory to write output."""
    if null_row_sums_dirpath[-1] != '/':
        null_row_sums_dirpath += '/'
    if obs_row_sums_dirpath[-1] != '/':
        obs_row_sums_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    for i in os.listdir(null_row_sums_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        print '::::::::::::::::::::'
        print i
        print '::::::::::::::::::::'
        output_sampletype_dirpath = output_dirpath + i + '/'
        if not os.path.exists(output_sampletype_dirpath):
            os.makedirs(output_sampletype_dirpath)
        for j in os.listdir(null_row_sums_dirpath + i):
            if j[0]== '.' or j == 'README':
                continue
            print '----------', j, '----------'
            output_filepath = output_sampletype_dirpath + j
            null_row_sums_filepath = null_row_sums_dirpath + i + '/' + j
            obs_row_sums_filepath = obs_row_sums_dirpath + i + '/' + j
            null_row_sum_counts = get_row_sum_counts(null_row_sums_filepath)
            obs_row_sum_counts = get_row_sum_counts(obs_row_sums_filepath)
            p_value = multinom_likelihood_ratio_test(obs_row_sum_counts, null_row_sum_counts)
            write_output(null_row_sum_counts, obs_row_sum_counts, p_value, output_filepath)
            print 'p value:', p_value
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3])
