import sys
import os
home_dirpath = os.getenv("HOME") + '/'
sys.path.insert(0, home_dirpath + 'Scripts/Py_Scripts/Snyderome/')
import simulate_removing_reads
import subprocess
import re
import get_each_gene_significant_or_not_foreach_patient
import test_for_similar_significant_genes_across_patients
import G_test_for_independence_for_2X2_table

def get_day0_pdf(freq_prog_filepath):
    """This script simply converts the frequecies of the genes at day0 into a pdf. To do this it gets the first column from 'freq_prog_filepath', and then normalizes each element of this resulting list of floats by the total of that list. This essentially is an empirical pdf. Returns a list of floats, which is the empirical probability of each gene."""
    """Added something here, where if 'freq_prog_filepath' is actually a list, then it will treat that as a direct input of frequencies, and just turn that into a pdf. Added this so that the pdfs can be updated when one of the elements are removed from them below."""
    if type(freq_prog_filepath) is str:
        pdf = []
        total = 0.0
        filein = open(freq_prog_filepath, "r")
        filein.readline()
        for i in filein:
            day0_freq = float(i.split('\t')[0])
            pdf.append(day0_freq)
            total += day0_freq
        filein.close()
    else:
        pdf = freq_prog_filepath
        total = sum(pdf)
    if total == 0:
        return []
    pdf = [i/total for i in pdf]
    return pdf

def get_gene_names(gene_name_order_filepath):
    """This script simply gets the order of gene names that correspond to the pdf that is created above."""
    gene_names = []
    filein = open(gene_name_order_filepath, 'r')
    filein.readline()
    for i in filein:
        gene_names.append(i[:-1])
    filein.close()
    return gene_names

def get_num_diff_expr(diff_expr_genes_subdirpath):
    """This script cycles through the directory that contains the expression function plots for each of the genes that were deemed differentially expressed by the FPCAtest ('diff_expr_genes_subdirpath'). It gets the number of files in this directory that end in '.pdf', as this will correspond to the total number of these genes that are differentially expressed."""
    num_diff_expr = 0
    for i in os.listdir(diff_expr_genes_subdirpath):
        if i[-4:] == '.pdf':
            num_diff_expr += 1
    if num_diff_expr == 0:
        return num_diff_expr, []
    FPCAtest_output_table_filepath = diff_expr_genes_subdirpath + 'FPCAtest_output_as_table'
    filein = open(FPCAtest_output_table_filepath, "r")
    filein.readline()
    FPCAtest_output_genenames = []
    for i in filein:
        FPCAtest_output_genenames.append(i.split('\t')[0])
    filein.close()
    return num_diff_expr, FPCAtest_output_genenames

def make_temp_chosen_genes_filepaths(temp_sampled_genes_dirpath, pdf, cdf, num_diff_expr, gene_names, FPCAtest_output_genenames):
    """This script is super important for the module! It first makes a mock file that contains all of the gene names that where found in the 'FPCAtest_output_as_table' file from the differentially expressed genes directory. It then randomly samples genes (without replacement) from the cdf (which was made from the gene frequencies at day0). When a gene is selected, a mock '.pdf' file is created with its gene name as the name of the file."""
    fileout = open(temp_sampled_genes_dirpath + 'FPCAtest_output_as_table', 'w')
    fileout.write('gene_names\tetc\n')
    for i in FPCAtest_output_genenames:
        fileout.write(i + '\tblah\n')
    fileout.close()
    for l in xrange(num_diff_expr):
        chosen_index = simulate_removing_reads.inverse_transform_sampling(cdf)
        del pdf[chosen_index]
        pdf = get_day0_pdf(pdf)
        cdf = simulate_removing_reads.convert_pdf2cdf(pdf)
        chosen_gene = gene_names[chosen_index]
        del gene_names[chosen_index]
        chosen_gene = re.sub('/', '_', chosen_gene)
        temp_chosen_gene_filepath = temp_sampled_genes_dirpath + chosen_gene + '_plot.pdf'
        f = open(temp_chosen_gene_filepath, "w")
        f.close()
    return

def update_row_sums(genename_dic, row_sums):
    """This script uses the information within 'genename_dic' (which is the sig_or_not info for 1 trial) and updates the counts for each possible row sum from that interation."""
    for i in genename_dic:
        row_sum = sum([int(j) for j in genename_dic[i]])
        row_sums[row_sum] += 1
    return row_sums

def write_row_sum_counts(row_sums, sim_row_sums_filepath):
    """This script writes the row sum counts information to file."""
    fileout_row_sums = open(sim_row_sums_filepath, "w")
    fileout_row_sums.write('num_patients_that_shared_gene' + '\t' + 'num_of_genes_for_each_category' + '\n')
    for i in xrange(len(row_sums)):
        fileout_row_sums.write(str(i) + '\t' + str(row_sums[i]) + '\n')
    fileout_row_sums.close()
    return

def write_output(pvalues_dic_G, pvalues_dic_fisher, output_sampletype_dirpath, patient_ids, gene_class):
    """This script writes the output. It does this for each patient. And it also writes a seperate file for each of the tests run (G test for independence and Fisher's exact test). The output simply consists of a long list of p values that were obtained from each of the trials."""
    for i in patient_ids:
        output_filepath = output_sampletype_dirpath + i + '/' + gene_class + '_Gtest_pvalues'
        fileout = open(output_filepath, "w")
        fileout.write('p_values\n')
        for j in sorted(pvalues_dic_G[i]):
            fileout.write(str(j) + '\n')
        fileout.close()
        output_filepath = output_sampletype_dirpath + i + '/' + gene_class + '_Fishers_exact_pvalues'
        fileout = open(output_filepath, "w")
        fileout.write('p_values\n')
        for j in sorted(pvalues_dic_fisher[i]):
            fileout.write(str(j) + '\n')
        fileout.close()
    return

def run(freq_prog_dirpath, gene_name_order_dirpath, diff_expr_genes_dirpath, output_dirpath, sim_row_sums_dirpath, trials, gene_classes):
    """This script's purpose is to sample a group of genes from each patient, based upon the frequency of the genes at day0, for each of the patients. The number of genes selected for each patient is equal to the number of genes that were significantly differentially expressed from the FPCAtest. The idea is to simulate a set of genes that are analogous to the set of differentially expressed genes, except this simulated set are selected purely based upon their freq before the vaccine response (day0). Once this null set is selected the 'each_gene_sig_or_not' file is created based on this set, and then the tests for independence are run on this file (both G test, and Fisher's exact test). This 'run' script is quite long. The reason for this is that, within a sampletype, the script has to first cycle through each of the patients to get the gene frequencies at day0, and then use this information to make the simulated sets of genes. This creates kind of a complicated for loop structure. Not so elegent... This list of genes, and then downstream tests are run 'trials' number of times.
    'freq_prog_dirpath' = This is the directory that contains the freq progressions for each of the genes.
    'gene_name_order_dirpath' = This directory contains the files that are the gene name order that correspons to the freq prog files.
    'diff_expr_genes_dirpath' = This is the directory that contains the plots of gene expression functions for the differentially expressed genes
    'output_dirpath' = This is where the output goes, which is file that is a long list of each of the pvalues from the tests for independence
    'trials' = Number of trials to sample the genes and run the tests etc.
    'gene_classes' = a list of names, where each name is a gene class to run this module on.
    'sim_row_sums_dirpath' = This is an output directory where the simulated 'sig_or_not' information will be written. This variable (and output) was added later. Essentially, the information in this file will be a simulated null distribution, and we want to compare this null to the observed. So, we are saving the info to disk, and then will compare later."""
    if freq_prog_dirpath[-1] != '/':
        freq_prog_dirpath += '/'
    if gene_name_order_dirpath[-1] != '/':
        gene_name_order_dirpath += '/'
    if diff_expr_genes_dirpath != '/':
        diff_expr_genes_dirpath += '/'
    if output_dirpath != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)

    if sim_row_sums_dirpath[-1] != '/':
        sim_row_sums_dirpath += '/'
    if not os.path.exists(sim_row_sums_dirpath):
        os.makedirs(sim_row_sums_dirpath)

    temp_samplegenes_dirpath = os.path.dirname(output_dirpath[:-1]) + '/temp/'
    if os.path.exists(temp_samplegenes_dirpath):
        subprocess.call(['rm', '-r', temp_samplegenes_dirpath])
    os.makedirs(temp_samplegenes_dirpath)
    temp_is_sig_or_not_dirpath = os.path.dirname(output_dirpath[:-1]) + '/sig_or_not/'
    if os.path.exists(temp_is_sig_or_not_dirpath):
        subprocess.call(['rm', '-r', temp_is_sig_or_not_dirpath])
    os.makedirs(temp_is_sig_or_not_dirpath)
    #for each sample type
    for i in os.listdir(freq_prog_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        print '--------------------', i, '--------------------'
        output_sampletype_dirpath = output_dirpath + i + '/'
        if not os.path.exists(output_sampletype_dirpath):
            os.makedirs(output_sampletype_dirpath)

        sim_row_sums_sampletype_dirpath = sim_row_sums_dirpath + i + '/'
        if not os.path.exists(sim_row_sums_sampletype_dirpath):
            os.makedirs(sim_row_sums_sampletype_dirpath)

        #for each gene class
        for gene_class in gene_classes:
            print '============='
            print gene_class
            print '============='
            pdfs = []
            gene_name_lists = []
            num_diff_expr_list = []
            FPCAtest_output_genename_lists = []
            cdfs = []
            patient_ids = []
            no_genes = None
            #cycle through patients to get each patients cdf and such
            for j in os.listdir(freq_prog_dirpath + i):
                if j[0] == '.' or j == 'README':
                    continue
                output_patient_dirpath = output_sampletype_dirpath + j + '/'
                if not os.path.exists(output_patient_dirpath):
                    os.makedirs(output_patient_dirpath)
                freq_prog_filepath = freq_prog_dirpath + i + '/' + j + '/' + gene_class
                pdf = get_day0_pdf(freq_prog_filepath)
                if len(pdf) == 0:
                    no_genes = 'yup'
                    break
                pdfs.append(pdf)
                gene_name_order_filepath = gene_name_order_dirpath + i + '/' + j + '/' + gene_class
                gene_names = get_gene_names(gene_name_order_filepath)
                gene_name_lists.append(gene_names)
                diff_expr_genes_subdirpath = diff_expr_genes_dirpath + i + '/' + j + '/' + gene_class + '/'
                num_diff_expr, FPCAtest_output_genenames = get_num_diff_expr(diff_expr_genes_subdirpath)
                num_diff_expr_list.append(num_diff_expr)
                FPCAtest_output_genename_lists.append(FPCAtest_output_genenames)
                cdf = simulate_removing_reads.convert_pdf2cdf(pdf)
                cdfs.append(cdf)
                patient_ids.append(j)
            if no_genes:
                print 'No genes found for', gene_class, 'in', j, 'so whole gene class aborted for', i
                continue
            #Now it's time to sample from the cdfs (without replacement) to
            #create lists of genes sampled from day 0
            #First prime the p values dictionary for each patient
            pvalues_dic_G = {}
            pvalues_dic_fisher = {}
            for j in patient_ids:
                pvalues_dic_G[j] = []
                pvalues_dic_fisher[j] = []

            #the index number in this list corresponds to each possible row
            #sum, and the value corresponds to the count for that row sum
            row_sums = [0 for j in xrange(len(patient_ids) + 1)]

            #For each trial
            for j in xrange(trials):
                #make tree of temp input files for each patient
                temp_samplegenes_sampletype_dirpath = temp_samplegenes_dirpath + i + '/'
                os.makedirs(temp_samplegenes_sampletype_dirpath)
                #Cycle through patients to sample from each cdf
                for k in xrange(len(cdfs)):
                    temp_samplegenes_patient_dirpath = temp_samplegenes_sampletype_dirpath + patient_ids[k] + '/'
                    os.makedirs(temp_samplegenes_patient_dirpath)
                    temp_samplegenes_final_dirpath = temp_samplegenes_patient_dirpath + gene_class + '/'
                    os.makedirs(temp_samplegenes_final_dirpath)
                    #now sample genes from day0
                    pdf = pdfs[k][:]
                    cdf = cdfs[k][:]
                    gene_names = gene_name_lists[k][:]
                    #This is what samples the genes from the cdf
                    make_temp_chosen_genes_filepaths(temp_samplegenes_final_dirpath, pdf, cdf, num_diff_expr_list[k], gene_names, FPCAtest_output_genename_lists[k])
                    ds = os.listdir(temp_samplegenes_final_dirpath)
                #now asses independence b/t patients (i.e. get p values)
                
                genename_dic = get_each_gene_significant_or_not_foreach_patient.run(temp_samplegenes_dirpath, temp_is_sig_or_not_dirpath, [gene_class], 'yes')
                #update row sum counts
                row_sums = update_row_sums(genename_dic, row_sums)
                
                sig_or_not_filepath = temp_is_sig_or_not_dirpath + i + '/' + gene_class
                for k in xrange(len(patient_ids)): #get p values for each patient
                    pvalue_G, pvalue_fisher = test_for_similar_significant_genes_across_patients.run(sig_or_not_filepath, 0, 0, k)
                    pvalues_dic_G[patient_ids[k]].append(pvalue_G)
                    pvalues_dic_fisher[patient_ids[k]].append(pvalue_fisher)
                #now remove temp files
                subprocess.call(['rm', '-r', temp_samplegenes_dirpath + i])
                subprocess.call(['rm', '-r', temp_is_sig_or_not_dirpath])

            #now write the row sum counts to file for this gene_class
            sim_row_sums_filepath = sim_row_sums_sampletype_dirpath + gene_class
            write_row_sum_counts(row_sums, sim_row_sums_filepath)

            #These are the P values for the Gtest for independence results of the
            #genes that are actually differentially expressed
            write_output(pvalues_dic_G, pvalues_dic_fisher, output_sampletype_dirpath, patient_ids, gene_class)
    subprocess.call(['rm', '-r', temp_samplegenes_dirpath])
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7:])
