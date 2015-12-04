import sys
import os
import multiple_urn_multiple_pulls_problem

def get_obs_row_sum_genename_dic(SGS_table_filepath):
    """This script gets the observed number of patients that each gene was found to be significantly differentially expressed in. Does this by looping through each line (each gene) in the provided 'sig_or_not_filepath'. It returns a dictionary, where each genename is an entry, and has a definition of the number of patients that that gene was differentially expressed in (row sum). Also returns the number of patients in the study, for downstream stuff."""
    obs_row_sum_dic = {}
    filein = open(SGS_table_filepath, "r")
    #first get the number of patients
    header = filein.readline()
    num_of_patients = 0
    num_sig_genes_foreach_patient = []
    for i in header.split('\t')[1:-1]:
        num_of_patients += 1
        num_sig_genes_foreach_patient.append(0)
    #now get SGS value for each gene
    for i in filein:
        line = i[:-1].split('\t')
        genename = line[0]
        SGS = int(line[-1])
        obs_row_sum_dic[genename] = SGS
        patient_num = 0
        for j in line[1:-1]:
            num_sig_genes_foreach_patient[patient_num] += int(j)
            patient_num += 1
    filein.close()
    return obs_row_sum_dic, num_of_patients, num_sig_genes_foreach_patient

def prep_dic_of_gene_probs(obs_row_sum_dic, num_of_patients):
    """This script simply the dictionary of gene probabilities for each patient. This dictionary has the form of each entry is a gene name, and it is defined by a list of probabilities, where each prob in the list is the prob of selecting that gene in the corresponding patient (patient 1's prob corresponds to the 1st element in the list, patient 2 the 2nd, etc). This script, however, simply preps this dic by making each prob for each patient/gene equal to 0."""
    dic_of_gene_probs = {}
    definition = [0 for i in xrange(num_of_patients)]
    for i in obs_row_sum_dic:
        dic_of_gene_probs[i] = definition[:]
    return dic_of_gene_probs

def get_gene_probs(freq_prog_dirpath, dic_of_gene_probs):
    """This script actually fills in the probabilities for each gene in each patient into 'dic_of_gene_probs'. It does this by first cycling through patients, and then getting each genes frequency at day 0 for that patient. Once this is done it normalizes each of these genes frequencies by the total of the frequencies for that patient. This makes the probability of selecting each of the genes at day0. Fills in these values for each of the patients and each of the gene into 'dic_of_gene_probs'."""
    patient_num = 0
    total = 0
    #gets day0 freqs as they are in the freq prog file
    for i in os.listdir(freq_prog_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        freq_prog_filepath = freq_prog_dirpath + i
        filein = open(freq_prog_filepath, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            gene_name = line[0]
            day0_freq = float(line[1])
            dic_of_gene_probs[gene_name][patient_num] = day0_freq
            total += day0_freq
        #now normalize each of the day0 freqs by the total
        #if total != 0
        if total == 0:
            continue
        for j in sorted(dic_of_gene_probs):
            dic_of_gene_probs[j][patient_num] /= total
        patient_num += 1
    return dic_of_gene_probs

def convert_pdf2cdf(pdf):
    """This script takes the pdf and simply converts it to a cdf. Pretty straight forward."""
    total = 0
    cdf = []
    for i in pdf:
        cdf.append(total + i)
        total += i
    return cdf

def get_pvalues(dic_of_gene_probs, obs_row_sum_dic, num_sig_genes_foreach_patient):
    """This script uses the list of gene probs for each gene to get the pdf for the probability of selecting a given gene in X number of patients, where 0<=X<=5. Essentially, this pdf is our null distribution. This script simply cycles through each of the genes and get their list of probabilities in each patient, and it then passes this list to the script 'multiple_urn_problem.py' to actually get the pdf."""
    pvalues = []
    pdf_modes = {}
    for i in dic_of_gene_probs:
        if obs_row_sum_dic[i] == 0:
            pvalue = 1
            pdf = multiple_urn_multiple_pulls_problem.run(dic_of_gene_probs[i], num_sig_genes_foreach_patient)
            pdf_mode = pdf.index(max(pdf))
        else:
            pdf = multiple_urn_multiple_pulls_problem.run(dic_of_gene_probs[i], num_sig_genes_foreach_patient)
            pdf_mode = pdf.index(max(pdf))
            cdf = convert_pdf2cdf(pdf)
            pvalue = 1 - cdf[obs_row_sum_dic[i]-1]
        pvalues.append([pvalue, i])
        pdf_modes[i] = pdf_mode
    return pvalues, pdf_modes

def write_output(pvalues, obs_row_sum_dic, dic_of_gene_probs, pdf_modes, output_filepath):
    """This script writes the output. Each gene is a row in the file, and various info for the gene is reported (tab delimited), which is described ing the header."""
    fileout = open(output_filepath, "w")
    fileout.write('gene_name\tp_value\tobs_SGS\tmode_of_null_dstrb_SGS\tsum_of_day0_freqs\n')
    for i in sorted(pvalues):
        obs_row_sum_dic[i[1]]
        pdf_modes[i[1]]
        dic_of_gene_probs[i[1]]
        fileout.write(i[1] + '\t' + str(i[0]) + '\t' + str(obs_row_sum_dic[i[1]]) + '\t' + str(pdf_modes[i[1]]) + '\t' + str(sum(dic_of_gene_probs[i[1]])) + '\n')
    fileout.close()
    return

def run(SGS_table_filepath, freq_prog_dirpath, output_filepath):
    """This module's purpose is to run a gene by gene test to determine which genes are significantly differentially expressed in more patients than would be expected by chance. It does this by getting a probability of selecting a given gene for each patient, based on that gene's expression level at day0. Based on this it creates a null distribution for the expectation of the number of patients that the given gene will be selected in. This null distribution is then compared to the observed number of patients that the given gene was found to be differentially expressed, and a p value is then calculated.
    'sig_or_not_dirpath' = Path to the directory that has, for each gene, whether or not it was differentially expressed for each patient.
    'freq_prog_dirpath' = path to the directory that has the frequency progressions for each gene. The values in this file will be normalized by column (i.e. each value is normalized by the sum of its column), so there is not need to worry about how the file was normed.
    'gene_name_order_dirpath' = These are the files that have the gene names for each of the rows in the corresponding file in 'freq_prog_dirpath'.
    'output_dirpath' = Output is written here.
    'gene_classes' = List of the names of the gene classes so they can be cycled through.
    UPDATE: This module was changed a bit so that it runs the 'multiple_urn_multiple_pulls_problem.py' module, which should be a better null model for the data."""
    if freq_prog_dirpath[-1] != '/':
        freq_prog_dirpath += '/'
    obs_row_sum_dic, num_of_patients, num_sig_genes_foreach_patient = get_obs_row_sum_genename_dic(SGS_table_filepath)
    dic_of_gene_probs = prep_dic_of_gene_probs(obs_row_sum_dic, num_of_patients)
    dic_of_gene_probs = get_gene_probs(freq_prog_dirpath, dic_of_gene_probs)
    pvalues, pdf_modes = get_pvalues(dic_of_gene_probs, obs_row_sum_dic, num_sig_genes_foreach_patient)
    write_output(pvalues, obs_row_sum_dic, dic_of_gene_probs, pdf_modes, output_filepath)
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3])
