import os
import sys
import subprocess

def get_genename_dic(fpca_test_output_dipath):
    """This script gets the genes names for all the genes found in each patient. It does this by making a dictionary, where each index is a gene name and then it makes the definition for that index a list of '0''s of length equal to the number of patients. '0' denotes "not significant" in this case, so this script is essentially priming the genenames_dic, which will then be updated with which genes are significant later."""
    genename_dic = {}
    patients = [i for i in sorted(os.listdir(fpca_test_output_dipath)) if i[0] != '.']
    for i in os.listdir(fpca_test_output_dipath):
        if i[0] == '.' or i == 'README':
            continue
        FPCAtest_output_filepath = fpca_test_output_dipath + i
        filein = open(FPCAtest_output_filepath, "r")
        filein.readline()
        for j in filein:
            gene_name = j.split('\t')[0]
            try:
                genename_dic[gene_name]
            except KeyError:
                genename_dic[gene_name] = [0] * len(patients)
        filein.close()
    return genename_dic, patients

def get_SGS(fpca_test_output_dipath, genename_dic, patients, FDR):
    patient_num = 0
    for i in patients:
        FPCAtest_output_filepath = fpca_test_output_dipath + i
        filein = open(FPCAtest_output_filepath, "r")
        filein.readline()
        for j in filein:
            line = j.split('\t')
            gene_name = line[0]
            pval = float(line[1])
            if pval <= FDR:
                genename_dic[gene_name][patient_num] = 1
        filein.close()
        patient_num += 1
    return genename_dic

def sort_and_write_output(genename_dic, SGS_output_filepath, patients):
    #sort genes according to their SGS value
    genename_list = []
    for i in genename_dic:
        SGS = sum(genename_dic[i])
        entry = [SGS, i]
        entry.extend(genename_dic[i])
        genename_list.append(entry)
    genename_list = sorted(genename_list, reverse=True)
    fileout = open(SGS_output_filepath, "w")
    fileout.write('gene_name\t' + '\t'.join(patients) + '\tSGS\n')
    for i in genename_list:
        sigs = [str(j) for j in i[2:]]
        fileout.write(i[1] + '\t' + '\t'.join(sigs) + '\t' + str(i[0]) + '\n')
    fileout.close()
    return

def run(gene_usage_dirpath, output_dirpath, FDR=0.05):
    """'gene_usage_dirpath' = The path to the directory that contains the gene expression trajectories for each gene. More specifically, each file in this directory should be a gene usage matrix, where each row is a gene, each column is a time-point, and each element is the expression level for that gene at that time-point. The gene name is the first entry in each row, and each column has a header giving the unique, numeric, time-point ID. The name of each of these files should be a unique identifier for each individual patient in the study. Essentially, each of these files must be in the same format as the 'gene_usage' output from the 'pipeline.bash' script."""
    if gene_usage_dirpath[-1] != '/':
        gene_usage_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    #First identiy the differentially expressed genes
    fpca_test_output_dipath = output_dirpath + 'FPCA_test_results/'
    if not os.path.exists(fpca_test_output_dipath):
        os.makedirs(fpca_test_output_dipath)
    for i in os.listdir(gene_usage_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        gene_usage_filepath = gene_usage_dirpath + i
        fpca_test_output_filepath = fpca_test_output_dipath + i
        subprocess.call(['Rscript', 'call_FPCA_test.R', gene_usage_filepath, fpca_test_output_filepath])
    #Now get the number of patients each gene is significant in (SGS)
    genename_dic, patients = get_genename_dic(fpca_test_output_dipath)
    genename_dic = get_SGS(fpca_test_output_dipath, genename_dic, patients, FDR)
    SGS_output_filepath = output_dirpath + 'SGS_table'
    sort_and_write_output(genename_dic, SGS_output_filepath, patients)
    return

if __name__ == '__main__':
    if sys.argv[1:] == 3:
        run(sys.argv[1], sys.argv[2], float(sys.argv[3]))
    else:
        run(sys.argv[1], sys.argv[2])
