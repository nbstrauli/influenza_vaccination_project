import sys
import os
import re
import itertools

def get_evalue(query, write_evalue, count):
    """This script will return the e value from the IgBLAST alignment. It will only return the evalue of the V gene (or more precisely the e value of the gene that is first listed in the list of e values, which I think is always the V gene). It writes it to a tab deliminated file where each line is the count of queries followed by the e value of the V gene."""
    evalue_regex = '^\S+\s+\S+\s+(\S+)\s*\n$'#This is kinda a weak regex. Should only be used with the if statment below
    for i in xrange(len(query)):
        if query[i][-14:-1] == '(Bits)  Value':
            evalue = re.search(evalue_regex, query[i+2]).group(1)
            write_evalue.write(str(count) + '\t' + evalue + '\n')
            return
    write_evalue.write(str(count) + '\tNull\n')
    return

def get_evalues_and_names_for_VDandJgenes(query, write_vgene_heavy_evalue, write_dgene_heavy_evalue, write_jgene_heavy_evalue, write_vgene_lambda_evalue, write_jgene_lambda_evalue, write_vgene_kappa_evalue, write_jgene_kappa_evalue, write_vgene_heavy_name, write_dgene_heavy_name, write_jgene_heavy_name, write_vgene_lambda_name, write_jgene_lambda_name, write_vgene_kappa_name, write_jgene_kappa_name, gene_names, count):
    evalue_regex = '^(\S+)\s+\S+\s+(\S+)\s*\n$'#This is kinda a weak regex. Should only be used with the if statment below
    for i in xrange(len(query)):
        #if this conditional is satisfied then that means that
        #the query has evalues in it. So, (staying within this
        #conditional statement) we will get the necessary info
        #for evalues and gene names and such and then exit the
        #function altogether. If the conditional is never
        #satisfied then it means that there are no evalues in
        #the query and we will write the necessary output for
        #that case.
        if query[i][-14:-1] == '(Bits)  Value':
            gene_names_and_evalues = {}#this is a dictionary where the index is the name of the gene and the def is its evalue
            for j in query[i+2:]:
                if j == '\n':
                    best_gene_names = {'IGHV':'N/A', 'IGHD':'N/A', 'IGHJ':'N/A', 'IGLV':'N/A', 'IGLJ':'N/A', 'IGKV':'N/A', 'IGKJ':'N/A'}
                    best_evalues = {'IGHV':'N/A', 'IGHD':'N/A', 'IGHJ':'N/A', 'IGLV':'N/A', 'IGLJ':'N/A', 'IGKV':'N/A', 'IGKJ':'N/A'}
                    list_of_gene_classes = ['IGHV', 'IGHD', 'IGHJ', 'IGLV', 'IGLJ', 'IGKV', 'IGKJ']
                    #Here we get the gene name and evalue for the best aligning
                    #germline genes for each gene class. Note that this
                    #approach only works if igblast is set to report only 1
                    #alignment per gene class (ex: -num_alignments_V 1)
                    for k in gene_names_and_evalues:
                        for l in list_of_gene_classes:
                            if best_gene_names[l] == 'N/A':
                                try:
                                    gene_names[l][k]
                                    best_gene_names[l] = k
                                    best_evalues[l] = gene_names_and_evalues[k]
                                    continue
                                except KeyError:
                                    pass
                    gene_name_output_files = [write_vgene_heavy_name, write_dgene_heavy_name, write_jgene_heavy_name, write_vgene_lambda_name, write_jgene_lambda_name, write_vgene_kappa_name, write_jgene_kappa_name]
                    evalue_output_files = [write_vgene_heavy_evalue, write_dgene_heavy_evalue, write_jgene_heavy_evalue, write_vgene_lambda_evalue, write_jgene_lambda_evalue, write_vgene_kappa_evalue, write_jgene_kappa_evalue]
                    for k, l, m in itertools.izip(gene_name_output_files, evalue_output_files, list_of_gene_classes):
                        k.write(str(count) + '\t' + best_gene_names[m] + '\n')
                        l.write(str(count) + '\t' + best_evalues[m] + '\n')
                    return
                else:
                    match = re.search(evalue_regex, j)
                    if match.group(1)[:4] == 'lcl|': #sometimes the genenames have 'lcl|' in front of them and sometimes not
                        gene_names_and_evalues[match.group(1)[4:]] = match.group(2)
                    else:
                        gene_names_and_evalues[match.group(1)] = match.group(2)
    write_vgene_heavy_evalue.write(str(count) + '\tN/A\n')
    write_dgene_heavy_evalue.write(str(count) + '\tN/A\n')
    write_jgene_heavy_evalue.write(str(count) + '\tN/A\n')
    write_vgene_lambda_evalue.write(str(count) + '\tN/A\n')
    write_jgene_lambda_evalue.write(str(count) + '\tN/A\n')
    write_vgene_kappa_evalue.write(str(count) + '\tN/A\n')
    write_jgene_kappa_evalue.write(str(count) + '\tN/A\n')
    write_vgene_heavy_name.write(str(count) + '\tN/A\n')
    write_dgene_heavy_name.write(str(count) + '\tN/A\n')
    write_jgene_heavy_name.write(str(count) + '\tN/A\n')
    write_vgene_lambda_name.write(str(count) + '\tN/A\n')
    write_jgene_lambda_name.write(str(count) + '\tN/A\n')
    write_vgene_kappa_name.write(str(count) + '\tN/A\n')
    write_jgene_kappa_name.write(str(count) + '\tN/A\n')
    return

def get_cdr3_seq(query, write_cdr3_seq, count):
    DNAseq_regex = '.*Query.*?([0-9]+)\s+([ATGCN-]+)\s+([0-9]+)'
    #The regular expression below only works when '-domain_sysem' parameter
    #in IgBLAST is set to 'imgt'. Change this regex if you want to use a 
    #different domain system.
    FWR3_regex = 'FR3-IMGT\t[0-9]+\t([0-9]+)\t'
    cdr3End_regex = '(([TACG]{3})*)T[TG][TCG]GG[TACG][TACG][TACG][TACG]GG[TACG]'
    DNAseq = ''
    DNAseq_line = 1
    for i in query:
        if i[0:6] == "Query=":
            id = i[7:-1]
        match = re.search(FWR3_regex, i)
        if match:
            cdr3_start = int(match.group(1))
        match = re.search(DNAseq_regex, i)
        if match:
            DNAseq = DNAseq + match.group(2)
            if DNAseq_line == 1:
                DNAseq_start = int(match.group(1)) - 1
            DNAseq_line += 1
    try:
        cdr3ToEnd = DNAseq[cdr3_start - DNAseq_start:]
        match = re.search(cdr3End_regex, cdr3ToEnd)
        cdr3_seq = match.group(1)
        write_cdr3_seq.write(id + "\t" + cdr3_seq + "\n")
    except:
        cdr3_seq = "Null"
        cdr3_start = None
        write_cdr3_seq.write(id + "\t" + cdr3_seq + "\n")
    return id, cdr3_start, cdr3_seq

def parse_dat(igblast_output, ref_seqs_dirpath):
    """This script"""
    if ref_seqs_dirpath[-1] != '/':
        ref_seqs_dirpath += '/'
    filein = open(igblast_output, "r")
    write_cdr3_seq = open(igblast_output + "_cdr3_seq", "w")
    write_vgene_heavy_evalue = open(igblast_output + '_vgene_heavy_evalue', "w")
    write_dgene_heavy_evalue = open(igblast_output + '_dgene_heavy_evalue', "w")
    write_jgene_heavy_evalue = open(igblast_output + '_jgene_heavy_evalue', "w")
    write_vgene_lambda_evalue = open(igblast_output + '_vgene_lambda_evalue', "w")
    write_jgene_lambda_evalue = open(igblast_output + '_jgene_lambda_evalue', "w")
    write_vgene_kappa_evalue = open(igblast_output + '_vgene_kappa_evalue', "w")
    write_jgene_kappa_evalue = open(igblast_output + '_jgene_kappa_evalue', "w")
    write_vgene_heavy_name = open(igblast_output + '_vgene_heavy_name', "w")
    write_dgene_heavy_name = open(igblast_output + '_dgene_heavy_name', "w")
    write_jgene_heavy_name = open(igblast_output + '_jgene_heavy_name', "w")
    write_vgene_lambda_name = open(igblast_output + '_vgene_lambda_name', "w")
    write_jgene_lambda_name = open(igblast_output + '_jgene_lambda_name', "w")
    write_vgene_kappa_name = open(igblast_output + '_vgene_kappa_name', "w")
    write_jgene_kappa_name = open(igblast_output + '_jgene_kappa_name', "w")
    write_evalue = open(igblast_output + '_evalues', "w")
    #retrieve the names of all the reference genes
    #for each gene class
    gene_names = {}
    for i in os.listdir(ref_seqs_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        gene_names[i[:-6]] = {}
        filein_ref_seq = open(ref_seqs_dirpath + i, "r")
        for j in filein_ref_seq:
            if j[0] == '>':
                gene_names[i[:-6]][j[1:-1]] = None
        filein_ref_seq.close()
    query = []
    start = 0
    count = 0
    for i in filein:
        if i[0:6] == "Query=":
            start = 1
            count += 1
        if start:
            query.append(i)
        if i[0:16] == "Effective search":
            id, cdr3_start, cdr3_seq = get_cdr3_seq(query, write_cdr3_seq, count)
            get_evalues_and_names_for_VDandJgenes(query, write_vgene_heavy_evalue, write_dgene_heavy_evalue, write_jgene_heavy_evalue, write_vgene_lambda_evalue, write_jgene_lambda_evalue, write_vgene_kappa_evalue, write_jgene_kappa_evalue, write_vgene_heavy_name, write_dgene_heavy_name, write_jgene_heavy_name, write_vgene_lambda_name, write_jgene_lambda_name, write_vgene_kappa_name, write_jgene_kappa_name, gene_names, count)
            get_evalue(query, write_evalue, count)
            start = 0
            query = []
    write_cdr3_seq.close()
    write_vgene_heavy_evalue.close()
    write_dgene_heavy_evalue.close()
    write_jgene_heavy_evalue.close()
    write_vgene_lambda_evalue.close()
    write_jgene_lambda_evalue.close()
    write_vgene_kappa_evalue.close()
    write_jgene_kappa_evalue.close()
    write_vgene_heavy_name.close()
    write_dgene_heavy_name.close()
    write_jgene_heavy_name.close()
    write_vgene_lambda_name.close()
    write_jgene_lambda_name.close()
    write_vgene_kappa_name.close()
    write_jgene_kappa_name.close()
    write_evalue.close()
    filein.close()
    return

if __name__ == '__main__':
    parse_dat(sys.argv[1], sys.argv[2])
