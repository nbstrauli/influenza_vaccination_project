import sys
import os
import subprocess
import re

def run(input_fasta_dirpath):
    """This script is very simple, and not very important. It simply reformats the fasta files that were downloaded from the IMGT website, so that there are no white space lines, and that each sequence for each gene is only one line. Note that this script loads all the sequences into memory, so it would not be advised to use this on very large fasta files. Also, this script will over-write your input fasta files, so beware of that too."""
    if input_fasta_dirpath[-1] != '/':
        input_fasta_dirpath += '/'
    for i in os.listdir(input_fasta_dirpath):
        if i[-6:] == '.fasta':
            seq_data = {}
            filein = open(input_fasta_dirpath + i, "r")
            for j in filein:
                if j[0] == '>':
                    gene_name = j[1:-1].split('|')[1]
                    seq_data[gene_name] = ''
                elif j[0] != "\s":
                    #imgt sequences can sometimes have annoying '.' (gaps)
                    #in them so we remove these
                    seq = re.sub('\.', '', j[:-1])
                    seq_data[gene_name] += seq
            filein.close()
            subprocess.call(['rm', input_fasta_dirpath + i])
            fileout = open(input_fasta_dirpath + i, "w")
            for i in sorted(seq_data):
                fileout.write('>' + i + '\n' + seq_data[i] + '\n')
            fileout.close()
    return

if __name__ == '__main__':
    run(sys.argv[1])
