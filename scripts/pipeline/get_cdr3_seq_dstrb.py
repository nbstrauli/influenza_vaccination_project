import sys

def get_cdr3_seq_dstrb(cdr3_seqs_filepath):
    """This script gets the counts and sequences of each unique CDR3 sequence contained in 'cdr3_seqs_filepath'. The input file is one cdr3 seq per line, of the form of starting with some sort of sequence id (just can't contain a tab) followed by a tab then the CDR3 sequence. If the CDR3 seq is 'Null', 'N/A', 'rejected', or '' (NoneType) then it will not consider that sequence. The output is a file where each line starts with a unique cdr3 seq then tab then the count of that seq."""
    filein = open(cdr3_seqs_filepath, "r")
    unique_cdr3s = {}
    for i in filein:
        cdr3_seq = i[:-1].split('\t')[1]
        if cdr3_seq == 'N/A' or cdr3_seq == 'Null' or cdr3_seq == 'rejected' or cdr3_seq == '':
            continue
        try:
            unique_cdr3s[cdr3_seq] += 1
        except KeyError:
            unique_cdr3s[cdr3_seq] = 1
    filein.close()
    sorted_unique_cdr3s_withseq = []
    sorted_unique_cdr3s_numberonly = []
    for i in unique_cdr3s:
        sorted_unique_cdr3s_withseq.append([unique_cdr3s[i], i])
        sorted_unique_cdr3s_numberonly.append(unique_cdr3s[i])
    sorted_unique_cdr3s_withseq = sorted(sorted_unique_cdr3s_withseq)#this is now a list in ascending order, where each element is of the form [number_of_copies_of_CDR3_found, sequence_of_the_CDR3]    
    unique_cdr3s_fileout = open(cdr3_seqs_filepath + '_unique_seqs', "w")
    for i in sorted_unique_cdr3s_withseq:
        unique_cdr3s_fileout.write(str(i[0]) + '\t' + str(i[1]) + '\n')
    unique_cdr3s_fileout.close()
    return

if __name__ == '__main__':
    get_cdr3_seq_dstrb(sys.argv[1])
