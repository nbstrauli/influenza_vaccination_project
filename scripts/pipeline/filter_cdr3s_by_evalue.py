import sys
import itertools

def filter(cdr3_seq_filepath, evalue_filepath, output_filepath, max_evalue):
    """This script filters the list of CDR3 sequences contained in 'cdr3_seq_filepath' based upon the list of corresponding e values contained in 'evalue_filepath'. If the evalue is less than or equal to the 'max_evalue' then this script will write the corresponding line from 'cdr3_seq_filepath' to 'output_filepath'."""
    filein1 = open(evalue_filepath, "r")
    filein2 = open(cdr3_seq_filepath, "r")
    fileout = open(output_filepath, "w")
    for i, j in itertools.izip(filein1, filein2):
        evalue = i[:-1].split('\t')[1]
        if evalue == 'Null':
            continue
        elif float(evalue) <= max_evalue:
            fileout.write(j)
    filein1.close()
    filein2.close()
    fileout.close()
    return

if __name__ == '__main__':
    filter(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))
