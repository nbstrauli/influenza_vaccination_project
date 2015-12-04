import os
import sys

def fast_Q2A(fastq_filepath):
    """This script will convert a fastq file to a fasta file that has the same name as the original fastq file except the suffix to the file name will change from fastq to fasta. The input file must have a fastq suffix otherwise will generate a weird output file name."""
    filein = open(fastq_filepath, "r")
    fileout = open(fastq_filepath[:-5] + "fasta", "w")
    found_id = 0
    num_of_seqs = 0
    for i in filein:
        if i[0] == "@":
            seq_id = ">" + i[1:]
            found_id = 1
            num_of_seqs += 1
            continue
        if found_id == 1:
            seq = i
            found_id = 0
            fileout.write(seq_id + seq)
    filein.close()
    fileout.close()
    print num_of_seqs
    return os.path.abspath(fileout.name)

if __name__ == '__main__':
    fast_Q2A(sys.argv[1])
