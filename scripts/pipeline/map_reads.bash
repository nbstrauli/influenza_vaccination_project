#!/bin/bash
#$ -S /bin/bash
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

# $1 = directory path that contains fastq sequences
# $2 = directory to write output to
# $3 = path for the reference sequence database
# $4 = path for the germline Ig gene BLAST database
# $5 = path for the germline ref sequences in fasta format
# $6 = path to the directory that contains IgBLAST, as well as 
#      the 'internal_data' directory that needs to accompany it

#determine which chunk of the data we will be dealing
#with based upon the value of $SGE_TASK_ID
file_basename=$(ls $1 | sort | head -$SGE_TASK_ID | tail -1)
output_dir=$2${file_basename}
input_file=$1${file_basename}
mkdir -p ${output_dir}
echo "starting ${file_basename}"

#run tophat
tophat -o ${output_dir} --b2-very-fast $3 ${input_file}
echo "tophat done"

#convert the BAM file of unmapped reads to fastq using bedtools
bamToFastq -i ${output_dir}/unmapped.bam -fq ${output_dir}/unmapped.fastq
echo "fastq of unmapped reads made"

#convert fastq file to fasta file format using very simple 
#in-house parsing script
python fast_Q2A.py ${output_dir}/unmapped.fastq
echo "fasta made"

#run IgBLAST
#IgBLAST is easiest to run when in the directory that contains IgBLAST
working_dir=$(pwd)
cd $6
igblastn -germline_db_V ${4}V -germline_db_J ${4}J -germline_db_D ${4}D -organism human -domain_system imgt -query ${output_dir}/unmapped.fasta -show_translation -outfmt 3 -out ${output_dir}/unmapped_igblast_output -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1
cd ${working_dir}
echo "IgBLAST done"

#parse the alignments from IgBLAST
python parse_igblast_output.py ${output_dir}/unmapped_igblast_output $5
echo "IgBLAST output parsed"

#filter the CDR3 sequences based upon their V gene alignment e value
#the cutoff here is set to 1e-6. Try changing this cutoff to a higher
#value if you want more CDR3 seqs
python filter_cdr3s_by_evalue.py ${output_dir}/unmapped_igblast_output_cdr3_seq ${output_dir}/unmapped_igblast_output_evalues ${output_dir}/unmapped_igblast_output_cdr3_seq_filtered 1e-6
echo "CDR3 sequences filtered"

#get unique CDR3 seqs with counts
#MAY NEED TO REMOVE THIS!
python get_cdr3_seq_dstrb.py ${output_dir}/unmapped_igblast_output_cdr3_seq_filtered
echo "Identical CDR3 sequences pooled and counted"

exit
