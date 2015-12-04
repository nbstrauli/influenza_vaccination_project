#!/bin/bash

# $1 = fastq input files. This should be a directory where the only
#      contents are the fastq files resulting from RNAseq. Each file
#      should represent a unique time-point
# $2 = output directory. This is the directory where all output is
#      written (whether temporary or final output)
# $3 = path to the directory that contains IgBLAST as well as the 
#      'internal_data' folder that needs to accompany it

#get the full file paths for the reference databases for
#tophat and igblast
working_dir=$(pwd)
cd ../../reference_databases/human_genome_Ab_parts_masked_bowtie_indices/
masked_genome_database_path=$(pwd)"/UCSC_abparts_masked"
cd ../../reference_databases/igblast_databases/
iglbast_gl_database_path=$(pwd)"/human_gl_"
cd ../../reference_databases/germline_Ig_ref_seqs/human
ig_ref_seqs_path=$(pwd)
cd $working_dir

#partition input fastq files into more manageable fastq subfiles that
#have 1 million reads each.
#While we're at it, get the lowest read count across time-points
mkdir -p $2
echo "Partioning your input fastq files into 1,000,000 read sized chunks. If your input data is large, this may take awhile..."
temp_fastq_subdir="$2temp_fastq_subdir/"
lowest_read_count=$(python partition_fastq_files.py $1 ${temp_fastq_subdir} 1000000)

#map all the reads using tophat and IgBLAST
#This is an array job, where there is one job for each of the files
#created from the above step
file_count=0
for j in $(ls ${temp_fastq_subdir}); do
    file_count=$(( $file_count + 1 ))
done
temp_mapping_output_dir="$2temp_mapping_output/"
mkdir -p ${temp_mapping_output_dir}
#below is the directory that the stdout and error for the 
#'map_reads.bash' jobs will be saved
mkdir -p output_map_reads
output=$(qsub -t 1-${file_count} -o output_map_reads map_reads.bash ${temp_fastq_subdir} ${temp_mapping_output_dir} ${masked_genome_database_path} ${iglbast_gl_database_path} ${ig_ref_seqs_path} $3)
jobID_mapReads=$(echo $output | awk -F '[ .]' '{print $3}')

#combine gene names and evalues and CDR3 seqs within a time point
temp_gene_names_dirpath="$2temp_gene_names"
temp_evalues_dirpath="$2temp_evalues"
temp_cdr3_seqs_dirpath="$2temp_cdr3_seqs"
mkdir -p output_combine
output=$(qsub -hold_jid ${jobID_mapReads} -o output_combine combine_gene_names_evalues_cdr3s_data_foreach_tpoint.py ${temp_mapping_output_dir} ${temp_gene_names_dirpath} ${temp_evalues_dirpath} ${temp_cdr3_seqs_dirpath})
jobID_combine=$(echo $output | awk -F '[ .]' '{print $3}')

#get Ig expression over time
mkdir -p output_IgExpression
output=$(qsub -hold_jid $jobID_mapReads -o output_IgExpression get_Ig_expression_over_time.py ${temp_mapping_output_dir} $2"Ig_expression_overtime" 1e-20)
jobID_expression=$(echo $output | awk -F '[ .]' '{print $3}')

#perform gene usage analysis
mkdir -p output_geneUsage
qsub -t 1-7 -hold_jid ${jobID_combine}","${jobID_expression} -o output_geneUsage get_gene_usage_time_progression.py $2"temp_gene_names" $2"temp_evalues" $2"gene_usage" 1e-20 True ${ig_ref_seqs_path} $2"Ig_expression_overtime" 1000 vgene_heavy dgene_heavy jgene_heavy vgene_lambda jgene_lambda vgene_kappa jgene_kappa

#perform CDR3 usage analysis
mkdir -p output_cdr3Usage
qsub -hold_jid ${jobID_combine}","${jobID_expression} -o output_cdr3Usage get_gene_usage_time_progression.py $2"temp_cdr3_seqs" $2"cdr3_usage" $2"Ig_expression_overtime" 1000 

#format CDR3 sequences into CDR3 seq and count form
mkdir -p output_cdr3Seqs
output=$(qsub -hold_jid ${jobID_combine} -o output_cdr3Seqs make_cdr3_seqAndCount_foreach_tpoint.py $2"temp_cdr3_seqs" $2"cdr3_seqs_and_counts")
echo ${output}
jobID_seqCount=$(echo $output | awk -F '[ .]' '{print $3}')
echo ${jobID_seqCount}

#simulate removing reads (downsample)
echo "lowest read count:"
echo $lowest_read_count
mkdir -p output_simRm
mkdir -p $2"temp_downSampled_cdr3_seqAndCounts"
mkdir -p $2"temp_downSampled_cdr3_seqs"
tpoint_filenames=$(ls $1)
jobIDs_simRm=''
num_tpoints=0
for i in $tpoint_filenames; do
    let 'num_tpoints++'
    #get the time-point id
    tpoint=$(echo $i | awk -F '[_.]' '{print $(NF-1)}')
    in_cdr3_seq_and_count_filepath=$2"cdr3_seqs_and_counts/"${tpoint}
    tpoint_filepath="${1}${i}"
    out_cdr3_seq_and_count_dirpath=$2"temp_downSampled_cdr3_seqAndCounts/"${tpoint}
    mkdir -p $out_cdr3_seq_and_count_dirpath
    out_cdr3_seq_dirpath=$2"temp_downSampled_cdr3_seqs/"${tpoint}
    mkdir -p $out_cdr3_seq_dirpath
    #should remove the portion of this script that writes just the seqs cause it doesn't work
    output=$(qsub -t 1-10 -hold_jid ${jobID_seqCount} -o output_simRm simulate_removing_reads.py $in_cdr3_seq_and_count_filepath $tpoint_filepath $out_cdr3_seq_and_count_dirpath $out_cdr3_seq_dirpath ${lowest_read_count})
    jobID=$(echo $output | awk -F '[ .]' '{print $3}')
    jobIDs_simRm=$jobIDs_simRm",$jobID"
done
jobIDs_simRm=$(echo $jobIDs_simRm | cut -c 2-)

#calculate pairwise genetic distance
mkdir -p output_diversity
qsub -t 1-${num_tpoints} -hold_jid ${jobIDs_simRm} -o output_diversity calc_pairwise_dist.py $2"temp_downSampled_cdr3_seqAndCounts" $2"cdr3_diversity"

exit
