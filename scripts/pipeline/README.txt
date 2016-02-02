All the scripts within this directory together run the bioinformatic pipeline that isolates Ab sequences from RNAseq data. There are many scripts within this directory, however, the script 'pipeline.bash' will call all of these scripts in the appropriate order. Below is a description of 'pipeline.bash', its input data, its output data, and its dependencies. For descriptions of each of the other scripts within this directory, open them in the editor of your choice and read the descriptions that are included with the source code (it should be pretty comprehensive).

#########################

pipeline.bash:

This script essentially calls a whole bunch of other scripts/functions that carry out the steps of finding Ab encoding reads in the multitude of reads from RNAseq data. It starts out by using TopHat2 to align the RNAseq reads to a masked reference genome, where the loci that encode for Ab genes are masked out of this reference. It then uses IgBLAST to isolate the Ab encoding reads from the reads that DON'T map in the previous step. It then parses the information in the IgBLAST alignments to find the overall Ab expression, the Ab diversity (using CDR3 sequence data), the expression level for each Ab gene segment, and the expression level for each CDR3 sequence, for each time-point. Using IgBLAST to align each potential Ab encoding read is memory intensive, and would likely take too long on a personal computer. Because of this, we utilize the Sun Grid Engine (SGE) software with a computational cluster to distribute many of the steps in this pipeline across many computational nodes in order to speed things up. This pipeline will not work without SGE, and most likely will take too long without a computational cluster. Although one could alter the code with relative ease to get rid of the components that utilize SGE, should they have a sufficiently small amount of data.

Usage

>bash pipeline.bash [path_to_input_directory] [path_to_output_directory] [path_to_IgBLAST_directory]

'path_to_input_directory' - This is the path to the directory that contains the input fastq files. Each fastq file should represent all the RNAseq data for one time-point. The names of each fastq file should include their time-point information as follows: the end of the file name should end in '_X', where X is an integer that represents the time-point. So if the time-point for a given fastq file is day 77 then the name of this file should look something like this: 'blah_77.fastq'. All of the downstream analysis of the tools in this repository depend on time-series information. That is why time-series information is assumed in the input data. However, if one merely wants to get Ab information from their RNAseq data, and is not interested in the other aspects of this repository, then they can just arbitrarily name their files using the described format above and 'pipeline.bash' should run fine.

'path_to_output_directory' - This is the path to the directory that all output from 'pipeline.bash' will be written to. Their is quite a bit of temporary files that are created while pipeline.bash is running. These files can be quite large. They have a prefix of 'temp_'. Make sure that you have quite a bit of space available (i.e. 100's of GB's) available for space in your machine. The output from all the IgBLAST alignments are saved to disk temporarily, and they can take up quite a bit of space.

'path_to_IgBLAST_directory' - This is the path to the directory that contains IgBLAST in your environment. For some reason it is difficult to get IgBLAST to work on the command line without actually navigating to the directory that contains the program. This directory must contain 'igblastn' as well as the 'internal_data' directory that must accompany it in order for it to work.

Output

pipeline.bash creates many output files, most of which are 'temp_' files. We have divided up these output files based upon whether they are temporary or permanent (i.e. not deleted at the end of the script). Here are descriptions of all of them:

permanent files:

	Ig_expression_overtime - This file gives the overall Ig expression found foreach time-point. It is a tab delimited file. The first column gives the time-point. The second column gives the number of reads that successfully mapped to the Ab-masked genome using tophat, the third column gives the number of reads successfully mapped using IgBLAST, and the final column gives the ratio of the 2nd and 3rd columns.
	
	cdr3_seqs_and_counts - This directory contains the nucleotide sequences for each unique CDR3 seq that was found, accompanied by the number of times that sequence was found, for each time-point. Each file within this directory represents a time-point. The files are tab delimited, where each line represents a unique CDR3 sequence, the first column is the nucleotide sequence of that CDR3, and the 2nd column lists the number of reads that this CDR3 sequence was found in.
	
	gene_usage - This file contains the expression level for each gene that was detected in the IgBLAST alignments, and at each time-point. It has this information for each class of gene segments. Each file in this directory contains the gene usage information for a given gene class. The gene usage files are tab delimited, where the first column lists the names for each of the genes that were detected, the 2nd column lists the expression level for each of these genes for the 1st time-point, the 3rd column does the same for the 2nd time-point, etc. The genes are ordered by their absolute range in expression level over the time course, and allele information for the genes has been removed.
	
	cdr3_usage - This file contains the expression level for each unique CDR3 sequence that was detected. It is formatted exactly the same way as the files in gene_usage/ except each row represents a CDR3 sequence rather than a gene segment.
	
	cdr3_diversity - This directory contains the diversity information for each time-point. Each file within this directory contains the CDR3 diversity calculated for each downsampling trial. See the manuscript associated with this repository for a description about how we calculated CDR3 diversity (essentially, it is pi, the population genetics diversity statistic). Diversity was calculated for each down sampling trial and then the mean of each of these diversity estimates is listed in the final line.

temporary files:

	temp_fastq_subdir - pipeline.bash works by first separating each input fastq file into many smaller fastq files. It then distributes the task of mapping each of the reads in these smaller fastq files across nodes in a computational cluster. This directory contains each of the smaller fastq files.
	
	temp_mapping_output - This is a directory that contains the output for each of the mapping steps (i.e. using tophat, and then using IgBLAST). It contains bam files from tophat and alignment files from IgBLAST. It also contains the information that is parsed from the IgBLAST alignments. This includes CDR3 sequences, germline gene identity, and e values. Within temp_mapping_output/ , there are directories that contains all the information listed above for each of these smaller fastq files.
	
	temp_gene_names - This directory contains the gene segment names for each class of gene segments and for each time-point. This directory is structured as follows: each subdirectory within this directory is named after one of the time-points, within each time-point directory contain the gene segment name information for all the reads within that time-point. This includes all 7 classes of gene segments (IGHV, IGHD, IGHJ, IGLV, IGLJ, IGKV, and IGKJ). The files within these subdirectories are formatted such that each line represents one of the reads that WEREN'T able to be mapped using tophat. If the read successfully mapped to one of the gene segments of that class using IgBLAST then the name of that gene will be listed. If the read did not successfully map to any genes for a particular class then 'N/A' will be listed. Thus, there is one line for each read that was attempted to be aligned using IgBLAST in these files, and each line has either a gene name (if the read successfully mapped to a given gene segment) or 'N/A' (if the read did not map to any gene segments for a given class).
	
	temp_evalues - This directory gives the e values that correspond to each alignment attempted with IgBLAST. It gives the e values associated with the alignments for each of the different gene segment classes. The directory is structured the same way as temp_gene_names/ is. The files are also tab delimited, where each line represents a read that was given to IgBLAST for alignment. The 1st column is the number of the read, and the 2nd column is the e value associated with that read. The files within a subdirectory contain this information foreach gene segment class. For example the 33rd line in each of these files represent the same read, but give e values for aligning that read to different gene classes.
	
	temp_cdr3_seqs - This directory contains the CDR3 sequence information collected for each time-point. Each file within this directory represents a time-point. Each line in these tab delimited files here represents a CDR3 sequence that was collected from some read. The first column has some unique identifier for the read, and the 2nd column gives the nucleotide sequence of the CDR3 that was found in that read.
	
	temp_downSampled_cdr3_seqsAndCounts - These files are the result of downsampling the CDR3 information. Each subdirectory with in this directory represents a time-point. There are 10 files within each time-point subdirectory. Each of these files represents one trial of downsampling (to account for stochastic effects). They are formatted similarly to cdr3_seqs_and_counts except the CDR3 seq is in the 2nd column and the count is listed first. Also, the last line represents all of the reads for this time-point that did not yield a CDR3 sequence, or were unmapped, and it gives that number. Thus, if you add up all the numbers in the first column they should equal the total number of reads for that time-point.
	
Dependencies:

pipeline.bash is, as its name suggests, a bioinformatic pipeline that strings together a variety of tools that were made in-house as well as by others, so, it has a lot of dependencies. It should work with all of the latest versions of the dependencies listed below as of 12/16/2015. You will want to include the applications in these packages into your PATH, so that they are callable in any directory. That is with the exception of IgBLAST, as for some reason this does not work unless you are in the directory that actually houses IgBLAST and its dependent 'internal_data'. 

TopHat - This is what aligns the RNAseq reads to the masked genome. It can be downloaded here: http://ccb.jhu.edu/software/tophat/index.shtml

Bowtie 2 - TopHat needs Bowtie 2, which can be found here: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/

IgBLAST - A VERY important component of our pipeline! It can be found here: http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone

bedtools - For basic seq-file parsing. This can be found here: http://bedtools.readthedocs.org/en/latest/

SAMtools - To convert BAM files to fastq format. This can be found here: http://samtools.sourceforge.net/

EMBOSS - This is a large package of various bioinformatic applications. We use their 'needle' tool for the global alignment of pairwise sequences. The package can be found at: http://emboss.sourceforge.net/apps/#list
