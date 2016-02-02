These scripts together run the global gene usage test for convergence. More specifically, they run a statistical test that asks if the patients in the study are using more similar sets of Ab gene-segments to respond to the same stimulus than would be expected by chance. Most of the code, and the main script for this task is in 'convergent_gene_usage_test_global.py'. So, to run this test you will only need to run this script, provided that your input data is correct.

Usage:

>python convergent_gene_usage_test_global.py [path_to_gene_usage_directory] [path_to_FPCAtest_results_directory] [path_to_the_SGS_table] [filepath_for_output] [trials] [FDR]

	path_to_gene_usage_directory - This is the filepath to the directory that contains the gene usage files for each patient. This directory should contain one gene usage file for each individual in the study. A gene usage file is the file that is created by 'pipeline.bash' which lists the gene expression for each gene detected, and at each time-point. Essentially it is a matrix where each row represents a gene, each column represents a time-point, and each entry is the expression level for a given gene/time-point. 

	path_to_FPCAtest_results_directory - This is the directory that contains the results of the FPCA-based test for each patient. Specifically, this is the output 'FPCA_test_results' from 'test_for_differential_expr.py'. It contains one file foreach patient. Each file is tab delimited, where the first column lists the gene names, the 2nd column lists the FDR corrected p-value for a given gene, and the 3rd column gives the F-statistic for that gene.

	path_to_the_SGS_table - This is a table that tells whether or not each gene was found to be significant in each patient. It is also one of the output files of 'test_for_differential_expr.py'. It is a tab delimited table, where each row represents a gene, the 1st column lists that gene's name, the following columns (up to the last one) tell whether or not that gene was significant for a given patient. A gene that is significant is encoded as a '1' and is encoded as a '0' otherwise. The final column gives the sum across the row, or in other words it tells how many patients a given gene was found to be significant in. This value is called "sum of gene significances" or SGS.
	
	filepath_for_output - This is the filepath for which to write the output file to.
	
	trials - Integer. This gives the number of rounds to sample genes from the day 0 gene expression distribution, which is used to construct the null distribution. The more the better.
	
	FDR - Float, OPTIONAL, default 0.05 . This gives the threshold for considering a gene-segment as significantly differentially expressed. Anything above this threshold will not be considered significant.
	
Output:

This test outputs a single relatively small tab delimited table, where it gives the counts for each possible SGS value (for eaxample, if you have 5 individuals in the study than the SGS value could be [0,5]). It gives these counts for the observed data, as well as the null distribution. It also lists the p-value that results of you use a multinomial G-test to compare the observed counts to that of the null.
