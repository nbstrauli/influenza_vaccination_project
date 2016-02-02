These scripts run the convergent gene usage test for individual genes. That is, they test to see if a given gene is responding to the same stimulus in more patients than would be expected by chance, given the gene's day 0 expression level in each of the patients. The main script in this directory is 'convergent_gene_usage_test_indiv_genes.py'. You will want to run this script to run the convergence test.

Usage:

python convergent_gene_usage_test_indiv_genes.py [filepath_to_SGS_table] [path_to_gene_usage_directory] [path_for_output_file]

	filepath_to_SGS_table - This is the path to the file that contains the SGS table information. The SGS table is the output that was created by 'test_for_differential_expr.py'. It is a tab delimited table, where each row represents a gene, the 1st column lists that gene's name, the following columns (up to the last one) tell whether or not that gene was significant for a given patient. A gene that is significant is encoded as a 1 and is encoded as a 0 otherwise. The final column gives the sum across the row, or in other words it tells how many patients a given genes was found to be significant in. This value is called "sum of gene significances" or SGS.
	
	path_to_gene_usage_directory - This is the filepath to the directory that contains the gene usage files for each patient. This directory should contain one gene usage file for each individual in the study. A gene usage file is the file that is created by 'pipeline.bash' which lists the gene expression for each gene detected, and at each time-point. Essentially it is a matrix where each row represents a gene, each column represents a time-point, and each entry is the expression level for a given gene/time-point.
	
	path_for_output_file - This is the path for the file that the output of this script will be written to.
	
Output:

convergent_gene_usage_test_indiv_genes.py makes one output file that contains the results of the convergence test for each gene-segment tested. This file is tab delimited, where each row represents a unique gene, the 1st column lists the gene's name, the 2nd column gives the (uncorrected) p-value for each gene, the 3rd column gives the observed SGS value for each gene, the 4th column gives the mode of the null distribution of SGS values (that is, it gives the most likely SGS value under the null hypothesis), and the 5th column gives the sum of the day 0 expression levels across each of the patients for a given gene. The idea with this value, is that if this sum of day 0 expression levels is high, then the gene is most likely highly expressed across the patients, and thus will be likely to be selected by chance in many patients, under the null hypothesis.
