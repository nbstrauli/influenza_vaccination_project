The scripts in this directory form a wrapper for the FPCA based test for differential gene expression, formed by Wu and Wu, BMC Bioinformatics, 2013. The main script to use here is 'test_for_differential_expr.py'. It is a relatively simple python wrapper that calls the R code (written by Wu et al.) for the FPCA based test, and then processes some of its output. At the request of the authors one must download the R code from their website, which can be found at: http://www.imcportal.org/repository/software/r-code-for-significant-testing-for-time-course-gene-expression-data-using-functional-principal-component-analysis-approaches .

Usage:

>python test_for_differential_expr.py [gene_usage_directory_path] [output_directory_path] [FDR]

	'gene_usage_directory_path' - This is the path to the directory that contains the gene usage files for each patient for a given gene. A gene usage file is the file that is created by 'pipeline.bash' which lists the gene expression for each gene detected, and at each time-point. Essentially it is a matrix where each row represents a gene, each column represents a time-point, and each entry is the expression level for a given gene/time-point. Within this directory there should be exactly one gene usage file for each individual in the study, and each gene usage file should pertain to the same gene segment (ex: all the gene usage files should be for IGHV, or IGKV, etc...).
	
	'output_directory_path' - This is the path to the directory for which all output will be written. It can be whatever you like.
	
	'FDR' - OPTIONAL, default = 0.05. This sets the maximum value that the false discovery rate (FDR) can be for the results of the FPCA based test, in order for a gene to be considered significant. The FPCA based test corrects for multiple testing (by default) using the Benjamini–Hochberg procedure. This sets the upper bound for the resulting FDR's for each gene.
	
Output:

test_for_differential_expr.py does two things, it tests to see which genes are significantly differentially expressed (i.e. responding to some stimulus) in each individual, and then it compiles those results into an "sum of gene significances" (SGS) table.

	FPCA_test_results - This is the directory that contains the results from the FPCA-based test for differential gene expression. It contains one file foreach patient. Each file is tab delimited, where the first column lists the gene names, the 2nd column lists the FDR corrected p-value for a given gene, and the 3rd column gives the F-statistic for that gene.

	SGS_table - This file is the compilation of the results in FPCA_test_results/. It is a tab delimited table, where each row represents a gene, the 1st column lists that gene's name, the following columns (up to the last one) tell whether or not that gene was significant for a given patient. A gene that is significant is encoded as a 1 and is encoded as a 0 otherwise. The final column gives the sum across the row, or in other words it tells how many patients a given genes was found to be significant in. This value is called "sum of gene significances" or SGS.
	
Dependencies:

FPCAtest.R - As mentioned, one must download the FPCA-based test R code from http://www.imcportal.org/repository/software/r-code-for-significant-testing-for-time-course-gene-expression-data-using-functional-principal-component-analysis-approaches . Download this code and save it in this directory as 'FPCAtest.R'.

akima - The FPCA-based test requires this R package.

sp - The FPCA-based test requires this R package.
