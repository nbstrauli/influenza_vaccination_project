call_FPCAtest = function(freq_prog_filepath, output_filepath){
	t = read.table(freq_prog_filepath, header=TRUE)
	tpoints = c()
	for (i in colnames(t)){
		i_len = nchar(i, type='chars')
		tpoints = c(tpoints, as.numeric(substr(i, 2, i_len)))
	}
	m = as.matrix(t)
	
	#install dependencies
	if ('sp' %in% rownames(installed.packages()) == FALSE){
		install.packages('sp')
	}
	if ('akima' %in% rownames(installed.packages()) == FALSE){
		install.packages('akima')
	}
	
	#run the FPCA-based test created by Wu and Wu, 2013
	source(file='./FPCAtest.R')
	output = FPCAtest(m, tpoints, rep(1, length(tpoints)))
	#adjust p values using Benjamini Hochberg method
	adj_pvals = p.adjust(output$p, method="BH")
	
	#write results of the test to output file
	d = data.frame(gene_names=rownames(t), p_values_foreach_gene=adj_pvals, F_stats_foreach_gene=output$stat)
	d = d[order(d$p_values_foreach_gene),]
	write.table(d, file=output_filepath, sep='\t', row.names=FALSE, quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
call_FPCAtest(args[1], args[2])
