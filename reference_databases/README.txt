human_genome_Ab_parts_masked_bowtie_indices:

The files within this directory are the sequence libraries for tophat2. They correspond to the human reference genome from UCSC, where the Ab encoding regions are masked out.

We removed the fasta file that was the entire human genome with the Ab parts masked out. This file was in 'reference_databases/human_genome_Ab_parts_masked_bowtie_indices'. It makes the tophat mapping step faster, but the file was too large to upload to GitHub.

germline_Ig_ref_seqs:

The files within this directory are the reference fasta files for each germline Ab gene segment. We have included all these sequences for humans, however, reference seqs for other organisms could be added here as well using different subdirectories with the same file naming sceme (i.e. IGHV, IGHD, IGHJ, etc). The human data provided here was gathered by the following: We went to the IMGT/V-QUEST website, then clicked on "IMGT/V-QUEST reference directory sets". This brought us to a bunch of tables of links, and we clicked on the links that corresponded to human Ig sequences. These links (one for each gene class) brought us to a fasta page for all the seqs. Here is the url for the seqs: http://www.imgt.org/vquest/refseqh.html#VQUEST .

igblast_databases:

These are blast search databases (for IgBLAST) of the reference sequences for each of the Ab gene segments. These were made using makeblastdb from stand-alone BLAST, and using the fasta files contained in 'germline_Ig_ref_seqs' as input.
