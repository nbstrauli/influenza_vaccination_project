This directory contains the scripts and data related to the manuscript "Statistical Inference of a Convergent Antibody Repertoire Response to Influenza Vaccine". There are 2 main things that the scripts here do: isolate antibody (Ab) sequences from RNAseq data, and test for a convergent Ab response signal across patients to some antigenic stimulus. The scripts have been divided into 4 submodules that can be used independently, in case one is only interested in a specific aspect of the code.

There are 4 subdirectories contained within that are (meant) to contain different types of information. 

'input_data' - This is an empty directory that you are meant deposite your input RNAseq data into.
'output' - This is also an empty directory that all the output from the scripts are written to.
'reference_databases' - This directory contains all of the reference data and libraries that are necessary for the programs herein.
'scripts' - This directory contains the source code (mostly python) for all of the programs associated with this manuscript (excluding dependencies).

#######################

input_data:
ALL of the data that goes in here should be in fasta format. All the scripts in this module are designed to analyze time series RNAseq data. So, each fasta file should represent one timepoint, and these time-points should be labeled by a '_X' at the end of the file name. Where 'X' is an integer (ordered) identifier of the time-point. So if the fasta file 'shneeky.fasta' represents the RNAseq data for day 43, then its identifier at the end of the file would be '_43' and the full file name should be 'shneeky_43.fasta'. Essentially, the downstream analysis will rely on numerically sorting the integers that come after the last underscore of the file names in this input directory. If one is only interested in getting Ab information from RNAseq data, and (strangely) not interested in any of the time-series analysis downstream, then one does not need to worry about these numerical IDs.

#######################

For detailed explanations of the subdirectories 'output', 'reference_databases', and 'scripts', and the data that they contain, see the README files within each of them.
