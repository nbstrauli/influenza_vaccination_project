# $1 = cdr3 sequence #1
# $2 = cdr3 sequence #2
# $3 = output alignment filepath

PATH="$PATH:/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss"
needle -asequence $1 -bsequence $2 -gapopen 10.0 -gapextend 0.5 -outfile $3 -brief Y -sprotein1 Y -sprotein2 Y
