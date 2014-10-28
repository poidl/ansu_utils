function doit {
s1=$'\\begin{figure}
\\centering
\includegraphics[width=0.9\\textwidth]{'

s2=$'}
\\caption{}\label{fig:'

s3=$'}
\\end{figure}'
base=`basename $1`
s4="${base%.*}"
echo "$s1""$1""$s2""$s4""$s3"
}

# 2 is file
doit $1 >>$2

