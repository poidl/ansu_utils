n=16
#fname='/home/nfs/z3439823/Documents/emails_trevor/email20140821_pics/latex_figures_out.txt'
fname='/home/nfs/z3439823/Documents/emails_trevor/latex_figures_out.txt'

rm $fname

function doit {
s1=$'\\begin{figure}
\\centering
\includegraphics[width=0.9\\textwidth]{\\apgexp dp_omega_omega/'

s2=$'.png}
\\caption{}\label{fig:dp_omega_omega_'

s3=$'}
\\end{figure}'
echo "$s1""$i""$s2""$i""$s3"
}

for i in $(seq -f "%02g" 1 $n);
do
#echo "$s1">> $fname
doit $i >>$fname
done
#\begin{figure}
#\centering
#\includegraphics[width=0.9\textwidth]{\apgexp gins3d/01.png}
#\caption{}\label{fig:gins3d_01}
#\end{figure}
