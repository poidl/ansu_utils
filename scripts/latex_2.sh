mypath=/home/nfs/z3439823/mymatlab/poster/figures

docs=(error_depth_averaged error_iso_averaged_df error_iso_averaged_erms stats_pdf_error stats_pdf_values stats_percent_larger instabilities_dep_vs_lat stats_nonpositive_values std_on_omega)

for item in ${docs[*]}
do
	file="$item.txt"
	touch $file
	pdfs=($(ls $mypath/$item*))
	for pdf in ${pdfs[*]}
	do
		#echo $pdf
		./latex_dumpfigure.sh $pdf $file
	done
./latex_txt2tex.sh $file
done


