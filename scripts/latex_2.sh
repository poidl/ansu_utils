mypath=/home/nfs/z3439823/mymatlab/poster/figures

docs=(error_depth_averaged error_iso_averaged_df stats_pdf_error stats_pdf_values stats_percent_larger)

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


