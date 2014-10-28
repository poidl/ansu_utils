base=`basename $1`
fname="${base%.*}"".tex"

touch $fname

cat latex_header.txt >> $fname
cat $1 >> $fname
cat latex_footer.txt >>$fname

pdflatex $fname
rm *aux
rm *log
rm $1
rm *tex

