P='/home/z3439823/mymatlab/omega/ansu_utils/'
CPFROM='exp449'
for i in {450..548}
do
	rsync -avr --exclude-from '/home/z3439823/bin/exclude.txt' $CPFROM/* exp$i
done
