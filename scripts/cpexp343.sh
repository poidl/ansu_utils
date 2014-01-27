P='/home/z3439823/mymatlab/omega/ansu_utils/'
CPFROM='exp343'
for i in {344..442}
do
	rsync -avr --exclude-from '/home/z3439823/bin/exclude.txt' $CPFROM/* exp$i
done
