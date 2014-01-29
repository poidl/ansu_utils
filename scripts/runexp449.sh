P='/home/z3439823/mymatlab/omega/ansu_utils/'
for i in {450..548}
do
	matlab -nodesktop -nosplash -r "run /home/z3439823/mymatlab/omega/ansu_utils/exp$i/run.m; quit"
done
