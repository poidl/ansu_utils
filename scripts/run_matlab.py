#!/usr/bin/python
import os as os

P='/home/z3439823/mymatlab/omega/ansu_utils/'
for i in range(180,187):
	s='matlab -nodesktop -nojvm -r "run '+P+'exp'+str(i)+'/run.m; quit"'
	print s
	os.system(s)
#	rm -r $P/exp1$j
#	cp -r $P/exp1$i $P/exp1$j
#	rm $P/exp1$j/figures/*
#	rm $P/exp1$j/data/iteration_history.mat

