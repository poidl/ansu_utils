#!/usr/bin/python
import os as os

P='/home/z3439823/mymatlab/omega/ansu_utils/'
for i,j in zip(range(152,172), range(180,200)):
	s='cp -r '+P+'exp'+str(i)+' '+P+'exp'+str(j)
	print s
	os.system(s)
#	rm -r $P/exp1$j
#	cp -r $P/exp1$i $P/exp1$j
#	rm $P/exp1$j/figures/*
#	rm $P/exp1$j/data/iteration_history.mat

