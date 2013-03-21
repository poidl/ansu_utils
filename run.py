# -*- coding: utf-8 -*-
import os as os

for dirno in range(9,13):
	os.system('matlab -nodesktop  -nosplash -r \"cd omega/ansu_utils/exp'+\
			'{0:03d}'.format(dirno)+'; run; quit\"')

#  	os.system('rsync -avr --exclude-from \'/home/z3439823/bin/exclude_cpexp.txt\'\
#                   exp'+'{0:03d}'.format(srcdir)+'/ exp'+'{0:03d}'.format(destdir))
       
