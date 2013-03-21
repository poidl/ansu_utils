# -*- coding: utf-8 -*-
import os as os

for srcdir,destdir in zip(range(5,9),range(9,13)):
	os.system('rsync -avr --exclude-from \'/home/z3439823/bin/exclude_cpexp.txt\'\
                   exp'+'{0:03d}'.format(srcdir)+'/ exp'+'{0:03d}'.format(destdir))
#    os.system('him1 '+str(srcdir)+' '+str(destdir))
#    os.system('mkdir run'+str(destdir)+'/saves') 
