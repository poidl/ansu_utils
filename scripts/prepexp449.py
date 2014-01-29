#!/usr/bin/python
import os as os
import numpy as np

P='/home/z3439823/mymatlab/omega/ansu_utils/'
cpfrom=450
cpto=cpfrom+99

dirnr=range(cpfrom,cpto+1)

a=np.linspace(0.,2000.,101)
depths=a[1:-1]

#print size(dirnr)
#print size(depths)

for dirn,dep in zip(dirnr, depths):
    s='sed -i s/initial_pressure=1e3/initial_pressure='+str(dep) \
         +'/g '+P+'exp'+str(dirn)+'/params.m'
    #print s
    os.system(s)
