#!/usr/bin/python3
# vim: set cin ts=4 sw=4 tw=80:
# To run do the following:
# python3 recurse_n_avg.py

import sys, os
from math import *
import numpy as nm

edat  = 'E.dat'
emadat = 'EmA.dat'
embdat = 'EmB.dat'
i2dat = 'I2.dat'

lof_e  = list()
lof_ema = list()
lof_emb = list()
lof_i2 = list()

for root, dirs, files in os.walk('./'):
#	print(os.path.join(root, files[1]))
	fpathlist = {}
	for i in files:
		fpathlist[i] = os.path.join(root, i)
	
	if (edat in fpathlist and emadat in fpathlist and embdat in fpathlist 
    and i2dat in fpathlist):
		lof_e.append (fpathlist[edat] )
		lof_ema.append(fpathlist[emadat])
		lof_emb.append(fpathlist[embdat])
		lof_i2.append(fpathlist[i2dat])

def get_avg_arr(lof):
	dat  = nm.loadtxt(lof[0])
	esum = nm.zeros(nm.shape(dat[:,1]))
	for f in lof:
		dat = nm.loadtxt(f)
		esum = esum + dat[:,1]
		esum = esum / len(lof)
	return (dat[:,0], esum)

dat, avg = get_avg_arr(lof_e)
res = nm.array([dat, avg])
nm.savetxt("E-avg.dat", res.T, fmt='% .3f\t% .5E')

dat, avg = get_avg_arr(lof_ema)
res = nm.array([dat, avg])
nm.savetxt("EmA-avg.dat", res.T, fmt='% .3f\t% .5E')

dat, avg = get_avg_arr(lof_emb)
res = nm.array([dat, avg])
nm.savetxt("EmB-avg.dat", res.T, fmt='% .3f\t% .5E')

dat, avg = get_avg_arr(lof_i2)
res = nm.array([dat, avg])
nm.savetxt("I2-avg.dat", res.T, fmt='% .3f\t% .5E')
