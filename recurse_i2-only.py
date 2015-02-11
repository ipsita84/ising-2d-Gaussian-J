#!/usr/bin/python3
# vim: set cin ts=4 sw=4 tw=80:
# To run do the following:
# python3 recurse_n_avg.py

import sys, os
from math import *
import numpy as nm

i2dat = 'I2.dat'

lof_i2 = list()

for root, dirs, files in os.walk('./'):
#	print(os.path.join(root, files[1]))
	fpathlist = {}
	for i in files:
		fpathlist[i] = os.path.join(root, i)
	
	if (i2dat in fpathlist):
		lof_i2.append(fpathlist[i2dat])

def get_avg_arr(lof):
	dat  = nm.loadtxt(lof[0])
#	esum = nm.zeros(nm.shape(dat[:,1]))
	col2dat = list()
	for f in lof:
		dat = nm.loadtxt(f)
		col2dat.append(dat[:,1].T)

	mean = nm.mean(col2dat, axis=0, dtype=nm.float64)
	std  = nm.std(col2dat, axis=0, dtype=nm.float64)
	return (dat[:,0], mean, std)

dat, avg, std = get_avg_arr(lof_i2)
res = nm.array([dat, avg, std])
nm.savetxt("I2-avg.dat", res.T, fmt='% .3f\t% .5E\t% .5E',
           header='Beta\t<I2>\tStd.dev.')

print(avg)
print(std)
