#!/usr/bin/python3
# vim: set ts=4 sw=4 tw=80:
import matplotlib.pyplot as plt #import this to make plots

import sys, os
from math import *
import numpy as nm

#specify fonts for plots and enable using Latex:
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)


#edat  = 'E.dat'
#emdat = 'Em.dat'
i2dat = 'I2.dat'

#lof_e  = list()
#lof_em = list()
lof_i2 = list()

for root, dirs, files in os.walk('./'):
    fpathlist = {}
    for i in files:
        fpathlist[i] = os.path.join(root, i)

    if (i2dat in fpathlist):
#        lof_e.append (fpathlist[edat] )
#        lof_em.append(fpathlist[emdat])
        lof_i2.append(fpathlist[i2dat])

for filename in lof_i2:
	fin = open(filename, 'r')
	beta_list  = []
	EAvg_list  = []

	labflag=False
	#Loop through all lines of the file:
	for line in fin:
		beta,EAvg = line.split()
		beta_list.append(float(beta))
		EAvg_list.append(float(EAvg))
		if (float(EAvg) < 0):
			labflag=True
		
	if (labflag):
		plt.plot(beta_list, EAvg_list, label=filename)
	else:
		plt.plot(beta_list, EAvg_list)

#	plt.plot(beta_list, EAvg_list)

#Label the plot:
plt.xlabel(r'$\beta$', fontsize=20)              #label for the x-axis
plt.ylabel(r'$\langle E \rangle$', fontsize=20)  #label for the y-axis
#plt.legend(prop={'size':16},loc='lower right')   #info about legend
plt.axis([0, 0.42, plt.axis()[2], plt.axis()[3]])
#Save the plot as a pdf:
plt.savefig('I2VsBeta.pdf',bbox_inches='tight')

#Display the plot:
plt.show()
