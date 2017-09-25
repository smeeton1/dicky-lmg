import os
import time
import numpy 
import cmath 
import math 
import psutil
import subprocess
import matplotlib.pyplot as plt



#################################################################################

DoS = []
f1 = open('results/DoS.dat', 'r')
for line in f
  data=line.split()
  DoS.append(double(data[1]),double(data[2]))

f.close

plt.figure(1)
plt.plot(DoS,'.','b')
plt.savefig('DoS.eps')

del Dos

#######################################################################################

EiV = []
f1 = open('results/eigenval.dat', 'r')
for line in f
  data=line.split()
  EiV.append(double(data[1]))

f.close

mJz = []
f1 = open('results/mjz.dat', 'r')
for line in f
  data=line.split()
  mJz.append(double(data[1]))

f.close

plt.figure(1)
plt.plot(EiV,mJz,'.','b')
plt.savefig('PeresL.eps')

del EiV
del mJz

#########################################################################################