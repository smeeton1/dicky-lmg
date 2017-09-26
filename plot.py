import os
import time
import numpy 
import cmath 
import math 
import psutil
import subprocess
import matplotlib.pyplot as plt

##################################################################################
##################################################################################

Nmax=10 ##qubit ensemble dimension must be even
nmax=2*Nmax ##field dimension only even numbers 
Delta=1.0
eta=0.2
gamma=0.3
omega=1.0
omega0=1.0
en=0.8

os.system('./cDicky N '+str(Nmax)+' n '+str(nmax)+' W '+str(omega)+' w '+str(omega0)+' D '+str(Delta)+' E '+str(eta)+' G '+str(gamma)+' e '+str(en))

#################################################################################

DoS = []
f1 = open('results/DoS.dat', 'r')
for line in f1
  data=line.split()
  DoS.append(double(data[1]),double(data[2]))

f1.close

plt.figure(1)
plt.plot(DoS,'.','b')
plt.savefig('DoS.eps')

del DoS

#######################################################################################

EiV = []
f1 = open('results/eigenval.dat', 'r')
for line in f1
  data=line.split()
  EiV.append(double(data[1]))

f1.close

mJz = []
f1 = open('results/mjz.dat', 'r')
for line in f1
  data=line.split()
  mJz.append(double(data[1]))

f1.close

plt.figure(1)
plt.plot(EiV,mJz,'.','b')
plt.savefig('PeresL.eps')

del EiV
del mJz

#########################################################################################

##################################################################################
##################################################################################