import sys
import numpy 
import cmath 
import math 
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
en=0.5

##############################################################################

FileDoS='results/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
DoS1 = []
DoS2 = []
fh= open(FileDoS, 'r')
for line in fh:
  data=line.split()
  DoS1.append(float(data[0]))
  DoS2.append(float(data[1]))

fh.close

ImgDoS='images/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(1)
plt.plot(DoS1,DoS2,'b.')
plt.savefig(ImgDoS)

del DoS1
del DoS2

#######################################################################################


Fileval='results/eigenval_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
EiV = []
f1 = open(Fileval, 'r')
for line in f1:
  data=line.split()
  EiV.append(float(data[0]))

f1.close


Filemjz='results/mjz_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
mJz = []
f1 = open(Filemjz, 'r')
for line in f1:
  data=line.split()
  mJz.append(float(data[0]))

f1.close

ImgPeresL='images/PeresL_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(2)
plt.plot(EiV,mJz,'.')
plt.savefig(ImgPeresL)

del EiV
del mJz

#########################################################################################

##################################################################################
##################################################################################