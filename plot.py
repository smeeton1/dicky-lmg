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
en=0.8

##############################################################################

DoS1 = []
DoS2 = []
fh= open('results/DoS.dat', 'r')
for line in fh:
  data=line.split()
  DoS1.append(float(data[0]))
  DoS2.append(float(data[1]))

fh.close

plt.figure(1)
plt.plot(DoS1,DoS2,'.')
plt.savefig('DoS.eps')

del DoS1
del DoS2

#######################################################################################

#EiV = []
#f1 = open('results/eigenval.dat', 'r')
#for line in f1:
  #data=line.split()
  #EiV.append(float(data[0]))

#f1.close

#mJz = []
#f1 = open('results/mjz.dat', 'r')
#for line in f1:
  #data=line.split()
  #mJz.append(float(data[0]))

#f1.close

#plt.figure(1)
#plt.plot(EiV,mJz,'.','b')
#plt.savefig('PeresL.eps')

#del EiV
#del mJz

#########################################################################################

##################################################################################
##################################################################################