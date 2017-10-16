import sys
import numpy 
import cmath 
import math 
import matplotlib.pyplot as plt

##################################################################################
##################################################################################
### takeing the inputs for the graphs or using defults
if len(sys.argv)>1:
  Nmax=int(sys.argv[1]) ##qubit ensemble dimension must be even
  nmax=int(sys.argv[2]) ##field dimension only even numbers 
  omega=float(sys.argv[3])
  omega0=float(sys.argv[4])
  Delta=float(sys.argv[5])
  eta=float(sys.argv[6])
  gamma=float(sys.argv[7])
  en=float(sys.argv[8])
else:
  Nmax=10
  nmax=20
  omega=1.0
  omega0=1.0
  Delta=1.0
  eta=0.2
  gamma=0.3
  en=0.8
  

e1=-1+eta/2*Delta
e2=1+eta/2*Delta
emin=-0.5*((math.pow(4*gamma,2)+eta*Delta)/(omega*Delta)+(omega*Delta)/(math.pow(4*gamma,2)+eta*Delta))+eta/2*Delta
eNe=-Delta/2*eta
f=(4*pow(gamma,2)+eta*omega)/omega*omega0

##############################################################################

FileDoS='results/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
DoS1 = []
DoS2 = []
fh= open(FileDoS, 'r')
for line in fh:
  data=line.split()
  DoS1.append(float(data[0])/(Nmax/2))
  DoS2.append(float(data[1])/(Nmax))

fh.close
yex=[]
yex.append(max(DoS2)+1)
yex.append(min(DoS2)-1)

ImgDoS='images/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(1)
plt.plot(DoS1,DoS2,'b.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.9,yex[0]-0.9])
if f>=1 and eta<Delta:
  plt.plot([emin,emin],yex,'r--')
if f>=1 and eta>=Delta:
  plt.plot([eNe,eNe],yex,'r--')
plt.savefig(ImgDoS)
del DoS1
del DoS2

#######################################################################################


Fileval='results/eigenval_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
EiV = []
f1 = open(Fileval, 'r')
for line in f1:
  data=line.split()
  EiV.append(float(data[0])/(Nmax/2))

f1.close


Filemjz='results/mjz_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
mJz = []
f1 = open(Filemjz, 'r')
for line in f1:
  data=line.split()
  mJz.append(float(data[0])/(Nmax))

f1.close

yex=[]
yex.append(max(mJz)+1)
yex.append(min(mJz)-1)

ImgPeresL='images/PeresL_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(2)
plt.plot(EiV,mJz,'r.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.995,yex[0]-0.995])
if f>=1 and eta<Delta:
  plt.plot([emin,emin],yex,'r--')
if f>=1 and eta>=Delta:
  plt.plot([eNe,eNe],yex,'r--')
plt.savefig(ImgPeresL)
del EiV
del mJz

#########################################################################################

##################################################################################
##################################################################################