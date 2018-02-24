import sys
import numpy 
import cmath 
import math 
import matplotlib.pyplot as plt
import scipy.integrate as integrate

########################################################################################################

##############################################################################
################# Phi Function ###############

def PhiF(a,omega,omega0,eta,gamma,i):
  h1=((omega*omega0)*((eta/(2*omega0))*pow(a,2)+a-i))/(2*pow(gamma,2)*(1-pow(a,2)))
  if h1<0:
    h1=0
  h1=math.sqrt(h1)
  if h1>1:
    h1=1
  return numpy.arccos(h1)

#####################################################################################


##### Simple Numaricla intergration ################

######## need zero check the sqrt is crashing when near zero

def IntPhi(b,a,omega,omega0,eta,gamma,i):
  if b>a:
    phi=(b-a)*(PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2
  else:
    phi=(a-b)*(PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2
  return phi


def IntPhi2(b,a,omega,omega0,eta,gamma,i):
  if b>a:
    phi=(b-a)/3*((PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2+PhiF(a+((b-a)/3),omega,omega0,eta,gamma,i)+PhiF(a+2*((b-a)/3),omega,omega0,eta,gamma,i)+PhiF(a+3*((b-a)/3),omega,omega0,eta,gamma,i))
  else:
    phi=(a-b)/3*((PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2+PhiF(b+((a-b)/3),omega,omega0,eta,gamma,i)+PhiF(b+2*((a-b)/3),omega,omega0,eta,gamma,i)+PhiF(b+3*((a-b)/3),omega,omega0,eta,gamma,i))
  return phi

def IntPhiN(b,a,omega,omega0,eta,gamma,i,N):
  j=1
  if b>a:
    phi=(b-a)/N*(PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2
    while j<N+1:
        phi=phi+(b-a)/N*PhiF(a+j*((b-a)/N),omega,omega0,eta,gamma,i)
        j+=1
  else:
    phi=(a-b)/N*(PhiF(a,omega,omega0,eta,gamma,i)+PhiF(b,omega,omega0,eta,gamma,i))/2
    while j<N+1:
        phi=phi+(a-b)/N*PhiF(b+j*((a-b)/N),omega,omega0,eta,gamma,i)
        j+=1
  return phi

#############################################################################################################

############### classical DoS functions ###################################


########## Region 1 #####################
def DoSR1(xmin,xmax,e1,e2,omega,omega0,eta,gamma):
  CDoS = []
  xdos = []
  i=e1+0.01
  while i<e2:
    z1=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    phi=IntPhi(z2,z1,omega,omega0,eta,gamma,i)
    CDoS.append(0.5*((z2+1)/2+phi/math.pi))
    xdos.append(i)
    i+=0.01
  while e2<=i<xmax:
    CDoS.append(0.5*1)
    xdos.append(i)
    i+=0.01
  return xdos,CDoS

################# Region 2 ####################

def DoSR2(xmin,xmax,e1,e2,emin,omega,omega0,eta,gamma):
  CDoS = []
  xdos = []
  i=emin+0.001
  while i<e1:
    h1=omega0*(omega-2*i*eta)
    if h1<0:
      h1=0
    zn=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    phi=IntPhiN(zp,zn,omega,omega0,eta,gamma,i,10)
    CDoS.append(0.5*(phi/math.pi))
    xdos.append(i)
    i+=0.001
  while e1<=i<e2:
    h1=omega0*(omega-2*i*eta)
    if h1<0:
      h1=0
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(h1)-omega)/eta
    phi=IntPhiN(z2,zp,omega,omega0,eta,gamma,i,100)
    CDoS.append(0.4*((z2+1)/2+phi/math.pi)-0.12)
    xdos.append(i)
    i+=0.001
  while e2<=i<xmax:
    CDoS.append(0.6*(1))
    xdos.append(i)
    i+=0.001
  return xdos,CDoS


#################### Region 3 #########################

def DoSR3(xmin,xmax,e1,e2,emin,eNe,omega,omega0,eta,gamma):
  CDoS = []
  xdos = []
  i=emin+0.001
  while i<eNe:
    h1=omega0*(omega-2*i*eta)
    if h1<0:
      h1=0
    zn=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    phi=IntPhiN(zp,zn,omega,omega0,eta,gamma,i,10)
    CDoS.append((phi/math.pi))
    xdos.append(i)
    i+=0.001
  while eNe<=i<e1:
    h1=omega0*(omega-2*i*eta)
    if h1<0:
      h1=0
    z1=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(h1)-omega)/eta
    zn=-(math.sqrt(h1)+omega)/eta
    phi=IntPhiN(z2,zp,omega,omega0,eta,gamma,i,100)
    phi2=IntPhiN(zn,z1,omega,omega0,eta,gamma,i,100)
    CDoS.append(0.6*((z2-z1)/2+(phi+phi2)/math.pi)-0.29)
    xdos.append(i)
    i+=0.001
  while e1<=i<e2:
    h1=omega0*(omega-2*i*eta)
    if h1<0:
      h1=0
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(h1)-omega)/eta
    phi=IntPhiN(z2,zp,omega,omega0,eta,gamma,i,100)
    CDoS.append(0.4*((z2+1)/2+phi/math.pi)-0.1)
    xdos.append(i)
    i+=0.001
  while e2<=i<xmax:
    CDoS.append(0.4*(1))
    xdos.append(i)
    i+=0.01
  return xdos,CDoS

#################################################################################

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
  
f=(4*gamma*gamma+eta*omega)/omega*omega0
e1=-1+eta/2*Delta
e2=1+eta/2*Delta
eNe=-Delta/2*eta
emin=(-0.5*(f+1/f)+eta/(2*omega0))

############################################################################## add 1 to 1 ratio
### Ploting for the DoS 

FileDoS='results/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
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
xmax=max(DoS1)+0.1
xmin=min(DoS1)-0.1

ImgDoS='images/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(1)#,figsize=(3,2))
plt.plot(DoS1,DoS2,'b.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.9,yex[0]-0.9])
#axes.set_xlim(xmin,xmax)
if f<1:
  xdos,CDoS=DoSR1(xmin,xmax,e1,e2,omega,omega0,eta,gamma)
if f>=1 and eta<Delta:
  plt.plot([emin,emin],yex,'r--')
  xdos,CDoS=DoSR2(xmin,xmax,e1,e2,emin,omega,omega0,eta,gamma)
if f>=1 and eta>=Delta:
  plt.plot([emin,emin],yex,'r--')
  plt.plot([eNe,eNe],yex,'r--')
  xdos,CDoS=DoSR3(xmin,xmax,e1,e2,emin,eNe,omega,omega0,eta,gamma)
plt.plot(xdos,CDoS,'r-')
plt.savefig(ImgDoS)
del DoS1
del DoS2

#######################################################################################
## Ploting Peres lattices

Fileval='results/eigenval_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
EiV = []
f1 = open(Fileval, 'r')
for line in f1:
  data=line.split()
  EiV.append(float(data[0])/(Nmax/2))

f1.close


Filemjz='results/mjz_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
mJz = []
f1 = open(Filemjz, 'r')
for line in f1:
  data=line.split()
  mJz.append(float(data[0])/(Nmax/2))

f1.close

yex=[]
yex.append(max(mJz)+1)
yex.append(min(mJz)-1)
xmax=max(EiV)+0.3
xmin=min(EiV)-0.3

ImgPeresL='images/PeresL_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(2)#,figsize=(3,2))
plt.plot(EiV,mJz,'b.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.995,yex[0]-0.995])
axes.set_xlim(xmin,xmax)
if f>=1 and eta<Delta:
  plt.plot([emin,emin],yex,'r--')
if f>=1 and eta>=Delta:
  plt.plot([emin,emin],yex,'r--')
  plt.plot([eNe,eNe],yex,'r--')
plt.savefig(ImgPeresL)
del EiV
del mJz

#########################################################################################


####### Q-Function ######################################################################

Fileval='results/eigenvec_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
EiVe = [[]]#*(int(int(Nmax/2)*(nmax+1)))
f1 = open(Fileval, 'r')
j=0
for line in f1:
  data=line.split()
  for i in range(0,int(Nmax/2)*(nmax+1)):
    EiVe[j].append(float(data[i])/(Nmax/2))
  j=j+1

f1.close

VeB=[]
hold=0.0
l=1
for j in range(0,int(Nmax/2)):
  for i in range(0,(nmax+1)):
    hold=hold +EiVe[l,i+j*(nmax+1)]
  
VeB.append(hold)

xrang= range((min(VeB)).real-0.1,(max(VeB)).real+0.1,0.1)
yrang= range((min(VeB)).imag-0.1,(max(VeB)).imag+0.1,0.1)
rho=numpy.outer(VeB,VeB)
nrho=len(rho)

Qfun=[]
alpha=[0]*nrho
for i in xrang:
  for l in yrang:
    for k in range(0,nrho):
      alpha[k]=math.exp(-(pow(i,2)+pow(l,2))/2)*(pow((1+1j*l),k)*math.factorial(k))
      
    Qfun[i,j]=real(numpy.dot(alpha,numpy.dot(rho,alpha)))
 
QfunName=='images/Qfun_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(3)#,figsize=(3,2))
plt.matshow(Qfun)
plt.savefig(QfunName)
del VeB
del xrang
del yrang
del rho
del alpha
del Qfun
########################################################################################



##################################################################################
##################################################################################
