import sys
import numpy 
import cmath 
import math 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import statistics
import scipy.sparse as sps

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
  DoS1.append(float(data[0])/(Nmax)/2)
  DoS2.append(-float(data[1])/(Nmax))

fh.close
yex=[]
yex.append(max(DoS2)+0.01)
yex.append(min(DoS2)-0.01)
xmax=max(DoS1)+0.1
xmin=min(DoS1)-0.1

ImgDoS='images/DoS_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(1)
plt.plot(DoS1,DoS2,'b.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.01,yex[0]-0.01])
axes.set_xlim(xmin,xmax)
if f<1:
  xdos,CDoS=DoSR1(xmin,xmax,e1,e2,omega,omega0,eta,gamma)
if f>=1 and eta<Delta:
  plt.plot([emin,emin],yex,'r--')
  xdos,CDoS=DoSR2(xmin,xmax,e1,e2,emin,omega,omega0,eta,gamma)
if f>=1 and eta>=Delta:
  plt.plot([emin,emin],yex,'y--')
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
  EiV.append(float(data[0])/(Nmax))

f1.close


Filemjz='results/mjz_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
mJz = []
f1 = open(Filemjz, 'r')
for line in f1:
  data=line.split()
  mJz.append(float(data[0])/(Nmax/2))

f1.close

QI=[]
if nmax>900:
  i=0
  while emin<mJz[i]:
    i=i+1
  QI.append(i)
  i=0
  while eNe<mJz[i]:
    i=i+1
  QI.append(i)
  i=0
  while e1<mJz[i]:
    i=i+1
  QI.append(i) 

  E0=0.01

  while E0<EiV(QI[0]):
    QI[0]=QI[0]+1
else:
  QI.append(0)


yex=[]
yex.append(max(mJz)+0.01)
yex.append(min(mJz)-0.01)
xmax=max(EiV)+0.01
xmin=min(EiV)-0.01

ImgPeresL='images/PeresL_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
plt.figure(2)#,figsize=(3,2))
plt.plot(EiV,mJz,'b.')
plt.plot([e1,e1],yex,'r--')
plt.plot([e2,e2],yex,'r--')
axes=plt.gca()
axes.set_ylim([yex[1]+0.01,yex[0]-0.01])
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

Fileval='results/DJz_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.dat' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
EiVe = []
f1 = open(Fileval, 'r')
for line in f1:
  data=line.split()
  if data:
    data=[complex(i)/(Nmax/2) for i in data]
    EiVe.append(data)

f1.close


for j in QI:
  prob=[]
  hold=0.0
  l=0
  for i in range(0,len(EiVe)):
    prob.append((EiVe[int(i)][QI[j]]*EiVe[int(i)][QI[j]].conjugate()).real)
  rp=range(0,len(prob))
  ImgPeresL='images/Prob_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f_%i.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en,QI[j])
  plt.figure(5)#,figsize=(3,2))
  plt.plot(rp,prob,'b.')
  plt.savefig(ImgPeresL)
  del prob




for l in QI:
  VeB=numpy.empty([int(Nmax/2)-1],dtype=complex)
  hold=0.0
  l=0
  for j in range(0,int(Nmax/2)-1):
    for i in range(0,(nmax)):
      hold=hold +EiVe[int(i+j*(nmax+1))][QI[l]]
    VeB[j]=complex(hold)
    hold=0.0

  a = [ math.sqrt(x) for x in range(len(VeB)-1)]

  am =(sps.spdiags(a,-1, len(VeB),len(VeB))+sps.spdiags(a,1, len(VeB),len(VeB)))#*(sps.spdiags(a,-1, len(VeB),len(VeB))+sps.spdiags(a,1, len(VeB),len(VeB))
  dx= abs((numpy.transpose(VeB).dot((am*am).dot(VeB.conjugate()))).real)-abs((numpy.transpose(VeB).dot((am).dot(VeB.conjugate()))).real)

  am =(sps.spdiags(a,-1, len(VeB),len(VeB))-sps.spdiags(a,1, len(VeB),len(VeB)))
  dy= abs((numpy.transpose(VeB).dot((am*am).dot(VeB.conjugate()))).real)-abs((numpy.transpose(VeB).dot((am).dot(VeB.conjugate()))).real)

  del am
  del a


  if dx>dy:
    h=4*dx/100
    d=dx
  else:
    h=4*dy/100
    d=dy



  xrang= numpy.arange(statistics.median(VeB).real-4*d,statistics.median(VeB).real+4*d,h)
  yrang= numpy.arange((statistics.median(VeB)).imag-4*d,(statistics.median(VeB)).imag+4*d,h)

  rho=numpy.outer(VeB,VeB)
  nrho=len(rho)

  Qfun=[[0 for i in xrange(len(xrang))] for i in xrange(len(yrang))]

  alpha=numpy.empty([nrho],dtype=complex)
  h=0
  g=0
  for i in xrang:
    for l in yrang:
      for k in range(0,nrho):
        alpha[k]=cmath.exp(-pow(abs(i+1j*l),2)/2)*(pow((i+1j*l),k)/math.sqrt(math.factorial(k)))
      
      Qfun[g][h]=(numpy.dot(numpy.transpose(alpha),numpy.dot(rho,alpha.conjugate()))).real
      g=g+1
      
    h=h+1
    g=0
  
  
  X, Y = numpy.meshgrid(xrang, yrang)

  QfunName='images/Qfun_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
  fig = plt.figure(3)#,figsize=(3,2))
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X,Y,Qfun)
  plt.matshow(Qfun)
  plt.savefig(QfunName)

  QfunName3d='images/Qfun3D_%d_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f.eps' % (Nmax,nmax,omega,omega0,Delta,eta,gamma,en)
  fig = plt.figure(4)#,figsize=(3,2))
  ax = fig.add_subplot(111, projection='3d')
  ax.plot_surface(X,Y,Qfun)
  plt.savefig(QfunName3d)

  del xrang
  del yrang
  del rho
  del alpha
  del Qfun
  del X
  del Y
  del VeB






########################################################################################



##################################################################################
##################################################################################
