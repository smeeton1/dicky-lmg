import os
import time
import numpy 
import cmath 
import math 
import psutil
import subprocess
import matplotlib.pyplot as plt

Dickytake=[]
Dickymem=[]
Dickysize=[]

for i in range(1,10):
  mem=0
  taken=0.0
  Nmax=100*i
  nmax=2*Nmax
  Dickysize.append((Nmax/2)*(nmax+1)+nmax/2+1)
  
  for k in range(0,10):
    begin=time.clock()
    os.system('./Dicky N '+str(Nmax)+' n '+str(nmax))
    p = psutil.Process(os.getpid())
    end=time.clock()
    taken=taken+end-begin
    mem=mem+p.memory_info().rss/(1024.0*1024.0)
    
  Dickytake.append(taken/10)
  Dickymem.append(mem/10)

  
plt.figure(1)
plt.plot(Dickysize,Dickytake,'-')
plt.xlabel('Matrix size')
plt.ylabel('Time')
plt.title('Program run time')
plt.savefig('Time.png')

plt.figure(2)
plt.plot(Dickysize,Dickymem,'-')
plt.xlabel('Matrix size')
plt.ylabel('Memory(Mb)')
plt.title('Program memory use')
plt.savefig('Memory.png')