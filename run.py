import os
import time
import numpy 
import cmath 
import math 
import psutil
import subprocess


for i in range(1,10)
  mem=0
  taken=0.0
  Nmax=2*4#i
  nmax=2*Nmax

  for k in range(0,10)
    begin=time.clock()
    os.system('./Dicky N '+str(Nmax)+' n '+str(nmax))
    p = psutil.Process(os.getpid())
    end=time.clock()
    taken=taken+end-begin
    mem=mem+p.memory_info().rss/(1024.0*1024.0)
    
  Dickytake[i]=taken/10
  Dickymem[i]=mem/10