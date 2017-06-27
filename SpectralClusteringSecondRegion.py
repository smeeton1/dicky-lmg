# Basic numpy, complex math, special functions and linear algebra
import numpy as np
import cmath as cm
import math as nm
import scipy.linalg as spl

#Sparse matrices and linear algebra
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

#import scipy.misc as spm
import scipy.special as spsp

#Measuring time
import time

#Parallel processing
from multiprocessing import Pool

#itertools
import itertools 

#Starts clocks
tic = time.clock()

ticp = time.clock()

# System parameters
Nmax = 100;  #qubit ensemble dimension only even numbers
nmax = 4*Nmax;  #field dimension only even numbers

omega  = 1
omega0 = 1

Delta  = 1
gamma  = 0.6

alpha  = 2*gamma/(omega*np.sqrt(Nmax))

# Creates Suro's matrices as an exact copy from his notebook
n0Elem = [ 2*x for x in range(np.int(nmax/2)+1)]
n0 = sps.coo_matrix(sps.spdiags(n0Elem,0,np.int(nmax/2)+1,np.int(nmax/2)+1))
del n0Elem

Ic = sps.coo_matrix(sps.identity(np.int(Nmax/2)))

nElem = [x for x in range(nmax + 1)]
n = sps.coo_matrix(sps.spdiags(nElem,0,nmax + 1,nmax + 1))
del nElem

tn = sps.coo_matrix(sps.block_diag( (n0, sps.kron(Ic,n))))
del n0
del Ic

Ia0 = sps.coo_matrix(sps.bsr_matrix((np.int(nmax/2)+1,np.int(nmax/2)+1)))

Ia = sps.coo_matrix(sps.identity(nmax+1))


m2Elem = [ (x+1)*(x+1) for x in range( np.int(Nmax/2) )]
m2 = sps.coo_matrix(sps.spdiags(m2Elem,0, np.int(Nmax/2),np.int(Nmax/2)))
del m2Elem

dJx2 = sps.coo_matrix(sps.block_diag( (Ia0, sps.kron(m2,Ia))))
del Ia

adElem = [ np.sqrt(x+1) for x in range(nmax +1) ]
ad = sps.coo_matrix(sps.spdiags(adElem,-1, nmax +1, nmax+1))
del adElem

aElem = [ np.sqrt(x) for x in range(nmax +1) ]
a = sps.coo_matrix(sps.spdiags(aElem,1, nmax +1, nmax+1))
del aElem

p = ad + a
del ad
del a

mElem = [ x+1 for x in range( np.int(Nmax/2) )]
m = sps.coo_matrix(sps.spdiags(mElem,0,np.int(Nmax/2),np.int(Nmax/2)))
del mElem

tm = sps.coo_matrix(sps.block_diag( (Ia0, sps.kron(m,p))))
del Ia0
del p

jpElem = [np.sqrt( Nmax/2 *( Nmax/2 + 1) - x*(x - 1)) for x in range(2,np.int(Nmax/2)+1)]
jp = sps.coo_matrix(sps.block_diag( ([[0]], sps.spdiags(jpElem,-1,np.int(Nmax/2),np.int(Nmax/2)))))
del jpElem

jp0 = sps.coo_matrix( ( [np.sqrt( (Nmax/2*(Nmax/2 +1))/ 2)] ,([1],[0])), shape=(np.int(Nmax/2)+1,np.int(Nmax/2)+1)) 

tocp = time.clock()

print "I've just finished defining the simple matrices. It took me: ", tocp-ticp

ticp = time.clock()

#Defines the fock to coherent basis projection in terms of logarithms and back with exponentials
def _func(ip):
		i, j, k = ip
		return cm.exp( (i+j-2*k)*np.log(alpha) + (j-k)*cm.log(-1) + 0.5*( nm.lgamma(i+1) + nm.lgamma(j+1) ) - ( nm.lgamma(i-k+1) + nm.lgamma(j-k+1) + nm.lgamma(k+1) ))

def expr(i,j):	
	retList = pool.map(_func, [(i,j,k) for k in range(0, np.minimum(i,j)+1) ] )
	result = np.sum(retList)
	return result

pool = Pool(8)
auxmat= [ [ np.exp(-alpha*alpha/2)*expr(i,j) for j in range(0,nmax+1)] for i in range(0,nmax+1)]
Da  =  sps.coo_matrix(auxmat)
del auxmat

auxmat = [ [ np.exp(-alpha*alpha/2)*expr(i,j) for j in range(0,nmax+1,2)] for i in range(0,nmax+1)]
Da0 =  sps.coo_matrix(auxmat)

pool.close()
pool.join()

del auxmat
del expr

tocp = time.clock()
print "I've just finished defining the matrices that require special functions. It took me: ", tocp-ticp
ticp = time.clock()

tjp  = sps.coo_matrix(sps.kron(jp,Da))
tjp0 = sps.coo_matrix(sps.kron(jp0,Da0))
tjm  = sps.coo_matrix.transpose(tjp)
tjm0 = sps.coo_matrix.transpose(tjp0)

del jp
del jp0
del Da
del Da0


auxdata= tjm0.data
auxcol = tjm0.col - nmax/2
auxrow = tjm0.row
auxshape  = tn.shape
tjm0ca = sps.coo_matrix( (auxdata,(auxrow,auxcol)), shape=auxshape)
tjp0ca = sps.coo_matrix.transpose(tjm0ca)

del auxdata
del auxcol
del auxrow
del auxshape
del tjp0
del tjm0

auxdata= tjp.data
auxcol = tjp.col - nmax/2
auxrow = tjp.row - nmax/2
auxshape  = tn.shape
tjpc = sps.coo_matrix( (auxdata,(auxrow,auxcol)), shape=auxshape)
tjmc = sps.coo_matrix.transpose(tjpc)

del auxdata
del auxcol
del auxrow
del auxshape
del tjp
del tjm

Dad1 = tjm0ca + tjpc
Dad2 = tjp0ca + tjmc
Dad3 = tjm0ca + tjp0ca

del tjm0ca
del tjp0ca
del tjpc
del tjmc

dJz = -(Dad1 + Dad2 + Dad3) / 2
dn = tn + alpha*alpha*dJx2 - alpha*tm

del Dad1
del Dad2
del Dad3


dJz2 = dJz.dot(dJz)

tocp = time.clock()
print "I'm ready to enter the loop. It took me: ", tocp-ticp
ticp = time.clock()


def myRange(start, end, step):
	while start<=end:
		yield start
		start += step

etas  = np.zeros(shape=((0.91-0.3)/0.01, ))
en1s  = np.zeros(shape=((0.91-0.3)/0.01, ))
en2s  = np.zeros(shape=((0.91-0.3)/0.01,))
en3s = np.zeros(shape=((0.91-0.3)/0.01, ))
cen1 = np.zeros(shape=((0.91-0.3)/0.01,20))
cen2 = np.zeros(shape=((0.91-0.3)/0.01,20))
cen3 = np.zeros(shape=((0.91-0.3)/0.01,20))

j=0
for eta in myRange(0.3, 0.91, 0.01):
	
	etas[j] = eta

	H = sps.coo_matrix( omega*tn - (omega*alpha*alpha)*dJx2 + Delta*dJz + (eta/Nmax)*dJz2 )
	
	tocp = time.clock()
	print "I just finished defining the Hamiltonian. It took me: ", tocp-ticp
	ticp = time.clock()
	
	f = (4*gamma*gamma + eta*omega)/(omega*omega0)
	en1 = - omega0*Nmax*(f + 1 /f)/4 + eta*Nmax/4
	aeval, aevec = spsla.eigsh(H, k=20, sigma=en1, which='LM', tol = 1e-6)
	
	del aevec
		
	en1s[j] = 2*(en1 )/(omega0*Nmax)
	cen1[j,:] = 2*(aeval)/(omega0*Nmax) 
	
	en2 = -omega0*Nmax/2 + eta*Nmax/4
	aeval, aevec = spsla.eigsh(H, k=20, sigma=en2, which='LM', tol = 1e-6)
	del aevec
	en2s[j] = 2*(en2 )/(omega0*Nmax)
	cen2[j,:] = 2*(aeval)/(omega0*Nmax) 

	en3 = omega0*Nmax/2 + eta*Nmax/4
	aeval, aevec = spsla.eigsh(H, k=20, sigma=en3, which='LM', tol = 1e-6)
	del aevec
	en3s[j] = 2*(en3 )/(omega0*Nmax)
	cen3[j,:] = 2*(aeval)/(omega0*Nmax) 

	tocp = time.clock()
	print "I just finished calculating the eigenvalues. It took me: ", tocp-ticp
	ticp = time.clock()
	j+=1

#Saving to a file
np.savetxt("etasClusteringFirstRegionNmax100.csv", etas, delimiter=",")
np.savetxt("en1ClusteringFirstRegionNmax100.csv", en1s, delimiter=",")
np.savetxt("en2ClusteringFirstRegionNmax100.csv", en2s, delimiter=",")
np.savetxt("en3ClusteringFirstRegionNmax100.csv", en3s, delimiter=",")
np.savetxt("cen1ClusteringFirstRegionNmax100.csv", cen1, delimiter=",")
np.savetxt("cen2ClusteringFirstRegionNmax100.csv", cen2, delimiter=",")
np.savetxt("cen3ClusteringFirstRegionNmax100.csv", cen3, delimiter=",")

toc = time.clock()
print "I just finished everyting. It took me: ", toc-tic