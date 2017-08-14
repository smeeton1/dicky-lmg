
#include <iostream>
#include <armadillo>

main()
{
  SpMat n0, I, n;
  float Delta, eta, gamma, omega, omega0, alpha;
  int Nmax, nmax;
  int i,j;//counters
  
  // initializing variables 
  Nmax=100; //qubit ensemble dimension must be even
  nm=3*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  alpha = 2*gamma/(omega*np.sqrt(Nmax));
  /*----------------------------------------*/
  
  I.set_size(int(Nmax/2),int(Nmax/2));
  n0.set_size(int(nmax/2)+1,int(nmax/2)+1);
  n.set_size(int(nmax/2)+1,int(nmax/2)+1);
  for(i=0;i<Nmax/2;i++){ I(i,i)=1; }
  for(i=0;i<nmax/2+1;i++){ n0(i,i)=i*i; n(i,i)=i; }
  
}