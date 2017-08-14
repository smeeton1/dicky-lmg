
#include <iostream>
#include <armadillo>

main()
{
  SpMat td;//, hold, n0, I, n;
  float Delta, eta, gamma, omega, omega0, alpha;
  int Nmax, nmax;
  unsigned int i,j;//counters
  
  // initializing variables 
  Nmax=100; //qubit ensemble dimension must be even
  nm=3*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  alpha = 2*gamma/(omega*np.sqrt(Nmax));
  /*----------------------------------------*/
  
  /*I.eye(int(Nmax/2),int(Nmax/2));
  n0.set_size(int(nmax/2)+1,int(nmax/2)+1);
  n.set_size(nmax+1,nmax+1);
  for(i=0;i<nmax/2+1;i++){ n0(i,i)=2*i; }
  for(i=0;i<nmax+1;i++){ n(i,i)=i; }
  hold=kron(I,n);*/
  td.set_size(int(Nmax/2)*(nmax+1)+int(nmax/2)+1,int(Nmax/2)*(nmax+1)+int(nmax/2)+1);
  for(i=0;i<nmax/2+1;i++){ td(i,i)=2*i; }
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;i++){
     td(i*int(Nmax/2)+nmax/2+1+j,i*int(Nmax/2)+nmax/2+1+j)=j;
   }
  }
  
  
}