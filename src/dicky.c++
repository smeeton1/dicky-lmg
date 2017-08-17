
#include <iostream>
#include <complex.h>
#include <armadillo>

main()
{
  sp_cx_mat H,tn,dJx2,tm,jp,jp0,Da,Da0,dJz,ho1,ho2;//, hold, n0, I, n,jp,jp0;
  float Delta, eta, gamma, omega, omega0, alpha;
  int Nmax, nmax;
  unsigned int i,j,k;//counters
  
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
  tn.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//here we are setting up the matrix for a^dagger a
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;j++){
     tn(i*int(Nmax/2)+j,i*int(Nmax/2)+j)=j;
   }
  }
  
  
  dJx2.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//here we are setting up the matrix (a+a^dagger)Jx
  for(i=0;i<nmax+1;i++){
   for(j=0;j<Nmax/2;j++){
     dJx2(i*int(Nmax/2)+j,i*int(Nmax/2)+j)=(j+1)*(j+1);
   }
  }
  
  tm.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//move this to the end only used for post processing
  for(i=0;i<Nmax/2-1;i++){
   for(j=0;j<nmax;j++){
     tm(i*int(nmax)+j,i*int(nmax)+j+1)=(i+1)*sqrt(j);
     tm(i*int(nmax)+j+1,i*int(nmax)+j)=(i+1)*sqrt(j+1);
   }
  }
  //tm(int(Nmax/2)*(nmax+1)+nmax/2+1,int(Nmax/2)*(nmax+1)+nmax/2)=(Nmax/2)*sqrt(nmax+1);
  //tm(int(Nmax/2)*(nmax+1)+nmax/2,int(Nmax/2)*(nmax+1)+nmax/2+1)=(Nmax/2+1)*sqrt(nmax+2);
  
  
   /* following is the setting up pf the matrix for Jz*/
  jp.set_size(int(Nmax/2),int(Nmax/2));
  for(i=0;i<Nmax/2;i++){ jp(i+1,i)= sqrt( Nmax/2 *( Nmax/2 + 1) - i*(i - 1)); }
  jp0.set_size(int(Nmax/2),int(Nmax/2));
  jp0(1,0)=sqrt( (Nmax/2*(Nmax/2 +1))/ 2;
  
  Da.set_size(nmax+1,nmax+1); 
  for(i=0;i<nmax+1;i++){
    for(j=0;j<nmax+1;j++){
      Da(i,j)=0;
      for(k=0;k<min(i,j)+1;k++){
	Da(i,j)+=cexp((i+j-2*k)*log(alpha) + (j-k)*clog(-1) + 0.5*(lgamma(i+1) + lgamma(j+1))-(lgamma(i-k+1) + lgamma(j-k+1) + lgamma(k+1)));
  }}}
  
  Da0.set_size(nmax+1,nmax+1);
  for(i=0;i<nmax+1;i++){
    for(j=0;j<(nmax+1)/2;j++){
      Da0(i,j)=0;
      for(k=0;k<min(i,j)+1;k++){
	Da0(i,j)+=cexp((i+j*2-2*k)*log(alpha) + (j*2-k)*clog(-1) + 0.5*(lgamma(i+1) + lgamma(j*2+1))-(lgamma(i-k+1) + lgamma(j*2-k+1) + lgamma(k+1)));
  }}}
  
  ho1=kron(jp,Da);
  ho2=kron(jp0,Da0);
  
  dJz=-(2*ho2.t+2*ho2+ho1+ho1.t)/2;
  
  
  H=omega*tn-omega*alpha*alpha*dJx2+d*tJz+eta/Nmax*dJz*dJz;
  
}