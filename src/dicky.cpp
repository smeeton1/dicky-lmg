
#include <iostream>
#include <fstream>
#include <complex.h>
#include <armadillo>


using namespace std;
using namespace arma;

main()
{

  //,ho1,ho2;//, hold, n0, I, n,tn,dJx2,tm;
  //cx_mat Da,Da0;
  float Delta, eta, gamma, omega, omega0, alpha;
  int Nmax, nmax;
  unsigned int i,j,k,l;//counters
  
  // initializing variables 
  Nmax=10; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  alpha = 2*gamma/(omega*sqrt(Nmax));
  sp_cx_mat H(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//,jp,jp0,dJz;
  sp_cx_mat jp(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  sp_cx_mat jp0(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  sp_cx_mat dJz(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  /*----------------------------------------*/
  //H.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  //here we are setting up the matrix for a^dagger a
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;j++){
     H(i*int(Nmax/2)+j,i*int(Nmax/2)+j)=omega*j-omega*alpha*alpha*(j+1)*(j+1);
   }
  }
  
  
  //dJx2.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//here we are setting up the matrix (a+a^dagger)Jx
//   for(i=0;i<nmax+1;i++){
//    for(j=0;j<Nmax/2;j++){
//      H(i*int(Nmax/2)+j,i*int(Nmax/2)+j)+=-omega*alpha*alpha*(j+1)*(j+1);
//    }
//   }
  
  /*tm.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//move this to the end only used for post processing
  for(i=0;i<Nmax/2-1;i++){
   for(j=0;j<nmax;j++){
     tm(i*int(nmax)+j,i*int(nmax)+j+1)=(i+1)*sqrt(j);
     tm(i*int(nmax)+j+1,i*int(nmax)+j)=(i+1)*sqrt(j+1);
   }
  }*/
  //tm(int(Nmax/2)*(nmax+1)+nmax/2+1,int(Nmax/2)*(nmax+1)+nmax/2)=(Nmax/2)*sqrt(nmax+1);
  //tm(int(Nmax/2)*(nmax+1)+nmax/2,int(Nmax/2)*(nmax+1)+nmax/2+1)=(Nmax/2+1)*sqrt(nmax+2);
  
  
   /* following is the setting up pf the matrix for Jz*/
  //for(i=0;i<Nmax/2;i++){ jp(i+1,i)= sqrt( Nmax/2 *( Nmax/2 + 1) - i*(i - 1)); }
  //jp0(1,0)=sqrt((Nmax/2*(Nmax/2 +1)))/ 2;
  
  //Da.set_size(nmax+1,nmax+1);
  //jp.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  for(l=0;l<Nmax/2;l++){
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	jp(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))=0;
	for(k=0;k<min(i,j)+1;k++){
	  jp(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))+=cexp((i+j-2*k)*log(alpha)+(j-k)*log(-1)+0.5*(lgamma(i+1)+lgamma(j+1))-(lgamma(i-k+1)+lgamma(j-k+1)+lgamma(k+1)));
	}
	jp(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))*=sqrt(Nmax/2*(Nmax/2+1)-l*(l-1));
    }}}
  
  //Da0.set_size(nmax+1,nmax+1);
  //jp0.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  for(i=0;i<nmax+1;i++){
    for(j=0;j<(nmax+1)/2;j++){
      jp0(i+int(Nmax/2),j)=0;
      for(k=0;k<min(i,j)+1;k++){
	jp0(i+int(Nmax/2),j)+=cexp((i+j*2-2*k)*log(alpha)+(j*2-k)*log(-1)+0.5*(lgamma(i+1)+lgamma(j*2+1))-(lgamma(i-k+1)+lgamma(j*2-k+1)+lgamma(k+1)));
      }
      jp0(i+int(Nmax/2),j)*=(sqrt((Nmax/2*(Nmax/2+1)))/ 2);
  }}
  
  //ho1=kron(jp,Da);
  //ho2=kron(jp0,Da0);
  //dJz.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  dJz=-(2*jp0.t()+2*jp0+jp+jp.t())/2;
  
  
  H=H+Delta*dJz+eta/Nmax*dJz*dJz;
  
  cx_vec eigval;
  cx_mat eigvac;
  eigs_gen(eigval, eigvac, H,int(0.75*(Nmax/2)*(nmax+1)));
  
  ofstream fileeva("eigenval.dat");
  fileeva << eigval;
  fileeva.close();
  
  ofstream fileeve("eigenvec.dat");
  for(i=0;i<eigvac.n_rows;i++){
    for(j=0;j<eigvac.n_cols;j++){
      fileeve << conj(eigvac(i,j))*eigvac(i,j)<< " ";
    }
    fileeve << "\n";
  }
  fileeve.close();
  
}