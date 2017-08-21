
#include <iostream>
#include <fstream>
#include <complex.h>
#include <armadillo>


using namespace std;
using namespace arma;

main()
{
  //sp_cx_mat H,dJz;//, hold, n0, I, n,jp,jp0,tn,dJx2,tm;
  //cx_mat Da;
  float Delta, eta, gamma, omega, omega0, alpha;
  int Nmax, nmax;
  complex<double> hold;
  unsigned int i,j,k,l;//counters
  
  // initializing variables 
  Nmax=10; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  alpha = 2*gamma/(omega*sqrt(Nmax));
  /*----------------------------------------*/
  sp_cx_mat H(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//,jp,jp0,dJz;
  sp_cx_mat dJz(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  //here we are setting up the matrix for a^dagger a and (a+a^dagger)Jx
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
  
 /* tm.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//move this to the end only used for post processing
  for(i=0;i<Nmax/2-1;i++){
   for(j=0;j<nmax;j++){
     tm(i*int(nmax)+j,i*int(nmax)+j+1)=(i+1)*sqrt(j);
     tm(i*int(nmax)+j+1,i*int(nmax)+j)=(i+1)*sqrt(j+1);
   }
  }*/
  //tm(int(Nmax/2)*(nmax+1)+nmax/2+1,int(Nmax/2)*(nmax+1)+nmax/2)=(Nmax/2)*sqrt(nmax+1);
  //tm(int(Nmax/2)*(nmax+1)+nmax/2,int(Nmax/2)*(nmax+1)+nmax/2+1)=(Nmax/2+1)*sqrt(nmax+2);
  
  /* following is the setting up pf the matrix for Jz*/
  
  dJz.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  //for(i=0;i<Nmax/2-1;i++){ jp(i+1,i)= sqrt( Nmax/2 *( Nmax/2 + 1)-(-Nmax/2+i)*(-Nmax/2 +i - 1)); }
  //jp0.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));
  /*jp0(i,i+1)=sqrt((Nmax/2*(Nmax/2 +1)))/ 2;
  for(i=0;i<Nmax/2-1;i++){ jp0(i,i+1)= sqrt( Nmax/2 *( Nmax/2 + 1)-(-Nmax/2+i)*(-Nmax/2 +i - 1)); }*/
  
  //Da.set_size(nmax+1,nmax+1); 
  for(l=0;l<Nmax/2;l++){ 
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	dJz(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))=0;
	dJz(j+l*int(Nmax/2),i+(l+1)*int(Nmax/2))=0;
	for(k=0;k<min(i,j)+1;k++){
	  hold=cexp((i+j-2*k)*log(alpha) + (j-k)*log(-1) + 0.5*(lgamma(i+1) + lgamma(j+1))-(lgamma(i-k+1) + lgamma(j-k+1) + lgamma(k+1)));
	  dJz(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))+=hold;
	  dJz(j+l*int(Nmax/2),i+(l+1)*int(Nmax/2))+=hold;
	}
	hold=sqrt(Nmax/2*(Nmax/2+1)-(-Nmax/2+l)*(-Nmax/2 +l-1));
	dJz(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))*=-hold/complex< double >(2,0);
	dJz(j+l*int(Nmax/2),i+(l+1)*int(Nmax/2))*=-hold/complex< double >(2,0);
      }      
    }
  }

  //dJz=-(jp+jp0);
  
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