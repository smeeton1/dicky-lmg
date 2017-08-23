
#include <iostream>
#include <fstream>
#include <complex>
#include <armadillo>
#include <string>
#include <sstream>


using namespace std;
using namespace arma;

main()
{
  double Delta, eta, gamma, omega, omega0, alpha,tol,f;
  int Nmax, nmax,size;
  complex<double> hold;
  unsigned int i,j,k,l;//counters
  
  // initializing variables 
  Nmax=4; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  alpha = 2*gamma/(omega*sqrt(Nmax));
  /*----------------------------------------*/
  size= int(Nmax/2)*(nmax+1)+int(nmax/2+1);
  sp_cx_mat H(size,size);//,jp,jp0,dJz;
  sp_cx_mat dJz(size,size);
  /*----------------------------------------*/
  f = (4*gamma*gamma + eta*omega)/(omega*omega0);
if (f < 1){
	f = - omega0*Nmax/2 + eta*Nmax/4;}
else{
	f = - omega0*Nmax*(f + 1 /f)/4 + eta*Nmax/4;}
	
  //here we are setting up the matrix for a^dagger a
  for(i=0;i<nmax/2+1;i++){H(i,i)=omega*2*i;}
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;j++){
     H(i*int(nmax+1)+j+int(nmax/2)+1,i*int(nmax+1)+j+int(nmax/2)+1)=omega*j-omega*alpha*alpha*(j+1)*(j+1)-f;
   }
  }

  
 /* tm.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//move this to the end only used for post processing
  for(i=0;i<Nmax/2-1;i++){
   for(j=0;j<nmax;j++){
     tm(i*int(nmax)+j,i*int(nmax)+j+1)=(i+1)*sqrt(j);
     tm(i*int(nmax)+j+1,i*int(nmax)+j)=(i+1)*sqrt(j+1);
   }
  }*/

  
  /* following is the setting up pf the matrix for Jz*/
  
  dJz.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));

  for(l=0;l<Nmax/2;l++){ 
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	dJz(i+(l+1)*int(Nmax/2),j+l*int(Nmax/2))=0;
	dJz(j+l*int(Nmax/2),i+(l+1)*int(Nmax/2))=0;
	for(k=0;k<min(i,j)+1;k++){
	  hold=complex<double>(exp((i+j-2*k)*log(alpha)))+complex<double>((j-k))*log(complex<double>(-1))+complex<double>(0.5*(lgamma(i+1)+lgamma(j+1))-(lgamma(i-k+1)+lgamma(j-k+1)+lgamma(k+1)));
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
  
  eigval=(2/Nmax)*eigval;
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