
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
  float tol;
  double Delta, eta, gamma, omega, omega0, alpha,f;
  int Nmax, nmax,size;
  unsigned int i,j,k,l;//counters
  
  // initializing variables 
  Nmax=6; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;
  tol=0.000001;
  alpha = 2*gamma/(omega*sqrt(Nmax));
  size= int(Nmax/2)*(nmax+1)+int(nmax/2+1);
  sp_cx_mat H(size,size);//,jp,jp0,dJz;
  sp_cx_mat jp(size,size);
  sp_cx_mat jp0(size,size);
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
     H(i*int(nmax+1)+j+int(nmax/2)+1,i*int(nmax+1)+j+int(nmax/2)+1)=omega*j-omega*alpha*alpha*(j+1)*(j+1);
   }
  }
  /*tm.set_size(int(Nmax/2)*(nmax+1),int(Nmax/2)*(nmax+1));//move this to the end only used for post processing
  for(i=0;i<Nmax/2-1;i++){
   for(j=0;j<nmax;j++){
     tm(i*int(nmax)+j,i*int(nmax)+j+1)=(i+1)*sqrt(j);
     tm(i*int(nmax)+j+1,i*int(nmax)+j)=(i+1)*sqrt(j+1);
   }
  }*/

cout<<H<<endl;
  for(l=0;l<Nmax;l++){
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	 if((i+(l+1)*int(Nmax)<size)&&(j+l*int(Nmax)<size)){
	  jp(i+(l+1)*int(Nmax),j+l*int(Nmax))=0;
	  for(k=0;k<min(i,j)+1;k++){
	    jp(i+(l+1)*int(Nmax),j+l*int(Nmax))+=exp(complex<double>((i+j-2*k)*log(alpha))+complex<double>(j-k)*log(complex<double>(-1))+complex<double>(0.5*(lgamma(i+1)+lgamma(j+1))-(lgamma(i-k+1)+lgamma(j-k+1)+lgamma(k+1)))); 
	  }
	  jp(i+(l+1)*int(Nmax),j+l*int(Nmax))*=complex<double>(exp(-alpha*alpha/2))*sqrt(complex<double>(Nmax/2*(Nmax/2+1)-double(l*(l-1))));
	 }
  }}}

  for(i=0;i<nmax+1;i++){
    for(j=0;j<(nmax+1)/2;j++){
      for(k=0;k<min(i,j)+1;k++){
	  jp0(i+int(Nmax),j)+=exp(complex<double>((i+j-2*k)*log(alpha))+complex<double>(j-k)*log(complex<double>(-1))+complex<double>(0.5*(lgamma(i+1)+lgamma(j+1))-(lgamma(i-k+1)+lgamma(j-k+1)+lgamma(k+1))));
	}
	jp0(i+int(Nmax),j)*=complex<double>(exp(-alpha*alpha/2)*sqrt((Nmax/2*(Nmax/2+1))/ 2));
  }}
  
  dJz=-(2*jp0.t()+2*jp0+jp+jp.t())/2;
  H=H+Delta*dJz+eta/Nmax*dJz*dJz;
  cout<<H<<endl;
  cx_vec eigval;
  cx_mat eigvac;
  eigs_gen(eigval, eigvac, H,int(10),"sr");
  
  eigval=2*eigval/Nmax;
  ofstream fileeva("eigenval.dat");
  fileeva << eigval<<endl;
  fileeva.close();
  
  ofstream fileeve("eigenvec.dat");
  for(i=0;i<eigvac.n_cols;i++){
    for(j=0;j<eigvac.n_rows;j++){
      fileeve << conj(eigvac(j,i))*eigvac(j,i)<< " ";
    }
    fileeve << "\n";
  }
  fileeve.close();
  
}