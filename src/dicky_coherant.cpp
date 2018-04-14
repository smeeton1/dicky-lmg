#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <armadillo>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <boost/math/special_functions/factorials.hpp>
#include <math.h>

using namespace std;
using namespace arma;


long double fac(int n){
  int i;
  long double sum=1;
  if(n>0){
  for(i=1;i<=n;i++){sum=sum*i;}
  return sum;}
  else{return 1;}
}

long double fac_st(int n){
  long double sum;
  if(n>0){
  sum=sqrt(2*M_PI*n)*pow(n/M_E,n);
  return sum;}
  else{return 1;}
}

double laguerre(double n,double m,int k){
  double sum,l1,l2,l3;
  int j;
  l1=1;l2=1+m-n;
  sum=l1+l2;
  for(j=1;j++;j=k){
    l3=(m+2*j+1-n)*l2/(j+1)-(j+m)*l1/(j+1);
    l1=l2;l2=l3;
    sum+=l3;
  }
return sum;
}

complex<long double> Kloop(int i,int j,int l,int Nmax,double alpha){
  complex<long double> hold=complex<double>(0, 0);
  long double lhold=0.0;
    #pragma omp parallel for reduction (+:lhold)
    for(int k=0;k<min(i,j)+1;k++){
      lhold=lhold+(pow(alpha,(i+j-2*k))*pow(-1,(j-k))*((sqrt(fac_st(i))*sqrt(fac_st(j)))/(fac_st(i-k)*fac_st(j-k)*fac_st(k)))); 
    }
    hold=complex<long double>(lhold);
    return hold*(-complex<long double>(exp(-alpha*alpha/2))*sqrt(complex<long double>(Nmax/2*(Nmax/2+1))-complex<long double>((-Nmax/2+l+1)*(-Nmax/2+l)))/complex<long double >(2,0));
}
	  

int main(int argc, char *argv[])
{
  long double Delta, eta, gamma, omega, omega0, alpha,tol,en,Dsum,Nsum;
  int Nmax, nmax,size;
  complex<double> hold;
  unsigned int i,j,k,l;//counters
  ostringstream osseva,osseve,ossdos,ossmjz,ossdjz;
  // initializing variables 
  Nmax=20; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;en=0.8;
  
    if(argc>1){//this is used to set the variables for the program from the command line using flags all can be changed or defults used
    for(i=1;i<argc;i=i+2){
      switch(*argv[i]) {
	case 'N':
	  if(isdigit(*argv[i+1]) ){
	    Nmax=atoi(argv[i+1]);
	  }
	break;
	case 'n':
	  if(isdigit(*argv[i+1]) ){
	    nmax=atoi(argv[i+1]);
	  }
	break;
	case 'W':
	  if(isdigit(*argv[i+1]) ){
	    omega=atof(argv[i+1]);
	  }
	break;
	case 'w':
	  if(isdigit(*argv[i+1]) ){
	    omega0=atof(argv[i+1]);
	  }
	break;
	case 'D':
	  if(isdigit(*argv[i+1]) ){
	    Delta=atof(argv[i+1]);
	  }
	break;
	case 'E':
	  if(isdigit(*argv[i+1]) ){
	    eta=atof(argv[i+1]);
	  }
	break;
	case 'G':
	  if(isdigit(*argv[i+1]) ){
	    gamma=atof(argv[i+1]);
	  }
	break;
	case 'e':
	  if(isdigit(*argv[i+1]) ){
	    en=atof(argv[i+1]);
	  }
	break;
	case 'h':
	  cout<<"Help"<<endl;
	  cout<<"This program is used to solve for the lowest states"<<endl;
	  cout<<"of the Hamiltonian W a^t a + w J_z + (g/N^(1/2))(a+a^t)J_x."<<endl;
	  cout<<"	N: Used to set the number of qubits. This must be even."<<endl;
	  cout<<"	n: Used to set the field dimension. This must be even."<<endl;
	  cout<<"	W: Used to set omega."<<endl;
	  cout<<"	w: Used to set omega0."<<endl;
	  cout<<"	D: Used to set Delta."<<endl;
	  cout<<"	E: Used to set eta."<<endl;
	  cout<<"	G: Used to set gamma."<<endl;
	  cout<<"	e: Used to set the percentage of eigenvalues desired must be between 0 and 1."<<endl;
	  cout<<"	h: Displays this help."<<endl;
	  return 0;
	break;
	Default :
	  cout<<"Not an option.";
	  return 0;
	break;
      }
    }  
  }
  
  
  // calculated values
  /*----------------------------------------*/   
  alpha = 2*gamma/(omega*sqrt(Nmax));
  size= int(Nmax/2)*(nmax+1);
  sp_cx_mat H(size,size);
  sp_cx_mat dJz(size,size);
  /*----------------------------------------*/
  // Creating the Hamiltonian
	
  //here we are setting up the matrix for a^dagger a
  //#pragma omp parallel shared(H)
  {
  //#pragma omp parallel for
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;j++){
     H(i*int(nmax+1)+j,i*int(nmax+1)+j)=omega*i-omega*alpha*alpha*(j+1)*(j+1);
   }
  }
  }
  
//   if(H.has_nan()){
//     cout<<"H dag has a NaN"<<endl;
//     return 0;
//   }
  /* following is the setting up pf the matrix for Jz*/
  
//
 // #pragma omp parallel shared(dJz)
  {
 // #pragma omp parallel for
  for(l=0;l<Nmax/2-1;l++){ 
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	  //if((i+(l+1)*int(nmax+1)<size)&&(j+l*int(nmax+1)<size)){
	    hold=complex<double>(Kloop(j,i,l,Nmax,alpha));
	    dJz(i+(l+1)*int(nmax+1),j+l*int(nmax+1))=hold;
	    dJz(j+(l)*int(nmax+1),i+(l+1)*int(nmax+1))=conj(hold);
	  //}
      }      
    }
  }
  }
  
//   if(dJz.has_nan()){
//     cout<<"dJz has a NaN"<<endl;
//     return 0;
//   }

  H=H+Delta*dJz+eta/Nmax*dJz*dJz;
  
  /*---------------------------------------------*/
  // testing if H has any NaN in it
  
//   if(H.has_nan()){
//     cout<<"H has a NaN"<<endl;
//     return 0;
//   }
  
  
  
  /*------------------------------------------------------------*/
  //getting Eigenvalues and Eigenvectors

  cx_vec eigval;
  cx_mat eigvac;
  eigs_gen(eigval, eigvac, H,int(en*size),"sr");//getting en% of the eigenvalues(eigval)  with smallest real part and corisponding eigenvectors(eigvac)
  
  //cout<<H*eigvac.col(1)-eigval(1)*eigvac.col(1)<<endl;//quick test of eigen value
  
  /*------------------------------------------------------------*/
  //Doing post processing and writing to file
  

  osseva<<"results/eigenval_"<<Nmax<<'_'<<nmax<<'_'<<showpoint<<setprecision(1)<<fixed<<omega<<'_'<<omega0<<'_'<<Delta<<'_'<<eta<<'_'<<gamma<<'_'<<showpoint<<setprecision(2)<<fixed<<en<<".dat";
  ofstream fileeva(osseva.str().c_str());
  fileeva << real(eigval);
  fileeva.close();
  
  osseve<<"results/eigenvec_"<<Nmax<<'_'<<nmax<<'_'<<showpoint<<setprecision(1)<<fixed<<omega<<'_'<<omega0<<'_'<<Delta<<'_'<<eta<<'_'<<gamma<<showpoint<<setprecision(2)<<fixed<<'_'<<en<<".dat";
  ofstream fileeve(osseve.str().c_str());
  for(i=0;i<eigvac.n_rows;i++){
    for(j=0;j<eigvac.n_cols;j++){
      fileeve << real(conj(eigvac(i,j))*eigvac(i,j))<< " ";
    }
    fileeve << "\n";
  }
  fileeve.close();
  
  ossmjz<<"results/mjz_"<<Nmax<<'_'<<nmax<<'_'<<showpoint<<setprecision(1)<<fixed<<omega<<'_'<<omega0<<'_'<<Delta<<'_'<<eta<<'_'<<gamma<<'_'<<showpoint<<setprecision(2)<<fixed<<en<<".dat";
  ofstream filemjz(ossmjz.str().c_str());
  for(i=0;i<int(en*size);i++){
   filemjz << real(eigvac.col(i).t()*dJz*eigvac.col(i))*2/Nmax; 
  }
  filemjz.close(); 
  
  ossdos<<"results/DoS_"<<Nmax<<'_'<<nmax<<'_'<<showpoint<<setprecision(1)<<fixed<<omega<<'_'<<omega0<<'_'<<Delta<<'_'<<eta<<'_'<<gamma<<'_'<<showpoint<<setprecision(2)<<fixed<<en<<".dat";
  ofstream filedos(ossdos.str().c_str());
  for(i=0;i<int(en*size)-21;i++){
    Dsum=0;
    Nsum=0;
    for(j=i;j<i+20;j++){
      Dsum+=real(eigval(j+1)+eigval(j));
      Nsum+=real(eigval(j+1)-eigval(j));
    }
    filedos << Dsum/42 << " " << 21/Nsum;
    filedos << "\n";
  }
  filedos.close();
  
  ossdjz<<"results/DJz_"<<Nmax<<'_'<<nmax<<'_'<<showpoint<<setprecision(1)<<fixed<<omega<<'_'<<omega0<<'_'<<Delta<<'_'<<eta<<'_'<<gamma<<showpoint<<setprecision(2)<<fixed<<'_'<<en<<".dat";
  
  for(i=0;i<eigvac.n_cols;i++){
    eigvac.col(i)=(dJz.t()*eigvac.col(i));
  }
  ofstream fileevem(ossdjz.str().c_str());
  for(i=0;i<eigvac.n_rows;i++){
    for(j=0;j<eigvac.n_cols;j++){
      fileevem << real(eigvac(i,j))<< " "<<imag(eigvac(i,j))<<"j ";
    }
    fileevem << "\n";
  }
  fileevem.close();
  
}