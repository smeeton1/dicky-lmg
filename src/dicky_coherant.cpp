
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
  double Delta, eta, gamma, omega, omega0, alpha,tol,en;
  int Nmax, nmax,size;
  complex<double> hold;
  unsigned int i,j,k,l;//counters
  
  // initializing variables 
  Nmax=4; //qubit ensemble dimension must be even
  nmax=2*Nmax; //field dimension only even numbers
  Delta=1.0;eta=0.2;gamma=0.3;omega=1.0;omega0=1.0;en=0.2;
  
    if(argc>1){//this is used to set the variables for the program from the command line using flags all can be changed or defults used
    for(i=1;i<argc;i=i+2){
      switch(*argv[i]) {
	case 'N':
	  if(isdigit(*argv[i+1]) ){
	    Nmax=int(*argv[i+1]-'0');
	  }
	break;
	case 'n':
	  if(isdigit(*argv[i+1]) ){
	    nmax=int(*argv[i+1]-'0');
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
  size= int(Nmax/2)*(nmax+1)+int(nmax/2+1);
  sp_cx_mat H(size,size);//,jp,jp0,dJz;
  sp_cx_mat dJz(size,size);
  /*----------------------------------------*/

	
  //here we are setting up the matrix for a^dagger a
  for(i=0;i<nmax/2+1;i++){H(i,i)=omega*2*i;}
  for(i=0;i<Nmax/2;i++){
   for(j=0;j<nmax+1;j++){
     H(i*int(nmax+1)+j+int(nmax/2)+1,i*int(nmax+1)+j+int(nmax/2)+1)=omega*j-omega*alpha*alpha*(i+1)*(i+1);
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
  
  dJz.set_size(size,size);

  for(l=0;l<Nmax/2+1;l++){ 
    for(i=0;i<nmax+1;i++){
      for(j=0;j<nmax+1;j++){
	if((i+(l+1)*int(nmax+1)+(nmax/2)+1<size)&&(j+l*int(nmax+1)+(nmax/2)+1<size)){
	  for(k=0;k<min(i,j)+1;k++){
	    hold=exp(complex<double>((i+j-2*k)*log(alpha))+complex<double>(j-k)*log(complex<double>(-1))+complex<double>(0.5*(lgamma(i+1)+lgamma(j+1))-(lgamma(i-k+1)+lgamma(j-k+1)+lgamma(k+1))));
	    dJz(i+(l+1)*int(nmax+1)+(nmax/2)+1,j+l*int(nmax+1)+(nmax/2)+1)+=hold;
	    dJz(i+(l)*int(nmax+1)+(nmax/2)+1,j+(l+1)*int(nmax+1)+(nmax/2)+1)+=hold;
	  }
	  hold=complex<double>(exp(-alpha*alpha/2))*sqrt(complex<double>(Nmax/2*(Nmax/2+1))-complex<double>((l+2)*(l+1)));
	  dJz(i+(l+1)*int(nmax+1)+(nmax/2)+1,j+l*int(nmax+1)+(nmax/2)+1)*=-hold/complex< double >(2,0);
	  dJz(i+(l)*int(nmax+1)+(nmax/2)+1,j+(l+1)*int(nmax+1)+(nmax/2)+1)*=-hold/complex< double >(2,0);
	}
      }      
    }
  }

  
  H=H+Delta*dJz+eta/Nmax*dJz*dJz;
  
  cx_vec eigval;
  cx_mat eigvac;
  eigs_gen(eigval, eigvac, H,int(en*size),"sr");
  
  eigval=(2/Nmax)*eigval;
  ofstream fileeva("results/eigenval.dat");
  fileeva << real(eigval);
  fileeva.close();
  
  ofstream fileeve("results/eigenvec.dat");
  for(i=0;i<eigvac.n_rows;i++){
    for(j=0;j<eigvac.n_cols;j++){
      fileeve << real(conj(eigvac(i,j))*eigvac(i,j))<< " ";
    }
    fileeve << "\n";
  }
  fileeve.close();
  
}