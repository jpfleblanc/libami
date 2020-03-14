//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami.hpp"
#include <iomanip>



void AmiCalc::write_S(S_t &s_array){

for (int i=0; i< s_array.size(); i++){
std::cout<<"Writing to file line 1"<<std::endl;
std::ofstream file;
std::stringstream filename;
filename<<"S_"<< i <<".dat";

file.open(filename.str());

for (int j=0; j< s_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<s_array[i][j].size();k++){

file<<s_array[i][j][k]<<" ";
}
file <<"]"<<std::endl;

}



file.close();

}

}

void AmiCalc::write_S_readable(S_t &s_array){

std::ofstream file;
std::stringstream filename;
filename<<"S_cpp.dat";

std::cout<<"Write_S has size "<< s_array.size()<<std::endl;
file.open(filename.str());
for (int i=0; i< s_array.size(); i++){




for (int j=0; j< s_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<s_array[i][j].size();k++){
//std::cout<<i<<" "<< j<<" "<<k<<std::endl;
file<<s_array[i][j][k]<<" ";
}
file <<"] ";

}

file << std::endl;



}

file.close();

}

void AmiCalc::write_R_readable(R_t &r_array){

std::ofstream file;
std::stringstream filename;
std::stringstream alphafile;
filename<<"R_mnta_cpp.dat";
alphafile<<"R_alpha_cpp.dat";	

file.open(filename.str());

std::cout<<"Write_R has size "<<r_array.size()<<std::endl;
std::cout<<"Write_Ri has size "<<r_array[r_array.size()-1].size()<<std::endl;
std::cout<<"Write_Rprod has size "<<r_array[r_array.size()-1][0].size()<<std::endl;

for( int i =0; i< r_array[r_array.size()-1].size(); i++){


for( int j=0; j< r_array[r_array.size()-1][i].size(); j++){
	
file<<"[ ";	
for( int k=0; k< r_array[r_array.size()-1][i][j].eps_.size(); k++){	

file<<r_array[r_array.size()-1][i][j].eps_[k]<<" ";

}
file<<"] ";
}

	
file<<std::endl;

}	
file.close();

// for alphas

file.open(alphafile.str());

std::cout<<"Write_R has size "<<r_array.size()<<std::endl;
std::cout<<"Write_Ri has size "<<r_array[r_array.size()-1].size()<<std::endl;
std::cout<<"Write_Rprod has size "<<r_array[r_array.size()-1][0].size()<<std::endl;

for( int i =0; i< r_array[r_array.size()-1].size(); i++){


for( int j=0; j< r_array[r_array.size()-1][i].size(); j++){
	
file<<"[ ";	
for( int k=0; k< r_array[r_array.size()-1][i][j].alpha_.size(); k++){	

file<<r_array[r_array.size()-1][i][j].alpha_[k]<<" ";

}
file<<"] ";
}

	
file<<std::endl;

}	
file.close();



	
	
}	

void AmiCalc::write_P_readable(P_t &p_array){

std::ofstream file;
std::stringstream filename;
std::stringstream alphafile;
filename<<"P_mnta_cpp.dat";
alphafile<<"P_alpha_cpp.dat";
file.open(filename.str());
std::cout<<"Write_P has size "<< p_array.size()<<std::endl;

for (int i=0; i< p_array.size(); i++){





for (int j=0; j< p_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<p_array[i][j].size();k++){
std::stringstream  pole_string;
file<<"[ ";
for (int r=0; r<p_array[i][j][k].eps_.size(); r++){
pole_string<<p_array[i][j][k].eps_[r]<<" ";

}

file<< pole_string.str();

file<<"] ";
}
file <<"] ";

}

file<<std::endl;



}
file.close();


file.open(alphafile.str());

for (int i=0; i< p_array.size(); i++){





for (int j=0; j< p_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<p_array[i][j].size();k++){
std::stringstream  pole_string;
file<<"[ ";
for (int r=0; r<p_array[i][j][k].alpha_.size(); r++){
pole_string<<p_array[i][j][k].alpha_[r]<<" ";

}

file<< pole_string.str();

file<<"] ";
}
file <<"] ";

}

file<<std::endl;



}
file.close();




}



void AmiCalc::write_P(P_t &p_array){

for (int i=0; i< p_array.size(); i++){

std::ofstream file;
std::stringstream filename;
filename<<"P_"<<i<<".dat";

file.open(filename.str());

for (int j=0; j< p_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<p_array[i][j].size();k++){
std::stringstream  pole_string;
file<<" ( ";
for (int r=0; r<p_array[i][j][k].alpha_.size(); r++){
pole_string<<p_array[i][j][k].alpha_[r]<<" ";

}

file<< pole_string.str();

file<<" ) ";
}
file <<"]"<<std::endl;

}



file.close();

}

}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



void AmiCalc::print_S( int dim, S_t &S_array){

std::cout<<"Printing S"<<std::endl;
 for (int i=0; i< S_array.size();i++){
std::cout<<"S["<<i<<"]=(";
  for (int j=0; j< S_array[i].size(); j++){
    std::cout<<"{";
    for (int k=0; k< S_array[i][j].size(); k++){

    std::cout<<S_array[i][j][k]<<" ";
    }
    std::cout<<"}";
  }
std::cout<<")";
std::cout<<std::endl;
 }
std::cout<<std::endl;



}

void AmiCalc::print_P( int dim, P_t &P_array){

for (int i=0; i< P_array.size(); i++){
std::cout<< "P["<< i<<"]=";

print_Pi(dim, P_array[i]);

}


}

void AmiCalc::print_Pi( int dim, Pi_t &Pi_array){

for (int i=0; i< Pi_array.size(); i++){
std::cout<<"[";
print_pole_array(Pi_array[i]);
std::cout<<"]";
std::cout<<std::endl;
}
std::cout<<std::endl;

}



void AmiCalc::print_R( int dim, R_t &R_array){

for (int i=0; i< R_array.size(); i++){
std::cout<< "R["<< i<<"]=";

print_g_prod_array(R_array[i]);

}

}

void AmiCalc::print_g_prod_array(g_prod_array_t g_array){

for (int i=0; i< g_array.size(); i++){
print_g_prod_info(g_array[i]);
}


}

void AmiCalc::print_final( int dim, R_t &R_array, P_t &P_array, S_t &S_array){

print_S(dim, S_array);
print_P( dim, P_array);
print_R( dim, R_array);


}




void AmiCalc::print_g_prod_info(AmiCalc::g_prod_t g){

std::cout<<"----------Printing G_product---------- " <<std::endl;
for (std::vector<AmiCalc::g_struct>::iterator it= g.begin(); it != g.end(); ++it){

std::cout<<"----------Printing next---------- " <<std::endl;
print_g_struct_info(it[0]);

}


}


void AmiCalc::print_pole_array(AmiCalc::pole_array_t g){

for (int i=0; i< g.size(); i++){
print_pole_struct_info(g[i]);

}


}


void AmiCalc::print_g_struct_info(AmiCalc::g_struct g){

std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}

void AmiCalc::print_pole_struct_info(AmiCalc::pole_struct g){

std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}

void AmiCalc::print_sign_array(AmiCalc::sign_array_t signs){

std::cout<<"Printing signs"<<std::endl;
std::cout<<"[";
 for (int i=0; i< signs.size();i++){
std::cout<<"(";
  for (int j=0; j< signs[i].size(); j++){
    std::cout<<signs[i][j]<<" ";

  }
std::cout<<")";
 }

std::cout<<"]";


}



void AmiCalc::print_signs(AmiCalc::sign_t signs){

std::cout<<"(";
 for (int i=0; i< signs.size();i++){

     std::cout<<signs[i]<<" ";

  }
std::cout<<")";


}



void AmiCalc::print_epsilon_info(AmiCalc::epsilon_t eps){


//for (std::vector<signed char>::iterator it= eps.begin(); it != eps.end(); ++it){
for (std::vector<int>::iterator it= eps.begin(); it != eps.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}


void AmiCalc::print_alpha_info(AmiCalc::alpha_t alpha){

for (std::vector<int>::iterator it= alpha.begin(); it != alpha.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}

void AmiCalc::print_complex_array(std::vector<std::complex<double>> &array){
	
for(int i=0; i< array.size();i++	){
std::cout<<array[i];	


}
	std::cout<<std::endl;
}

void AmiCalc::print_array(std::vector<std::vector<double>> &array){
	
for(int i=0; i<array.size();i++){
	std::cout<<"( ";
	for( int j=0; j< array[i].size(); j++){
	std::cout<<array[i][j]<<" ";
}
std::cout<<") ";
}
std::cout<<std::endl;
}

void AmiCalc::print_int_vector(std::vector<int> &array){
	
for(int i=0; i<array.size();i++){
	std::cout<<"( ";
	
	std::cout<<array[i]<<" ";

std::cout<<") ";
}
std::cout<<std::endl;
}


double AmiCalc::hf_energy(double kk){
double E_kk;
double mass=0.5;
double bolv=get_hf_sigma(kk);
double amu=hf_mu;
	
E_kk=0.5*(kk*kk/mass)+bolv-amu 	;
	
return E_kk;	
	
}

double AmiCalc::get_hf_sigma(double kk){

double amom=kk;

int ip_ind=kk/hf_kstep;
int ipc=ptoi[ip_ind];

int Npg=pgrid.size()-1;

if(ipc<0 || ipc> Npg){
throw std::runtime_error("Something wrong in get_hf_mu function");
}	

int ipc1, ipc2, ipc3;
double pc1,pc2,pc3;
double Ap1, Ap2, Ap3;
double Fp1, Fp2, Fp3;

double output=0.0;

if(ipc==0){
ipc1=0; ipc2=1; ipc3=2;	
}else if(ipc==Npg){
ipc1=Npg-2; ipc2=Npg-1; ipc3=Npg;	
}else {
ipc1=ipc-1; ipc2=ipc; ipc3=ipc+1;	
}

pc1=pgrid[ipc1]; pc2=pgrid[ipc2]; pc3=pgrid[ipc3];

Fp1=sigma_hf[ipc1]; Fp2=sigma_hf[ipc2]; Fp3=sigma_hf[ipc3];

Ap1=(amom-pc2)*(amom-pc3)/((pc1-pc2)*(pc1-pc3)) ;
Ap2=(amom-pc1)*(amom-pc3)/((pc2-pc1)*(pc2-pc3)); 
Ap3=(amom-pc1)*(amom-pc2)/((pc3-pc1)*(pc3-pc2)); 
    
output = Fp1*Ap1+Fp2*Ap2+Fp3*Ap3 ;

return output;
	
}

