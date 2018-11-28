//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami.hpp"
#include <iomanip>




std::complex<double> AmiCalc::evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, external_vars &external){

std::cout<<"Evaluating Result for construction: ";

int dim=parms.N_INT_;
SorF_t SorF_result;

for (int i=0; i< dim-1; i++){

SorF_t SF_left, SF_right;
SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=dot(S_array[i], fermi(parms,P_array[i], external));
std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot
SF_right=dot(S_array[i+1], fermi(parms,P_array[i+1], external));



SorF_result=cross(SF_left,SF_right);
std::cout<<"xS["<<i+1<<"].f(P["<<i+1<<"])";


//std::cout<<"After i "<<i<<"steps, K contains "<<std::endl;
//for( int x=0; x< SorF_result[0].size(); x++){
//std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
//}



} 


std::complex<double> final_result;

std::cout<<"*R"<< std::endl;
final_result=star(parms, SorF_result, R_array[dim], external);


return final_result;


}

/// All the special function structures.


AmiCalc::SorF_t AmiCalc::dot(Si_t Si, SorF_t fermi){

SorF_t output;

for( int i=0; i< Si.size(); i++){

std::vector<double> line;

for(int j=0; j< Si[i].size(); j++){

line.push_back(Si[i][j]*fermi[i][j]);

}

output.push_back(line);

}




return output;
}

AmiCalc::SorF_t AmiCalc::cross(SorF_t left, SorF_t right){
// should check that size() of left==1 or else this won't work. Throw an error.
// also, left[0].size()==right.size()

//std::cout<<"Size of left is "<< left.size() <<std::endl;
//std::cout<<"left[0] and right is "<< left[0].size()<< " "<< right.size() <<std::endl;

SorF_t output;


std::vector<double> line;

for( int i=0; i< left[0].size(); i++){
for( int rj=0; rj< right[i].size(); rj++){

line.push_back(left[0][i]*right[i][rj]);

}
}


//std::cout<<"Lengths are "<< line.size()<<std::endl;

output.push_back(line);


return output;


}

AmiCalc::SorF_t AmiCalc::fermi(ami_parms &parms, Pi_t Pi, external_vars external){

SorF_t output;


for (int i=0; i< Pi.size(); i++){

std::vector<double> group;

  for (int j=0; j< Pi[i].size(); j++){

  group.push_back( fermi_pole(parms, Pi[i][j], external));

}

output.push_back(group);


}



return output;

}

// TODO: Remove filewriting, this is debugging info.

std::complex<double> AmiCalc::star(ami_parms &parms, SorF_t K, Ri_t R, external_vars external){

std::cout<<"Size of left is "<< K.size() <<std::endl;
std::cout<<"Size of left array is "<< K[0].size() <<std::endl;
std::cout<<"Size of right is "<< R.size() <<std::endl;

std::complex<double> output=0;
std::complex<double> term;
std::complex<double> gprod;

std::ofstream file;
file.open("outfile.dat");


for( int i=0; i< K[0].size(); i++)
{

gprod=eval_gprod(parms, R[i], external);
term=K[0][i]*gprod;
//std::isnan(std::real(term))
if( true ){output+= term;

//std::cout<<"Term apparently is finite? "<< term << std::endl;
}

//print_g_prod_info(R[i]);




file <<i<<" "<< K[0][i] <<" "<< std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;

}

file.close();



//std::cout<<"Returning value of "<< output <<std::endl;

return output;

}

std::complex<double> AmiCalc::eval_gprod(ami_parms &parms, g_prod_t g_prod, external_vars external){
std::complex<double> output=0;

std::complex<double> denom_prod=1;
double prefactor=1.0;

double E_REG=parms.E_REG_;

for(int i=0; i< g_prod.size(); i++){
std::complex<double> denom=0;

for(int a=0; a< g_prod[i].alpha_.size(); a++){
denom+=double(g_prod[i].alpha_[a])*external.frequency_[a];

//std::cout<<"Alpha's and frequencies " << g_prod[i].alpha_[a] <<" " << external.frequency_[a] << std::endl;

}


for(int a=0; a< g_prod[i].eps_.size(); a++){
denom+=g_prod[i].eps_[a]*external.energy_[a];
}

//denom+=get_energy_from_g(parms, g_prod[i], external);

//if (denom==std::complex<double>((0,0))){denom=E_REG;}//  prefactor=1.0;}//0.0;}


denom_prod=denom_prod*denom;

}

output=1.0/denom_prod*prefactor;




/////////////////// Debugging stuff 
/*
if (std::isinf(std::real(output))){
std::cout<<"GProd gave "<<output<<std::endl;

for(int i=0; i< g_prod.size(); i++){
std::complex<double> sum1=0;
std::complex<double> sum2=0;
std::complex<double> zero=(0,0);

for(int a=0; a< g_prod[i].alpha_.size(); a++){
sum1+=double(g_prod[i].alpha_[a])*external.frequency_[a];
}


//std::cout<<"Alpha's and frequencies " << sum1 << std::endl;

for(int a=0; a< g_prod[i].eps_.size(); a++){
sum2+=g_prod[i].eps_[a]*external.energy_[a];
}
//std::cout<<"epsilons "<< sum2 <<std::endl;

if (sum1==zero && sum2==zero){

std::cout<<"Both sums were zero "<<std::endl;


for(int a=0; a< g_prod[i].alpha_.size(); a++){


std::cout<<"Alpha's and frequencies " << g_prod[i].alpha_[a] <<" " << external.frequency_[a] << std::endl;

}

for(int a=0; a< g_prod[i].eps_.size(); a++){
std::cout<<"eps and frequencies " << g_prod[i].eps_[a] <<" " << external.energy_[a] << std::endl;
}



}

}


}
//
*/








return output;
}


double AmiCalc::fermi_pole(ami_parms &parms, pole_struct pole, external_vars external){

double output;
int eta=0;

double beta=parms.BETA_;
double E_REG=parms.E_REG_;




for (int i=0; i< pole.alpha_.size(); i++){
eta+= pole.alpha_[i];
}

double E= get_energy_from_pole(pole,external);


double sigma= pow(-1.0, double(eta));

//val = 1.0/(sigma*exp(B*z[0]-gamma)+1.0)
//double gamma=0.01;

output=1.0/(sigma*exp(beta*E-E_REG)+1.0);

//std::cout<< "Fermi output is "<< output << " for beta E and sigma eta "<< beta <<" "<< E<<" "<< sigma <<" "<< eta << std::endl;
//print_pole_struct_info(pole);

return output;
}


double AmiCalc::get_energy_from_pole( pole_struct pole, external_vars external){

double output=0;


for (int i=0; i< pole.eps_.size(); i++){

output+= pole.eps_[i]*external.energy_[i];

}




return output;

}


double AmiCalc::get_energy_from_g( g_struct g, external_vars external){

double output=0;

for (int i=0; i< g.eps_.size(); i++){

output+= g.eps_[i]*external.energy_[i];

}



return output;

}


std::complex<double> AmiCalc::evaluate_multi_random(int NDAT, ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, std::mt19937 &rng){

std::ofstream file;
std::stringstream filename;
filename<<"results_NDAT.dat";

file.open(filename.str());


for (int i=0; i< NDAT;i++){

external_vars ext_f=construct_random_example_J(rng);
std::complex<double> result=AmiCalc::evaluate(parms, R_array, P_array, S_array,  ext_f);

file << i<<" "<<result.real()<< " "<< result.imag() << std::endl;

}

file.close();


}












