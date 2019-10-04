//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami.hpp"
#include <iomanip>



// energy_t energy_;
// frequency_t frequency_;
// double prefactor;
// double BETA_;

void AmiCalc::evaluate_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars){

Re_results.clear(); 
Im_results.clear();
Re_results.reserve(ami_eval_vars.size());
Im_results.reserve(ami_eval_vars.size());
for(int i=0; i<ami_eval_vars.size(); i++){	

// std::cout<<"Evaluating with ami_vars "<<std::endl;
// std::cout<<ami_eval_vars[i].BETA_<<" "<<ami_eval_vars[i].prefactor<<std::endl;
// print_complex_array(ami_eval_vars[i].energy_);
// print_complex_array(ami_eval_vars[i].frequency_);

std::complex<double> calc_result=evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);

// TODO: this seems really dangerous...
// if(calc_result.imag()==0 && std::abs(calc_result.real())>1 ){
	// if(flatten_warning){std::cerr<<"WARNING: Calculation returning large values - Flattening "<< calc_result <<std::endl; flatten_warning=false;}
// calc_result=0.0;
// }		

Re_results.push_back(calc_result.real());
Im_results.push_back(calc_result.imag());	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	


// TODO does this need to be the full AMI_MATRIX?


void AmiCalc::evaluate_solutions(std::vector<std::complex<double>> &results, solution_set &AMI, ami_vars_list &ami_eval_vars){

results.clear();
results.reserve(ami_eval_vars.size());
for(int i=0; i<ami_eval_vars.size(); i++){	

// std::cout<<"Evaluating with ami_vars "<<std::endl;
// std::cout<<ami_eval_vars[i].BETA_<<" "<<ami_eval_vars[i].prefactor<<std::endl;
// print_complex_array(ami_eval_vars[i].energy_);
// print_complex_array(ami_eval_vars[i].frequency_);

std::complex<double> calc_result=evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);

//if(std::abs(calc_result.real())>1){
//	std::cout<<"Calculation returning large values - Flattening "<< calc_result <<std::endl;//}	

results.push_back(calc_result);	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	


std::complex<double> AmiCalc::evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external){


 //std::cout<<"Evaluating Result for construction: ";

int dim=parms.N_INT_;

//std::cout<<"dim="<<dim<< std::endl;
SorF_t SorF_result;

for (int i=0; i< dim-1; i++){

SorF_t SF_left, SF_right;
SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=dot(S_array[i], fermi(parms,P_array[i], external));
//  std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot
//std::cout<<"Before i="<<i<<std::endl;
SF_right=dot(S_array[i+1], fermi(parms,P_array[i+1], external));

// std::cout<<i<<std::endl;

 SorF_result=cross(SF_left,SF_right);
 // std::cout<<"xS["<<i+1<<"].f(P["<<i+1<<"])";
 

// std::cout<<"After i "<<i<<"steps, K contains "<<std::endl;
// for( int x=0; x< SorF_result[0].size(); x++){
// std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
// }



} 


std::complex<double> final_result;

//std::cout<<"*R for dim="<<dim<< std::endl;
final_result=star(parms, SorF_result, R_array[dim], external);


return final_result;


}

/// All the special function structures.


AmiCalc::SorF_t AmiCalc::dot(Si_t Si, SorF_t fermi){

SorF_t output;
output.reserve(Si.size());

for( int i=0; i< Si.size(); i++){

std::vector<std::complex<double> > line;
line.reserve(Si[i].size());
for(int j=0; j< Si[i].size(); j++){

//std::cout<<"Dot term 
line.push_back(Si[i][j]*fermi[i][j]);

}

output.push_back(line);

}


return output;
}

AmiCalc::SorF_t AmiCalc::cross(SorF_t left, SorF_t right){
// should check that size() of left==1 or else this won't work. Throw an error.
// also, left[0].size()==right.size()

// std::cout<<"Size of left is "<< left.size() <<std::endl;
// std::cout<<"left[0] and right is "<< left[0].size()<< " "<< right.size() <<std::endl;

SorF_t output;


std::vector<std::complex<double> > line;
line.reserve(left[0].size()*right.back().size());// just guess that the last one is the biggest for reservation

for( int i=0; i< left[0].size(); i++){
for( int rj=0; rj< right[i].size(); rj++){
//std::cout<<"i and rj are "<<i<<" "<<rj<<std::endl ;
line.push_back(left[0][i]*right[i][rj]);

}
}


// std::cout<<"Lengths are "<< line.size()<<std::endl;

output.push_back(line);


return output;


}

AmiCalc::SorF_t AmiCalc::fermi(ami_parms &parms, Pi_t Pi, ami_vars external){

SorF_t output;
output.reserve(Pi.size());


for (int i=0; i< Pi.size(); i++){

std::vector<std::complex<double> > group;
group.reserve(Pi[i].size());

  for (int j=0; j< Pi[i].size(); j++){

  group.push_back( fermi_pole(parms, Pi[i][j], external));

}

output.push_back(group);


}



return output;

}

// TODO: Remove filewriting, this is debugging info.

std::complex<double> AmiCalc::star(ami_parms &parms, SorF_t K, Ri_t R, ami_vars external){

// std::cout<<"Size of left is "<< K.size() <<std::endl;
// std::cout<<"Size of left array is "<< K[0].size() <<std::endl;
// std::cout<<"Size of right is "<< R.size() <<std::endl;

std::complex<double> output=0;
std::complex<double> term;
std::complex<double> gprod;

// std::ofstream file;
// file.open("outfile.dat",  std::ofstream::out | std::ofstream::app);


for( int i=0; i< K[0].size(); i++)
{

gprod=eval_gprod(parms, R[i], external);
term=K[0][i]*gprod;
//std::isnan(std::real(term))

// if(abs(term)<1){
output+= term;
// }
// if( true ){output+= term;

// std::cout<<"Term apparently is finite? "<< term << std::endl;
// }
// if(std::real(term)>10){
// print_g_prod_info(R[i]);
// }



 // file <<i<<" "<< K[0][i] <<" "<< std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;

}

// file.close();



// std::cout<<"Returning value of "<< output <<std::endl;

return output;

}

std::complex<double> AmiCalc::eval_gprod(ami_parms &parms, g_prod_t g_prod, ami_vars external){
std::complex<double> output(0,0);

std::complex<double> denom_prod(1,0);
double prefactor=external.prefactor;

double E_REG=parms.E_REG_;

// std::cout<<"Eval Gprod"<<std::endl;

for(int i=0; i< g_prod.size(); i++){
std::complex<double> alphadenom(0,0);
std::complex<double> epsdenom(0,0);

for(int a=0; a< g_prod[i].alpha_.size(); a++){
alphadenom+=double(g_prod[i].alpha_[a])*external.frequency_[a];

// std::cout<<"Alpha's and frequencies " << g_prod[i].alpha_[a] <<" " << external.frequency_[a] << std::endl;

}
// TODO: Right here, if the denom==0 still, then the R entry was empty, so regulate the next section, eps -> eps+i0+

std::complex<double> zero(0,0);
std::complex<double> im(0,1);

if(alphadenom==zero){
//alphadenom+=E_REG*im+E_REG;
alphadenom+=E_REG;	
}

// std::cout<<"Energies"<<std::endl;
// Unsure. should this be -=? given that my epsilon is the positive quantity?
for(int a=0; a< g_prod[i].eps_.size(); a++){
epsdenom+=double(g_prod[i].eps_[a])*external.energy_[a];
//denom-=double(g_prod[i].eps_[a])*external.energy_[a];
// std::cout<<g_prod[i].eps_[a]<<" ext "<< external.energy_[a]<<" "<<std::endl;
}
// std::cout<<"End"<<std::endl;

// if(std::abs(epsdenom)< E_REG){
	// std::cout<<"Small epsilon denominator detected"<<std::endl;
// }

//denom+=get_energy_from_g(parms, g_prod[i], external);

//if (denom==std::complex<double>((0,0))){denom=E_REG;}//  prefactor=1.0;}//0.0;}
// std::cout<<"This denom gives term "<<alphadenom+epsdenom<<std::endl;

// if(std::abs(epsdenom)>E_REG){
denom_prod=denom_prod*(alphadenom+epsdenom);
// }


}

// std::cout<<"Denominator product is "<< denom_prod<<std::endl;
// std::cout<<"returning value "<< 1.0/denom_prod*prefactor<<std::endl;

output=1.0/denom_prod*prefactor;









return output;
}


std::complex<double>  AmiCalc::fermi_pole(ami_parms &parms, pole_struct pole, ami_vars external){

std::complex<double>  output;
int eta=0;

double beta=external.BETA_;
double E_REG=parms.E_REG_;




for (int i=0; i< pole.alpha_.size(); i++){
//eta+= pole.alpha_[i];
if(pole.alpha_[i]!=0){
	eta++;
}
}

std::complex<double>  E= get_energy_from_pole(pole,external);


double sigma= pow(-1.0, double(eta));
// if(eta==0){sigma=-1;}

std::complex<double> zero(0,0);
std::complex<double> im(0,1);

// TODO: this might be an actual error. in practice hopefully very small

// oct 4: TODO: comment or not?
if(eta%2==1 && abs(E.real())<abs(E_REG)){//E==zero){
	// std::cout<<"E was "<< E<<std::endl;
// return zero;
E+=E_REG*sgn(E.real());
}	


//if(eta==0){sigma=-1;}

//val = 1.0/(sigma*exp(B*z[0]-gamma)+1.0)
//double gamma=0.01;
//oct 4:  output=1.0/(sigma*std::exp(beta*E-E_REG*im)+1.0);
output=1.0/(sigma*std::exp(beta*(E))+1.0);
// output=1.0/(sigma*exp(beta*E)+1.0);

// if(eta%2==1 && abs(E)<abs(E_REG)){
	// std::cout<<"triggered catch  "<<std::endl;
// return 0.0;
// }	

// std::cout<<"E was  "<< E<<" and output "<<output <<" for eta of "<<eta<<std::endl;

// std::cout<< "Fermi output is "<< output << " for ereg beta E and sigma eta "<<E_REG<<" "<< beta <<" "<< E<<" "<< sigma <<" "<< eta << std::endl;
// print_pole_struct_info(pole);

// if(std::real(output)<0){return zero;}

return output;
}


std::complex<double> AmiCalc::get_energy_from_pole( pole_struct pole, ami_vars external){

std::complex<double> output(0,0);

// std::cout<<"Evaluating energies"<<std::endl;
for (int i=0; i< pole.eps_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.eps_[i])*external.energy_[i];

}


// std::cout<<"Output is "<<output<<std::endl;

return output;

}


std::complex<double>  AmiCalc::get_energy_from_g( g_struct g, ami_vars external){

std::complex<double>  output=0;

for (int i=0; i< g.eps_.size(); i++){

output+= double(g.eps_[i])*external.energy_[i];

}



return output;

}


void AmiCalc::evaluate_multi_random(int NDAT, ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, std::mt19937 &rng){

std::ofstream file;
std::stringstream filename;
filename<<"results_NDAT.dat";

file.open(filename.str());


for (int i=0; i< NDAT;i++){

ami_vars ext_f=construct_random_example_J(rng);
std::complex<double> result=AmiCalc::evaluate(parms, R_array, P_array, S_array,  ext_f);

file << i<<" "<<result.real()<< " "<< result.imag() << std::endl;

}

file.close();


}




AmiCalc::energy_t AmiCalc::construct_energy(AmiCalc::g_prod_t &R0, AmiCalc::internal_state &state, AmiCalc::ext_vars &external){

AmiCalc::energy_t result;	
AmiCalc::k_vect_list_t k_list;

k_list=state.internal_k_list_;
k_list.push_back(external.external_k_vector_);

// std::cout<<"Momentum list is "<<std::endl;
// print_array(k_list);

int count=0;

result.resize(R0[0].eps_.size(),0);
for(int i=0; i< R0.size(); i++){
	for(int j=0; j<R0[i].eps_.size();j++){
		if(R0[i].eps_[j]==1){
result[j]=eval_epsilon( construct_k(R0[i].alpha_ , k_list) , external.MU_);
//std::cout<<count<<" "<< result[j].real()<<std::endl;
count++; 
		}
	}
}	

if(count != R0[0].eps_.size()){
	std::cout<<count<<" "<< R0[0].eps_.size();
	throw std::runtime_error("Something wrong with epsilon");}
	
return result;
	
}

void AmiCalc::construct_ami_vars_list(AmiCalc::g_prod_t &R0, double prefactor, AmiCalc::internal_state &state, AmiCalc::external_variable_list &external,AmiCalc::ami_vars_list &vars_list){
vars_list.clear();
vars_list.reserve(external.size());
for(int i=0; i<external.size(); i++){

vars_list.push_back(construct_ami_vars(R0, prefactor, state, external[i]));
}	
	
}

std::complex<double> AmiCalc::eval_epsilon(AmiCalc::k_vector_t k, std::complex<double> mu ){
	
	std::complex<double> output(0,0);
	// print_kvector(k);
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*cos(k[i]);	
		
	}
	
	output -= mu;
	
	
	return -output;
}

AmiCalc::k_vector_t AmiCalc::construct_k(AmiCalc::alpha_t alpha, AmiCalc::k_vect_list_t &k){
	
AmiCalc::k_vector_t kout(k[0].size(),0);	

for(int j=0; j<kout.size(); j++){
for(int i=0;i<k.size(); i++){
	
kout[j]+= alpha[i]*k[i][j];	
	
}
}

return kout;	
	
}

AmiCalc::ami_vars AmiCalc::construct_ami_vars(AmiCalc::g_prod_t &R0, double prefactor, AmiCalc::internal_state &state, AmiCalc::ext_vars &external){
	
//energy_t energy={-4,1,-1};
// std::cout<<"Beta value is "<<external.BETA_<<std::endl;

AmiCalc::energy_t energy=construct_energy(R0, state, external);

AmiCalc::frequency_t frequency;
frequency.reserve(state.order_+1);

for(int i=0;i<state.order_;i++){ frequency.push_back(std::complex<double>(0,0));}

// TODO : this doesn't work with multiple external frequencies 
frequency.push_back(external.external_freq_[0]); // some number of external frequencies


AmiCalc::ami_vars final_out(energy, frequency);
final_out.BETA_=external.BETA_;
final_out.prefactor=prefactor; //state.prefactor_;
return final_out;

	
}







