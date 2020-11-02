#include "ami_calc.hpp"



void NewAmiCalc::evaluate_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars){

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_eval_vars.size(),0);
Im_results.resize(ami_eval_vars.size(),0);
for(int i=0; i<ami_eval_vars.size(); i++){	

// std::cout<<"Evaluating with ami_vars "<<std::endl;
// std::cout<<ami_eval_vars[i].BETA_<<" "<<ami_eval_vars[i].prefactor<<std::endl;
// print_complex_array(ami_eval_vars[i].energy_);
// print_complex_array(ami_eval_vars[i].frequency_);

std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);

// TODO: this seems really dangerous...
// if(calc_result.imag()==0 && std::abs(calc_result.real())>1 ){
	// if(flatten_warning){std::cerr<<"WARNING: Calculation returning large values - Flattening "<< calc_result <<std::endl; flatten_warning=false;}
// calc_result=0.0;
// }		

Re_results[i]=calc_result.real();
Im_results[i]=calc_result.imag();	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	


// TODO does this need to be the full AMI_MATRIX?


void NewAmiCalc::evaluate_solutions(std::vector<std::complex<double>> &results, solution_set &AMI, ami_vars_list &ami_eval_vars){

results.clear();
results.reserve(ami_eval_vars.size());
for(int i=0; i<ami_eval_vars.size(); i++){	


std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);


results.push_back(calc_result);	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	



//

void NewAmiCalc::construct_ami_vars_list(AmiBase::g_prod_t &R0, double prefactor, internal_state &state, external_variable_list &external,ami_vars_list &vars_list){
vars_list.clear();
vars_list.reserve(external.size());
for(int i=0; i<external.size(); i++){
// std::cout<<"Making ami vars on external "<< i<<std::endl;
vars_list.push_back(construct_ami_vars(R0, prefactor, state, external[i]));
}	
	
}

//

// TODO: Does this correctly handle susceptibilities?
// this could be part of the ...?
AmiBase::ami_vars NewAmiCalc::construct_ami_vars(AmiBase::g_prod_t &R0, double prefactor, internal_state &state, ext_vars &external){
	
//energy_t energy={-4,1,-1};
// std::cout<<"Beta value is "<<external.BETA_<<std::endl;
// std::cout<<"Frequency value is "<< external.external_freq_[0]<<std::endl;

AmiBase::energy_t energy=construct_energy(R0, state, external);

// the state 'order_' is actually just the internal k-length - or number of independent variables 
AmiBase::frequency_t frequency;
frequency.reserve(state.order_+1);

for(int i=0;i<state.order_;i++){ frequency.push_back(std::complex<double>(0,0));}

// TODO : this doesn't work with multiple external frequencies 
// if(external.external_freq_.size()!=0){
// frequency.push_back(external.external_freq_[0]); // some number of external frequencies
// }

// This should address above todo: should allow multiple external frequencies.
// TODO: need a check somewhere that the frequency length matches the alpha length 
for(int i=0; i< external.external_freq_.size(); i++){
	
frequency.push_back(external.external_freq_[i]);	
	
}


if(frequency.size()!= R0[0].alpha_.size()){
	throw std::runtime_error("Frequency size does not match alpha");
}


AmiBase::ami_vars final_out(energy, frequency);
final_out.BETA_=external.BETA_;
final_out.prefactor=prefactor; //state.prefactor_;
return final_out;

	
}



NewAmiCalc::k_vector_t NewAmiCalc::construct_k(AmiBase::alpha_t alpha, NewAmiCalc::k_vect_list_t &k){
	
NewAmiCalc::k_vector_t kout(k[0].size(),0);	

for(int j=0; j<kout.size(); j++){
for(int i=0;i<k.size(); i++){
	
kout[j]+= alpha[i]*k[i][j];	
	
}
}

return kout;	
	
}




// TODO: This function should be external to the library 
AmiBase::energy_t NewAmiCalc::construct_energy(AmiBase::g_prod_t &R0, internal_state &state, ext_vars &external){

AmiBase::energy_t result;	
NewAmiCalc::k_vect_list_t k_list;

k_list=state.internal_k_list_;
// if(external.external_k_vector_.size()!=0){
	
// print_g_prod_info(R0);	


for(int i=0; i<external.external_k_list_.size(); i++){
k_list.push_back(external.external_k_list_[i]);
}

// }

// std::cout<<"Momentum list is "<<std::endl;
// print_array(k_list);

int count=0;

result.resize(R0[0].eps_.size(),0);
for(int i=0; i< R0.size(); i++){
	for(int j=0; j<R0[i].eps_.size();j++){
		if(R0[i].eps_[j]==1){
			// std::cout<<"On energy item "<<j<<std::endl;
			// std::cout<<"t list entry is "<<state.t_list_[j]<<std::endl;
if(not_molecule){			
result[j]=eval_epsilon(state.t_list_[j], construct_k(R0[i].alpha_ , k_list) , R0[i].species_, external.MU_, external.H_, state.disp_);
}else{
	
result[j]=-global_hii[R0[i].species_];	
}

// std::cout<<"energy "<<count<<" "<< result[j].real()<<" "<<result.size()<<std::endl;
count++; 
		}
	}
}	

// std::cout<<count<<" "<< R0[0].eps_.size();
if(count != R0[0].eps_.size()){
	
	if(density_warning){
	std::cout<<count<<" "<< R0[0].eps_.size();
	std::cerr<<"Warning - Mismatch in size of R0 and epsilon: If this is not a density calculation then something is wrong!"<<std::endl;
	density_warning=false;
	
	}
}
	// throw std::runtime_error("Something wrong with epsilon");}
	
return result;
	
}


std::complex<double> AmiCalc::eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, disp_type disp){
	
	std::complex<double> output(0,0);
	// std::cout<<"In eval_epsilon with dispersion type "<< disp <<std::endl;
	// print_kvector(k);
	
if(disp==AmiCalc::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		// std::cout<<"Evaluated tb "<< -2.0*t*cos(k[i]) <<" with momentum "<< k[i]<<" and hopping "<< t <<std::endl;
	}
}

// Units are of rydbergs:  We used the atomic Rydberg units. Please see the attachment for the details. In this unit, the length scale is the Bohr radius a_0, and the energy scale is the Rydberg energy e^2/2a_0. Equivalently,  you may set \hbar=1, m=1/2 and e^2/2=1. This is why the dispersion becomes \epsilon_k=k^2/2, and the Coulomb replusion =8*pi/q^2. 
if(disp==AmiCalc::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

// assuming that spin=0 is up, and spin=1 is down. then spin-1/2, gives -1/2 for up and 1/2 for down.

if(disp==AmiCalc::hf){
	double q=0.0;
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	q+=std::pow(k[i],2);	
		
	}
	q=std::sqrt(q);

// std::cout<<"Calling hf_energy qith q "<< q<<std::endl;
output=hf_energy(q);
}else{	
	
output+=H*(spin-0.5);	
	
output -= mu;

}

//TODO we don't subtract mu if the hf dispersion is given. it contains its own mu value. 
	// if(std::abs(output)<0.1){
// std::cout<<"Returning epsilon of "<< output<<std::endl;
	// }	
return -output;	
	
	
	
}

std::complex<double> AmiCalc::eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, std::complex<double> mu , disp_type disp){
	
	std::complex<double> output(0,0);
	std::cout<<"In eval_epsilon with dispersion type "<< disp <<std::endl;
	// print_kvector(k);
	
if(disp==AmiCalc::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		
	}
}
if(disp==AmiCalc::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

if(disp==AmiCalc::hf){
	double q=0.0;
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	q+=std::pow(k[i],2);	
		
	}
	q=std::sqrt(q);

output=hf_energy(q);
}else{	
	
output -= mu;

}
	
	return -output;
}

