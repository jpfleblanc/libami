#include "ami_calc.hpp"



void NewAmiCalc::evaluate_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars){

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_eval_vars.size(),0);
Im_results.resize(ami_eval_vars.size(),0);
for(int i=0; i<ami_eval_vars.size(); i++){	

// std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);

std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list );

// if(std::abs(calc_result.real())>1  || std::abs(calc_result.imag())>1 ){
// std::cout<<"Large value returned at ext var "<<i<<std::endl;
// std::cout<<calc_result<<std::endl;

// std::complex<double> softened=soften_divergence(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list );
	
// }


Re_results[i]=calc_result.real();
Im_results[i]=calc_result.imag();	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	


void NewAmiCalc::evaluate_mirror_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, internal_state &state, external_variable_list &external){

Re_results.clear(); 
Im_results.clear();
Re_results.resize(external.size(),0);
Im_results.resize(external.size(),0);
for(int i=0; i<external.size(); i++){	

// std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);
NewAmiCalc::ami_vars_list vars_list;
	// graph.ami.construct_ami_vars_list(GG_AMI_MATRIX[ord][num][safest_graph[ord][num]].R0_,GG_AMI_MATRIX[ord][num][safest_graph[ord][num]].prefactor_, graph.current_state, extern_list, vars_list);
construct_ami_mirror_vars_list(AMI, state, external[i], vars_list); 

std::vector<std::complex<double>> calc_results;
calc_results.resize(vars_list.size());
std::complex<double> calc_result(0,0);
double norm=vars_list.size();

for(int j=0; j< vars_list.size(); j++){


calc_results[j]=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  vars_list[j], AMI.Unique, AMI.Rref, AMI.Eval_list );

// if(vars_list.size()>1){
std::cout<<"on mirror "<<j<<" of norm="<<norm<< std::endl;
std::cout<<"Returned value "<< calc_results[j]<<std::endl;
// }
// std::cout<<"Var list size is "<<vars_list.size()<<std::endl;
calc_result+=calc_results[j]/norm;//double(vars_list.size());
}



// if(std::abs(calc_result.real())>1  || std::abs(calc_result.imag())>1 ){
// std::cout<<"Large value returned at ext var "<<i<<std::endl;
// std::cout<<calc_result<<std::endl;

// std::complex<double> softened=soften_divergence(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list );
	
// }



Re_results[i]=calc_result.real();
Im_results[i]=calc_result.imag();	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	

// normalization issue if the number of mirror cases is not constant
void NewAmiCalc::construct_ami_mirror_vars_list(solution_set &AMI, internal_state &state, ext_vars &external, NewAmiCalc::ami_vars_list &vars_list){
	
// first just put the one ami_vars in place for the regular case 	
vars_list.push_back(construct_ami_vars(AMI.R0_, AMI.prefactor_, state, external));	
	
// internal_state ranstate=state;

// for(int i=0; i< ranstate.internal_k_list_.size(); i++){

// ranstate.internal_k_list_[i][0]=random_real(0,2.0*M_PI);
// ranstate.internal_k_list_[i][1]=random_real(0,2.0*M_PI);
	
// }


// vars_list.push_back(construct_ami_vars(AMI.R0_, AMI.prefactor_, ranstate, external));

	
// for each Unique_g generate a sign flipped state ?  Then for each Rref make sure that one of the states is relevant? There will always be an odd number of G's for self energy diagrams and so flipping all will always guarantee a sign change 	
for(int i=0; i< AMI.Unique.size(); i++){

internal_state mirror_state;	
bool success=get_mirror_state(AMI.R0_, AMI.Unique[i], state,external, mirror_state);

if(success){
vars_list.push_back(construct_ami_vars(AMI.R0_, AMI.prefactor_, mirror_state, external));
}
	
	
}

if(vars_list.size()==1){
	
internal_state ranstate=state;

for(int i=0; i< ranstate.internal_k_list_.size(); i++){

ranstate.internal_k_list_[i][0]=random_real(0,2.0*M_PI);
ranstate.internal_k_list_[i][1]=random_real(0,2.0*M_PI);
	
}


vars_list.push_back(construct_ami_vars(AMI.R0_, AMI.prefactor_, ranstate, external));

}
	
	
}

bool NewAmiCalc::get_mirror_state(AmiBase::g_prod_t &R0, AmiBase::g_struct &unique_g, internal_state &state, ext_vars &external, internal_state &mirror_state){
	
// std::cout<<"Looking for mirror state"<<std::endl;	
	
	// check if unique_g is actually the external leg 
bool is_ext=false;	
if(unique_g.alpha_.back()!=0 && unique_g.eps_.back()!=0){
is_ext=true;
for(int i=0; i< unique_g.alpha_.size()-1; i++){
	
if(unique_g.alpha_[i]!=0){
is_ext=false;
}	
	
}

for(int i=0; i< unique_g.eps_.size()-1; i++){
	
if(unique_g.eps_[i]!=0){
is_ext=false;
}	
	
}

}	

if(is_ext){return false;}//it is an external line so can't do anything
	

// std::cout<<"Got here"<<std::endl;

AmiBase::energy_t energy=construct_energy(R0, state, external);


// std::cout<<"Energies are "<<std::endl;
// for(int i=0; i< energy.size(); i++){
	// std::cout<<energy[i]<<" ";	
// }
// std::cout<<std::endl;


std::vector<int> indep_index;
std::vector<int> indep_alpha_index;
// collect list of indep alphas(epsilons actually)
for(int i=0; i< R0.size(); i++){

int ac=0;
int alpha_ind=-1;
for(int a=0; a< R0[i].alpha_.size()-1; a++){

if(R0[i].alpha_[a]!=0){
ac++;
alpha_ind=a;
}	
	
}

if(ac==1){
	
indep_alpha_index.push_back(alpha_ind);
// indep_index.push_back(i);	// this is the wrong index need to know what epsilon this is 

for(int j=0; j< R0[i].eps_.size(); j++){
	if(R0[i].eps_[j]!=0){
		
	
		indep_index.push_back(j);
		
	}
}

}
	
}

// std::cout<<"Independent indexes are "<<std::endl;

// for(int i=0; i< indep_index.size(); i++){
// std::cout<<indep_index[i]<<" "<<indep_alpha_index[i]<<std::endl;	
	
	
// }


double target_val=0;

// find an indep_index inside unique_g 
int ind;
int alpha_ind;
bool found=false;



for(int i=0; i< indep_index.size(); i++){

if(unique_g.eps_[indep_index[i]]!=0){
ind=indep_index[i];
alpha_ind=indep_alpha_index[i];
found=true;
break;
}	
	
}
// int tries=0;
// do{
// tries++;
// int r=random_int(0, indep_index.size()-1);

// ind=indep_index[r];
// alpha_ind=indep_alpha_index[r];

// if(unique_g.eps_[ind]!=0){
	// found=true;
// }
// }while(!found && tries<10);

// std::cout<<"Indexes are "<<ind<<" "<<alpha_ind<<std::endl;

// std::cout<<"Assigning target val"<<std::endl;

if(!found){throw std::runtime_error("Never found an independent epsilon");}

for( int i=0; i< energy.size(); i++){
	// target_val-=unique_g.eps_[i]*energy[i];
	if(i!=ind){
		target_val-=(2.0*unique_g.eps_[i]*energy[i]).real();
	}
}

// std::cout<<"External frequency is "<<external.external_freq_[0].real()<<std::endl;
target_val-=2.0*external.external_freq_[0].real()*unique_g.alpha_.back();

target_val=target_val*unique_g.eps_[ind]; // this flips the sign if necessary

target_val-=unique_g.eps_[ind]*(energy[ind]).real();

// then subtract twice the 

// in principle this is now my target.  I want to manipulate epsilon_ind to get this value 
// std::cout<<"Found target of "<< target_val<<std::endl;

// if target is from -1 to 1 then can find a k that flips this 
if(std::abs(target_val)<2.0){

mirror_state=state;

mirror_state.internal_k_list_[alpha_ind][0]=std::acos( -target_val/4.0 );
mirror_state.internal_k_list_[alpha_ind][1]=state.internal_k_list_[alpha_ind][0];


AmiBase::energy_t this_energy=construct_energy(R0, mirror_state, external);

// std::cout<<"Mirrored energies are "<<std::endl;
// for(int i=0; i< this_energy.size(); i++){
	// std::cout<<this_energy[i]<<" ";	
// }
// std::cout<<std::endl;


return true;

}	
	
	return false;
	
}




// TODO does this need to be the full AMI_MATRIX?

std::complex<double> NewAmiCalc::soften_divergence(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list){

// first check if any unique g has  no alphas in it.  if so save it's epsilon.  And we will generate a state with negative of that energy. 
// since we don't have momenta, lets assume that there is a transformation that maps epsilon->-epsilon 

// first check indep variables.  then check dependent ones?
// could flip the signs of the independent energies only as a first attempt. but without knowing the state can't really do this...
int index;
bool found=false;
for(int i=0; i< unique_g.size(); i++){
	
	int size=unique_g[i].alpha_.size();
	
	int zeros=count(unique_g[i].alpha_.begin(), unique_g[i].alpha_.end(), 0);
	
	if(size==zeros){
		
		std::cout<<"Found zero alphas on unique g "<<i<<std::endl;
		found=true;
		index=i;
		break;
	}
	
	
}

// now I know I want to flip the energies of unique_g

std::cout<<"Current energies are "<<std::endl;
for(int i=0 ; i< external.energy_ .size(); i++){
std::cout<<external.energy_[i]<<" ";
}

std::cout<<std::endl;

k_vect_list_t k_vec;
for(int i=0 ; i< external.energy_ .size(); i++){

k_vector_t this_k;

double kx=std::acos(-external.energy_[i].real()/4.0);
double ky=kx;

this_k.push_back(kx);
this_k.push_back(ky);

k_vec.push_back(this_k);

}


std::cout<<"Current restored "<<std::endl;
for(int i=0 ; i< external.energy_ .size(); i++){
	double val=-2.0*(std::cos(k_vec[i][0])+std::cos(k_vec[i][1]));
std::cout<<val<<" ";
}
std::cout<<std::endl;











}	


void NewAmiCalc::evaluate_solutions(std::vector<std::complex<double>> &results, solution_set &AMI, ami_vars_list &ami_eval_vars){

results.clear();
results.reserve(ami_eval_vars.size());
for(int i=0; i<ami_eval_vars.size(); i++){	


std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list );


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
	throw std::runtime_error(std::string("Frequency size does not match alpha:")+ std::to_string(frequency.size())+" "+std::to_string(R0[0].alpha_.size()));
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
// result[j]=eval_epsilon(state.t_list_[j], construct_k(R0[i].alpha_ , k_list) , R0[i].species_, external.MU_, external.H_, state.disp_);
result[j]=eval_epsilon(state.t_list_[j], state.tp_list_[j], construct_k(R0[i].alpha_ , k_list) , R0[i].species_, external.MU_, external.H_, state.disp_);
}else{

// std::cout<<"Assigning energy to result "<<j<<" "<<R0[i].species_<<std::endl;	
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
	std::cerr<<"Warning - Mismatch in size of R0 and epsilon: If this is not a density calculation (or something with counter-terms) then something might be wrong!"<<std::endl;
	density_warning=false;
	
	}
}
	// throw std::runtime_error("Something wrong with epsilon");}
	
return result;
	
}



std::complex<double> NewAmiCalc::eval_epsilon(hopping_t t, hopping_t tp, NewAmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, AmiBase::disp_type disp){
	
	std::complex<double> output(0,0);
	// std::cout<<"In eval_epsilon with dispersion type "<< disp <<std::endl;
	// print_kvector(k);
	
if(disp==AmiBase::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		// std::cout<<"Evaluated tb "<< -2.0*t*cos(k[i]) <<" with momentum "<< k[i]<<" and hopping "<< t <<std::endl;
	}
	
	double term=-4.0*tp;
	for(int i=0; i<k.size(); i++){
		
		term=term*cos(k[i]);
	}
	output+=term;
	
	
}

// Units are of rydbergs:  We used the atomic Rydberg units. Please see the attachment for the details. In this unit, the length scale is the Bohr radius a_0, and the energy scale is the Rydberg energy e^2/2a_0. Equivalently,  you may set \hbar=1, m=1/2 and e^2/2=1. This is why the dispersion becomes \epsilon_k=k^2/2, and the Coulomb replusion =8*pi/q^2. 
if(disp==AmiBase::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	output+=std::pow(k[i],2);//hf_mu;	
		
	}
}

// assuming that spin=0 is up, and spin=1 is down. then spin-1/2, gives -1/2 for up and 1/2 for down.

if(disp==AmiBase::hf){
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

std::complex<double> NewAmiCalc::eval_epsilon(hopping_t t, NewAmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, AmiBase::disp_type disp){
	
	std::complex<double> output(0,0);
	// std::cout<<"In eval_epsilon with dispersion type "<< disp <<std::endl;
	// print_kvector(k);
	
if(disp==AmiBase::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		// std::cout<<"Evaluated tb "<< -2.0*t*cos(k[i]) <<" with momentum "<< k[i]<<" and hopping "<< t <<std::endl;
	}
}

// Units are of rydbergs:  We used the atomic Rydberg units. Please see the attachment for the details. In this unit, the length scale is the Bohr radius a_0, and the energy scale is the Rydberg energy e^2/2a_0. Equivalently,  you may set \hbar=1, m=1/2 and e^2/2=1. This is why the dispersion becomes \epsilon_k=k^2/2, and the Coulomb replusion =8*pi/q^2. 
if(disp==AmiBase::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

// assuming that spin=0 is up, and spin=1 is down. then spin-1/2, gives -1/2 for up and 1/2 for down.

if(disp==AmiBase::hf){
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

std::complex<double> NewAmiCalc::eval_epsilon(hopping_t t, NewAmiCalc::k_vector_t k, std::complex<double> mu , AmiBase::disp_type disp){
	
	std::complex<double> output(0,0);
	// std::cout<<"In eval_epsilon with dispersion type "<< disp <<std::endl;
	// print_kvector(k);
	
if(disp==AmiBase::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		
	}
}
if(disp==AmiBase::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

if(disp==AmiBase::hf){
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


double NewAmiCalc::hf_energy(double kk){
// std::cout<<"HF q arguement is "<<kk<<std::endl;
double E_kk;
double mass=0.5;
double bolv=get_hf_sigma(kk);
double amu=hf_mu;
	
E_kk=(0.5*(kk*kk/mass)+bolv-amu);//*1.91916*1.91916;

// std::cout<<"Returning "<<E_kk<< " for kk, amu "<< kk<<" "<<amu<<std::endl;
	
return E_kk;	
	
}

double NewAmiCalc::get_hf_sigma(double kk){

double amom=kk;

// std::cout<<"In get_hf_sigma with "<<kk<<std::endl;

int ip_ind=kk/hf_kstep;
int ipc=ptoi[ip_ind];

// std::cout<<"ptoi gave "<< ipc<<std::endl;

int Npg=pgrid.size()-1;

if(ipc<0 || ipc> Npg){
std::cout<<"ptoi gave "<< ipc<<" for kk="<<kk<<std::endl;
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

 