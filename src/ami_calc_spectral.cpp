#include "ami_calc.hpp"

//TODO: when constructing the spectral representation need to reset all the energies so that everything is a simple-pole 


	std::vector<int> NewAmiCalc::toBinary(int n,int length)
{
	  std::vector<int> r;
    while(n!=0) {
			
			if(n%2==0){ 
			
			r.push_back(0);}else{
				
				r.push_back(1);}
			
			n/=2;
			
			}
		
std::vector<int> temp(length-r.size(),0);
// if(r.size()!=length){
std::reverse(r.begin(), r.end());
temp.insert(temp.end(), r.begin(), r.end());

// }	
			
			
    return temp;
}

void NewAmiCalc::get_pp_comb(int length, std::vector< std::vector<int> > &ppv){
ppv.clear();

int max = std::pow(2,length);	

for(int n=0; n<max ;n++){

ppv.push_back(toBinary(n,length));	
	
}
	
	
}

// TODO: Possible place for error due to negative in energy 
std::complex<double> NewAmiCalc::get_xb_from_pole( AmiBase::pole_struct pole, AmiBase::ami_vars external){

std::complex<double> output(0,0);

// std::cout<<"Evaluating energies"<<std::endl;
for (int i=0; i< pole.eps_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.eps_[i])*external.energy_[i];

}

for (int i=0; i< pole.alpha_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.alpha_[i])*external.frequency_[i];

}


// std::cout<<"Output is "<<output<<std::endl;

return output;

}

void NewAmiCalc::remap_external(AmiBase::ami_vars &external,AmiBase::ami_vars &this_external,AmiBase::Pi_t &Unique_poles, std::vector<int> &pp){

// given pp - for each zero in pp, need to select from Unique_poles one of the possible poles. it doesn't matter so long as they are different indices 

AmiBase::pole_array_t pole_list;
std::vector<int> used;

for(int i=0; i< pp.size(); i++){
	
if(pp[i]==0){

bool found=false;

for(int pole=0; pole<Unique_poles[i].size(); pole++){

int this_xb_index=Unique_poles[i][pole].index_;

if(std::find(used.begin(), used.end(), this_xb_index) == used.end()){
	// was not found in used 
	used.push_back(this_xb_index);
	found=true;
	pole_list.push_back(Unique_poles[i][pole]);
	break; // found a pole so exit the 'pole' for loop 
}



}

// if made it here without finding an appropriate pole then throw an error
// This likely needs to be generalized in some way to find a pole combination when this fails. 
if(!found){ throw std::runtime_error("Failed to find a pole combination: in remap_external for spectral representation");}

}	
	
}

// at this point should have a list of poles in pole_list.  now we want to recursively reassign values to this_external based on replacing the values in the pole.index_ (xb) with the appropriate value 

// the minus is already in this external 
this_external=external;

for(int i=0; i< pole_list.size(); i++){
	
	std::complex<double> this_xb_val=get_xb_from_pole(pole_list[i],this_external);
	this_external.energy_[pole_list[i].index_]=this_xb_val;
	
}
	
		
	return;
	
}




std::complex<double> NewAmiCalc::evaluate_simple_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list){
	
std::complex<double> A_prod=eval_spectral_product(external.energy_, xi_list, external.gamma_);

// std::cout<<"Aprod is "<<A_prod<<" with gamma"<<external.gamma_<<std::endl;

// now duplicate the external and set energy to xi and evaluate the normal function 

AmiBase::ami_vars this_external=external;

// std::cout<<"Evaluating with energies"<<std::endl;
for(int i=0; i< this_external.energy_.size(); i++){
 
this_external.energy_[i]=-xi_list[i];	
	// std::cout<<this_external.energy_[i]<<" ";
	
	
}
// std::cout<<std::endl;


// std::cout<<"Evaluating with xi"<<std::endl;
// for(int i=0; i< this_external.energy_.size(); i++){

// this_external.energy_[i]=-xi_list[i];	
	// std::cout<<xi_list[i]<<" ";
	
	
// }
// std::cout<<std::endl;

// evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external,g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list)

std::complex<double> normal=amibase.evaluate(	parms, R_array, P_array, S_array, this_external,unique_g, Rref, Eval_list);
	
	
return A_prod*normal;	
	
	
}


std::complex<double> NewAmiCalc::eval_spectral_product(std::vector<std::complex<double>> &Ei, std::vector<double> &xi, double &gamma){
	
std::complex<double> output(1,0);

// the energies are the negative of the epsilon values 
for(int i=0; i< Ei.size(); i++){
	// std::cout<<"On i="<<i<<std::endl;
std::complex<double> this_A=gamma/(std::pow(xi[i] +Ei[i],2) + std::pow(gamma,2))/M_PI;	
output=output*this_A;
}
	

return output;	
	
}


std::complex<double> NewAmiCalc::evaluate_spectral_v2(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list, AmiBase::Pi_t &Unique_poles){
	

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Called the evaluate spectral function without optimization");
	
}	


AmiBase::ami_vars gprod_external=external;
// AmiBase::ami_vars aprod_external=external;


for(int i=0; i< gprod_external.energy_.size(); i++){
 
gprod_external.energy_[i]=-xi_list[i];	
	
}



std::complex<double> final_result;
int dim=parms.N_INT_;


for(int opt_term=0; opt_term<Rref.size(); opt_term++){

// first look at what is in the opt_term and now iterate over the different combinations of delta functions 
int ndelta=0;

std::vector<int> delta_indices;

// std::cout<<"Unique

for(int i=0; i< Rref[opt_term].size(); i++){
// std::cout<<"For opt_term the number of pole choices is "<<Unique_poles[Rref[opt_term][i].first].size()<<std::endl;
if( Unique_poles[Rref[opt_term][i].first].size()!=0){
ndelta++;
delta_indices.push_back(i);
}	
	
	
}
// these are the delta combinations to be applied to delta_indices 
std::vector<std::vector<int>> pp_v;
get_pp_comb(ndelta,pp_v);

// std::cout<<"ppv size is "<<pp_v.size()<<" for ndelta="<<ndelta<<std::endl;

for(int delta_term=0; delta_term< pp_v.size(); delta_term++){
	// std::cout<<"On delta_term "<<delta_term<<std::endl;
	
	
	std::vector<int> pp;
	pp.resize(Rref[opt_term].size(),1);
	// expand the pp_comb to include the extra Rref entries 
	// lets define 1 to be normal, and zero to be a delta 
	for(int i=0; i< delta_indices.size(); i++){
		
		pp[delta_indices[i]]=pp_v[delta_term][i];
		
	}
	// now pp should be same length as Rref and tell which are deltas and which are not 

	int ndelta = count(pp.begin(), pp.end(), 0);// need this later for prefactor
  // if(ndelta!=1){continue;}	
	int n_xi_int=R_array[0][0].size()-ndelta;

// Next we need to manipulate this_external to address the mapping - which has not been figured out yet!


// std::cout<<"dim="<<dim<< std::endl;
AmiBase::SorF_t SorF_result;

AmiBase::ami_vars remapped_gprod_external; //=external;
// AmiBase::ami_vars remapped_aprod_external;

remap_external( gprod_external, remapped_gprod_external, Unique_poles, pp);

std::vector<double> this_xi_list(xi_list.size());

// push the remapping onto the xi's 
for(int i=0; i< gprod_external.energy_.size(); i++){
 
this_xi_list[i]=-remapped_gprod_external.energy_[i].real();// TODO: is there any situation in which this would NOT be real?	
	
}


// remap_external( aprod_external, remapped_aprod_external, Unique_poles, pp);


// now we have the external value for this pp combination.  so generate everything as done before 
// EXCEPT - exclude the delta G's 
// add a prefactor (-i pi)^(ndelta)


// QUestion why evaluate gprod if some are not allowed to be there...
// This part probably doesn't work since even for no integral would need to return Aprod*gprod

if(dim==0){
	
	throw std::runtime_error("This part not working yet");
// std::complex<double> gprod;

// gprod=amibase.eval_gprod(parms, R_array[0][0], remapped_gprod_external);

// return gprod;

}

if (dim==1){
AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;	
	
SF_left=amibase.dot(S_array[0], amibase.fermi(parms,P_array[0], remapped_gprod_external));	

SorF_result=SF_left;
	
}

for (int i=0; i< dim-1; i++){

AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=amibase.dot(S_array[i], amibase.fermi(parms,P_array[i], remapped_gprod_external));
 // std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot

SF_right=amibase.dot(S_array[i+1], amibase.fermi(parms,P_array[i+1], remapped_gprod_external));

 SorF_result=amibase.cross(SF_left,SF_right);
 



}
 


std::complex<double> this_result=optimized_spectral_star(parms, SorF_result, unique_g, Rref[opt_term],Eval_list[opt_term], remapped_gprod_external,pp);


std::complex<double> imag(0.,1.0);
 
	
std::complex<double> A_prod=eval_spectral_product(external.energy_, this_xi_list, external.gamma_);

final_result=final_result+this_result*A_prod*std::pow(-imag*M_PI, ndelta);


} // end delta_term 

}// end opt_term loop 




return final_result;

	
	
}



std::complex<double> NewAmiCalc::evaluate_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g, AmiBase::Pi_t &Unique_poles, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, internal_state &state, ext_vars &ext_var, double &xi_cut){
	
	

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Called the evaluate spectral function without optimization");
	// return evaluate(parms,R_array, P_array, S_array, external);
}


std::complex<double> final_result;


 // std::cout<<"Evaluating Result in spectral notation : ";

int dim=parms.N_INT_;


for(int opt_term=0; opt_term<Rref.size(); opt_term++){

// std::cout<<"On opt term "<<opt_term<<std::endl;

// first look at what is in the opt_term and now iterate over the different combinations of delta functions 
int ndelta=0;

std::vector<int> delta_indices;

// std::cout<<"Unique

for(int i=0; i< Rref[opt_term].size(); i++){
// std::cout<<"For opt_term the number of pole choices is "<<Unique_poles[Rref[opt_term][i].first].size()<<std::endl;
if( Unique_poles[Rref[opt_term][i].first].size()!=0){
ndelta++;
delta_indices.push_back(i);
}	
	
	
}
// these are the delta combinations to be applied to delta_indices 
std::vector<std::vector<int>> pp_v;
get_pp_comb(ndelta,pp_v);

// std::cout<<"ppv size is "<<pp_v.size()<<" for ndelta="<<ndelta<<std::endl;

for(int delta_term=0; delta_term< pp_v.size(); delta_term++){
	// std::cout<<"On delta_term "<<delta_term<<std::endl;
	
	
	std::vector<int> pp;
	pp.resize(Rref[opt_term].size(),1);
	// expand the pp_comb to include the extra Rref entries 
	// lets define 1 to be normal, and zero to be a delta 
	for(int i=0; i< delta_indices.size(); i++){
		
		pp[delta_indices[i]]=pp_v[delta_term][i];
		
	}
	// now pp should be same length as Rref and tell which are deltas and which are not 

	int ndelta = count(pp.begin(), pp.end(), 0);// need this later for prefactor
  // if(ndelta!=1){continue;}	
	int n_xi_int=R_array[0][0].size()-ndelta;

// Next we need to manipulate this_external to address the mapping - which has not been figured out yet!


// std::cout<<"dim="<<dim<< std::endl;
AmiBase::SorF_t SorF_result;

AmiBase::ami_vars this_external; //=external;

remap_external( external, this_external, Unique_poles, pp);
// now we have the external value for this pp combination.  so generate everything as done before 
// EXCEPT - exclude the delta G's 
// add a prefactor (-i pi)^(ndelta)




if(dim==0){
std::complex<double> gprod;

gprod=amibase.eval_gprod(parms, R_array[0][0], this_external);

return gprod;

}

if (dim==1){
AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;	
	
SF_left=amibase.dot(S_array[0], amibase.fermi(parms,P_array[0], this_external));	

SorF_result=SF_left;
	
}

for (int i=0; i< dim-1; i++){

AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=amibase.dot(S_array[i], amibase.fermi(parms,P_array[i], this_external));
 // std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot
//std::cout<<"Before i="<<i<<std::endl;
SF_right=amibase.dot(S_array[i+1], amibase.fermi(parms,P_array[i+1], this_external));

// std::cout<<i<<std::endl;

 SorF_result=amibase.cross(SF_left,SF_right);
 // std::cout<<"xS["<<i+1<<"].f(P["<<i+1<<"])";
 

// std::cout<<"After i "<<i<<"steps, K contains "<<std::endl;
// for( int x=0; x< SorF_result[0].size(); x++){
// std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
// }



}
 


std::complex<double> this_result=optimized_spectral_star(parms, SorF_result, unique_g, Rref[opt_term],Eval_list[opt_term], this_external,pp);

// std::cout<<"This result returned "<<this_result<<std::endl;

std::complex<double> imag(0.,1.0);
 // internal_state &state, ext_vars &ext_var
std::complex<double> A_prod=eval_spectral_product(R_array[0][0], state, ext_var,this_external);

final_result=final_result+this_result*A_prod*std::pow(-imag*M_PI, ndelta)*std::pow(2.0*xi_cut,n_xi_int);




} // end delta_term 

}// end opt_term loop 




return final_result;


	
	
}


std::complex<double> NewAmiCalc::eval_spectral_product(AmiBase::g_prod_t &R0,  internal_state &state, ext_vars &external, AmiBase::ami_vars &external_xi){
	
if(R0.size()!=external_xi.energy_.size()){ throw std::runtime_error("Mismatch in R0 and xi list size for spectral evaluation");}	
	
std::complex<double> output(1,0);

// recombine internal and external state momenta 
NewAmiCalc::k_vect_list_t k_list;

k_list=state.internal_k_list_;
// if(external.external_k_vector_.size()!=0){
	
// print_g_prod_info(R0);	


for(int i=0; i<external.external_k_list_.size(); i++){
k_list.push_back(external.external_k_list_[i]);
}



for(int i=0; i< R0.size(); i++){
	
// k_vector_t this_k=construct_k(R0[i].alpha_, k_list);	

std::complex<double> this_E=-eval_epsilon(state.t_list_[i], state.tp_list_[i], construct_k(R0[i].alpha_ , k_list) , R0[i].species_, external.MU_, external.H_, state.disp_);	

// std::cout<<"Evaluating A for i="<<i<<" and gamma "<<external.gamma_<<std::endl;

// double gamma=0.1;
// std::complex<double> this_A=external.gamma_/(std::pow(external_xi.energy_[i] - this_E,2) + std::pow(external.gamma_,2))/M_PI;

std::complex<double> this_A=external.gamma_/(std::pow(external_xi.energy_[i] - this_E,2) + std::pow(external.gamma_,2))/M_PI;

// std::cout<<this_A<<" "<<external_xi.energy_[i]<<" "<<this_E<<" "<<external.gamma_<<" "<<std::pow(external_xi.energy_[i] - this_E,2)<<std::endl;


output=output*this_A;
	
	
}
	
return output;	
}

// std::complex<double> NewAmiCalc::eval_A(k_vector_t &k, std::complex<double> omega){

// return 	
	
	
// }


std::complex<double> NewAmiCalc::optimized_spectral_star(AmiBase::ami_parms &parms, AmiBase::SorF_t K, AmiBase::g_prod_t &unique_g, AmiBase::ref_v_t &Rref,AmiBase::ref_v_t &Eval_list, AmiBase::ami_vars external, std::vector<int> &pp){


std::complex<double> output=0;
std::complex<double> term;
// std::complex<double> gprod;


std::vector< std::complex<double>> unique_vals;

for(int i=0; i< unique_g.size(); i++){
// we pretend each G is a g_prod and use the same function as the normal star operation for evaluating products of G's 
AmiBase::g_prod_t this_term;
this_term.push_back(unique_g[i]);

unique_vals.push_back(amibase.eval_gprod(parms,this_term,external)/external.prefactor); // This removes the overall prefactor for each g. we then add that back later . recent change should not matter. this prefactor is +-1 not to be confused with the prefactor from 
		
}

// now I have the numerical values of our unique G's 

bool verbose=0;

// for(int i=0; i< Eval_list.size(); i++){
	
	std::complex<double> ksum(0,0);
	for(int j=0; j< Eval_list.size(); j++){
		
		ksum+=K[0][Eval_list[j].first]*double(Eval_list[j].second);
		
	}
	// in principle, every entry in eval_list has the same Rref terms 
	AmiBase::ref_v_t pair_vec= Rref;  
	std::complex<double> this_gprod(1,0);
	
	// std::cout<<"Term comprised of unique G indexes ";
	for(int j=0; j< pair_vec.size(); j++){
	
	// we exclude the delta function cases from being multiplied here 
	if(pp[j]!=0){
	this_gprod=this_gprod*unique_vals[pair_vec[j].first];
	}
	
	}
	

	
	term=ksum*this_gprod*external.prefactor; // add back the overall prefactor for this term 
	

	output+=term;
	// There is only one term processed at a time! 


	
// }





return output;

}





// We would feed this the UniqueG vector 
void NewAmiCalc::collect_spectral_poles(AmiBase::g_prod_t &gprod, AmiBase::Pi_t &pa){
	// clear the array 
pa.clear();
pa.resize(gprod.size());	

// std::cout<<"Collecting spectral poles"<<std::endl;
	
int size=gprod[0].eps_.size();

// first lets identify which G's have external frequencies - 
std::vector<int> pole_indices;

for(int i=0; i<gprod.size(); i++){
	
	if(gprod[i].alpha_.back()!=0){
		pole_indices.push_back(i);
	}
	
}

// std::cout<<"Pole indices found with size "<<pole_indices.size()<<std::endl;

// once we have a set of G's that need to be treated we can decide what to do

// Challenge - in order to optimally use delta's we need to decide how to associate each pole with an delta(xi) and make sure we can find a pole for each of these G's.  
// Try 1: Lets assume that we can find one set that works for the all delta case, and use that for the intermediates. 

// first assign poles to anything with a single epsilon 
// std::vector<int> used(size,0);

for(int i=0; i< pole_indices.size(); i++){
	
	// if(std::count(gprod[pole_indices[i]].eps_.begin(),gprod[pole_indices[i]].eps_.end(),0) == size-1){
		
		for(int xb=0; xb<size; xb++){
			// if(used[xb]==1){continue;}
			AmiBase::pole_struct this_pole;
			bool result=false;
			
			// if(used[xb]==0){
			result=find_spectral_pole(xb, gprod[pole_indices[i]], this_pole);
			if(result){
				// used[xb]=1;
				// std::cout<<"Found a pole for xb="<<xb<<std::endl;
				pa[pole_indices[i]].push_back(this_pole);
				// continue; // only need one pole for each xb choice so move to next 
			}
			
		}
		
	// }
}
// now assign poles indescriminatedly






	
	
}

// This function looks at g and obtains the pole relative to the epsilon index xb that must be non-zero. If it is zero then returns false.  If true then the pole is returned 
bool NewAmiCalc::find_spectral_pole(int xb, AmiBase::g_struct &g, AmiBase::pole_struct &pole){

//
AmiBase::pole_struct pole_out;

if(g.alpha_.back()==0){ return false;}	
if(g.eps_[xb]==0){ return false;}

	
pole_out.eps_=g.eps_;
pole_out.alpha_=g.alpha_;
pole_out.index_=xb;

// rather than check for xb just flip all the signs then set the xb one to zero 
for(int i=0; i< pole_out.eps_.size(); i++){
	if(i==xb){continue;} // skip flipping sign of xb 
	pole_out.eps_[i]=-pole_out.eps_[i]*pole_out.eps_[xb];
}

for(int i=0; i< pole_out.alpha_.size(); i++){
	pole_out.alpha_[i]=-pole_out.alpha_[i]*pole_out.eps_[xb];
}
pole_out.eps_[xb]=0;

pole=pole_out;

return true; 	
	
}


void NewAmiCalc::evaluate_simple_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars,   std::vector<double> &xi_list){

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_eval_vars.size(),0);
Im_results.resize(ami_eval_vars.size(),0);

// std::cout<<"Evaluating ext parameters with xi_list"<<std::endl;
// std::cout<<"(";
// for(int i=0; i< xi_list.size(); i++){
	// std::cout<<xi_list[i]<<",";
// }
// std::cout<<")"<<std::endl;


for(int i=0; i<ami_eval_vars.size(); i++){

// std::complex<double> NewAmiCalc::evaluate_simple_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list)

std::complex<double> calc_result=evaluate_simple_spectral(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list,  xi_list);
// graph.ami.evaluate_simple_spectral(test_amiparms,ggm[2][0].ss_vec[0].R_, ggm[2][0].ss_vec[0].P_, ggm[2][0].ss_vec[0].S_,  external, unique, rref, eval_list, XI_LIST);

// double norm=std::pow(2.0*xi_cut,xi_list.size()

Re_results[i]=calc_result.real();
Im_results[i]=calc_result.imag();	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	

void NewAmiCalc::evaluate_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars, internal_state &state,  external_variable_list &ext_var_list, std::vector<double> &xi_list, double &xi_cut){

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_eval_vars.size(),0);
Im_results.resize(ami_eval_vars.size(),0);

// std::cout<<"Evaluating ext parameters with xi_list"<<std::endl;
// std::cout<<"(";
// for(int i=0; i< xi_list.size(); i++){
	// std::cout<<xi_list[i]<<",";
// }
// std::cout<<")"<<std::endl;


for(int i=0; i<ami_eval_vars.size(); i++){

ami_eval_vars[i].frequency_.back()=ami_eval_vars[i].frequency_.back().real();	

if(ami_eval_vars[i].energy_.size()!=xi_list.size()){
throw std::runtime_error("Mismatch in sizes");
}	

for(int j=0; j<ami_eval_vars[i].energy_.size(); j++){
ami_eval_vars[i].energy_[j]=xi_list[j];
}

// for(int j=0; j<ami_eval_vars[i].energy_.size(); j++){
	// std::cout<<ami_eval_vars[i].energy_[j];
// }
// std::cout<<std::endl;


// for(int j=0; j<ami_eval_vars[i].frequency_.size(); j++){
	// std::cout<<ami_eval_vars[i].frequency_[j];
// }
// std::cout<<std::endl;

// print_complex_array(ami_eval_vars[i].energy_);
// print_complex_array(ami_eval_vars[i].frequency_);

// std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i]);

// std::complex<double> spec_calc_result=evaluate_spectral(test_amiparms,ggm[2][0].ss_vec[0].R_, ggm[2][0].ss_vec[0].P_, ggm[2][0].ss_vec[0].S_,  external, unique, unique_poles, rref, eval_list, this_state, extern_list[0]);
// std::cout<<"Calling evaluate function"<<std::endl;
std::complex<double> calc_result=evaluate_spectral(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique,AMI.Unique_poles, AMI.Rref, AMI.Eval_list, state, ext_var_list[i], xi_cut);

// std::cout<<"Got result "<<calc_result<<std::endl;

// std::complex<double> calc_result=amibase.evaluate(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list );


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






/* 
// TODO do this for the OPT notation 
void AmiCalc::produce_spectral_solutions(AmiCalc::solution_set &ss, AmiCalc::solution_set_vec_t &ssv){

// clear the vector just in case 
ssv.clear();	

// presume that the ss should have a spectral notation for every G in the R0.
if(ss.INT_XI.size==0){
	ss.INT_XI.resize(ss.R0_.size(),1);
}

// Now, for R[last].  for each term R[last][i] get multiplied by K[i] in the end.  Here we are going to do the non-optimized notation first.  Where we generate multiple solution set vectors.  Later we can used the optimized notation to substantially simplify this - but it means that different terms have different INT_XI's so this is simpler as a starting point. 





	
	
	
	
} */

