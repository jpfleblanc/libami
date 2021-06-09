#include "ami_calc.hpp"


// std::complex<double> NewAmiCalc::evaluate_terms_sp(ami_parms &parms, terms &ami_terms, ami_vars &external,std::vector<double> &xi_list,  AmiBase::Pi_t &Unique_poles, double &xi_cut){

	
	
	
	
	
// }

std::complex<double> NewAmiCalc::evaluate_sp_term(AmiBase::ami_parms &parms, AmiBase::term &ami_term, AmiBase::ami_vars &external, std::vector<double> &xi_list, double &xi_cut){

std::complex<double> output(0,0);

// set energies of the G prod terms 	
AmiBase::ami_vars gprod_external=external;

for(int i=0; i< gprod_external.energy_.size(); i++){
 
gprod_external.energy_[i]=-xi_list[i];	
	
}	

// get unique poles 
AmiBase::Pi_t unique_poles;

int max_deltas=0;		
		collect_spectral_poles(ami_term.g_list, unique_poles);

// std::cout<<"Ami term is N_G="<<ami_term.g_list.size()<<std::endl;

// std::cout<<"Unique poles of this product "<<std::endl;

std::vector<int> delta_indices;
for(int i=0; i< unique_poles.size(); i++){

// std::cout<<i<<" n_poles="<<unique_poles[i].size()<<std::endl;
if(unique_poles[i].size()!=0){
	delta_indices.push_back(i);
	max_deltas++;
}

}	
	


std::vector<std::vector<int>> pp_v;
get_pp_comb(max_deltas,pp_v);	

// std::cout<<"pp_v has size "<<pp_v.size()<<std::endl;	

// for each pp_v entry we assign a new object pp that says	
for(int delta_term=0; delta_term<pp_v.size(); delta_term++){

  std::vector<int> pp;
	pp.resize(ami_term.g_list.size(),1);
	
	for(int i=0; i< delta_indices.size(); i++){
		// this part is tricky, length of pp is size of g_list. but pp_v is just the max possible deltas so could be smaller. 
		pp[delta_indices[i]]=pp_v[delta_term][i];
		
	}
	
	int ndelta = count(pp.begin(), pp.end(), 0);
	int n_xi_int=ami_term.g_list.size()-ndelta;
	
	// std::cout<<"For this term ndelta="<<ndelta<<" and n_xi_int="<<n_xi_int<<std::endl;

// TODO: start here 
	AmiBase::ami_vars remapped_gprod_external; 
	std::vector<int> used;

	remap_external( gprod_external, remapped_gprod_external, unique_poles, pp, used);
	double Nxi=delta_scale(xi_list,xi_cut, used);
	
	std::vector<double> this_xi_list(xi_list.size());
	for(int i=0; i< gprod_external.energy_.size(); i++){
 
	this_xi_list[i]=-remapped_gprod_external.energy_[i].real();
	
	}
	
// Question:  Does 'A' here get the remapped energy or no?
// No: Here the external energies are the epsilon_ki.  while the new xi_list 
// Question: Why is there a minus sign?

// std::cout<<"K state energies "<<std::endl;
// for(int m=0; m< external.energy_.size(); m++){
	
	// std::cout<<external.energy_[m]<<std::endl;
// }

std::complex<double> A_prod=eval_spectral_product(external.energy_, this_xi_list, external.gamma_);	

AmiBase::g_prod_t this_glist;
for(int m=0; m< ami_term.g_list.size(); m++){
	
	if(pp[m]==1){
		
		this_glist.push_back(ami_term.g_list[m]);
		
	}
	
}
std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, this_glist, remapped_gprod_external);

// Do we need to update the poles? just use the remapped energies I think 
// I think the fprod is skipping the xi that gets inserted
// need new function that does fprod but includes BOTH the energy from the 
std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, ami_term.p_list, remapped_gprod_external);

std::complex<double> term_val(0,0);
std::complex<double> norm(0,0);

term_val=ami_term.sign*gprod*fprod;

// std::cout<<"Term "<<delta_term<<" gave "<<term_val<<std::endl;	

std::complex<double> imag(0.,1.0);
// norm=std::pow(-imag*M_PI, ndelta);//
norm=std::pow(-imag*M_PI/(2.0*xi_cut), ndelta);
// i think Nxi function fixes this?
// or Nxi = 1/(2xi_cut)^m which is normalized over that integral 

term_val=term_val*A_prod*norm;//*Nxi;

// std::cout<<"A norm Nxi "<< A_prod<<" "<< norm<<" "<<Nxi<<std::endl;


// std::cout<<"Resulting term Term "<<delta_term<<" gave "<<term_val<<std::endl;	
	
	output+=term_val;

}	
	
	
	return output;
	
}





// TODO:  Start here and go through piece by piece from the start
std::complex<double> NewAmiCalc::evaluate_terms_spectral(AmiBase::ami_parms &parms, AmiBase::terms ami_terms, AmiBase::ami_vars &external, AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref, AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list,  AmiBase::Pi_t &Unique_poles, double &xi_cut){
	

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Called the evaluate spectral function without optimization");
	
}	


AmiBase::ami_vars gprod_external=external;
// AmiBase::ami_vars aprod_external=external;


for(int i=0; i< gprod_external.energy_.size(); i++){
 
gprod_external.energy_[i]=-xi_list[i];	
	
}
 


std::complex<double> final_result(0,0);
int dim=parms.N_INT_;


for(int opt_term=0; opt_term<Rref.size(); opt_term++){
//
// TO try to reduce the variance, pick a random opt_term unique g. and look at its epsilons
int rand_choice=random_int(0,Rref[opt_term].size()-1);

std::vector<double> xi_pair=xi_list;
AmiBase::ami_vars gprod_external_pair=gprod_external;

for(int i=0; i< unique_g[Rref[opt_term][rand_choice].first].eps_.size(); i++){
	
if	(unique_g[Rref[opt_term][rand_choice].first].eps_[i]!=0){
	
	gprod_external_pair.energy_[i]=-gprod_external_pair.energy_[i];
	xi_pair[i]=-xi_pair[i];
	
}
		
}

// first look at what is in the opt_term and now iterate over the different combinations of delta functions 
int ndelta_count=0;

std::vector<int> delta_indices;

// std::cout<<"Unique

for(int i=0; i< Rref[opt_term].size(); i++){
// std::cout<<"For opt_term the number of pole choices is "<<Unique_poles[Rref[opt_term][i].first].size()<<std::endl;
if( Unique_poles[Rref[opt_term][i].first].size()!=0){
ndelta_count++;
delta_indices.push_back(i);
}	
	
	
}
// these are the delta combinations to be applied to delta_indices 
std::vector<std::vector<int>> pp_v;
get_pp_comb(ndelta_count,pp_v);

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
	
	 
	// TODO: is this the best way to figure out the number of xi integrals?
	int n_xi_int=Rref[opt_term].size()-ndelta;

// Next we need to manipulate this_external to address the mapping - which has not been figured out yet!


// std::cout<<"dim="<<dim<< std::endl;
AmiBase::SorF_t SorF_result;
AmiBase::SorF_t SorF_result_pair;

AmiBase::ami_vars remapped_gprod_external; 

std::vector<int> used;
remap_external( gprod_external, remapped_gprod_external, Unique_poles, pp, used);

// double Nxi=delta_scale(xi_list,xi_cut, used);

// std::cout<<"Returned Nxi= "<< Nxi<<std::endl;

std::vector<double> this_xi_list(xi_list.size());

// push the remapping onto the xi's 
for(int i=0; i< gprod_external.energy_.size(); i++){
 this_xi_list[i]=-remapped_gprod_external.energy_[i].real();// 
}



std::complex<double> this_result=1; // need actual function to get the value 

//

std::complex<double> imag(0.,1.0);

std::complex<double> A_prod=eval_spectral_product(external.energy_, this_xi_list, external.gamma_);

final_result+=(this_result*A_prod)*std::pow(-imag*M_PI/(2.0*xi_cut), ndelta);//*Nxi;

} // end delta_term 

}// end 


return final_result;
	
	
}

std::complex<double> NewAmiCalc::evaluate_OPT_spectral_term(AmiBase::ami_parms &parms, AmiBase::terms &ami_terms, AmiBase::ami_vars &external,  AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<int> &pp){
	



if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Called spectral without optimized structure - exiting!");
}	


std::complex<double> output=0;
std::complex<double> term;


std::vector< std::complex<double>> unique_vals;

// std::cout<<"External prefactor is "<<external.prefactor<<std::endl;

for(int i=0; i< unique_g.size(); i++){
// we pretend each G is a g_prod and use the same function as the normal star operation for evaluating products of G's 
AmiBase::g_prod_t this_term;
this_term.push_back(unique_g[i]);

unique_vals.push_back(amibase.eval_gprod(parms,this_term,external)*external.prefactor); // This removes the overall prefactor for each g. we then add that back later 
		
}

// now evaluate each term in list 

for(int i=0; i< Eval_list.size(); i++){
	
	std::complex<double> ksum(0,0);
	for( int j=0; j< Eval_list[i].size(); j++){
		// std::cout<<j<<" "<< eval_fprod(parms, ami_terms[Eval_list[i][j].first].p_list, external)*ami_terms[Eval_list[i][j].first].sign<<" "<<double(Eval_list[i][j].second)<<std::endl;
		std::complex<double> kterm= amibase.eval_fprod(parms, ami_terms[Eval_list[i][j].first].p_list, external);
		
		ksum+= kterm*ami_terms[Eval_list[i][j].first].sign * double(Eval_list[i][j].second);
		
	}
	
	
	AmiBase::ref_v_t pair_vec=Rref[i];
	std::complex<double> this_gprod(1,0);
	for(int j=0; j< pair_vec.size(); j++){
		if(pp[j]!=0){
		this_gprod=this_gprod*unique_vals[pair_vec[j].first];
		}
		}
		
		term=ksum*this_gprod*external.prefactor; 
		
		// std::cout<<"In optimized star K[]*R"<<std::endl;
// std::cout<< std::setprecision(20)<< i<<" "<< ksum <<" "<< std::real(this_gprod)<<" "<<std::imag(this_gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;
		
		
		
		output+= term;
	
	
}


	
return output;	

	
	
	
}





double NewAmiCalc::delta_scale( std::vector<double> xi_list, double xi_cut, std::vector<int> &used){
	
double output=1.0;
double sigma=xi_cut/5.0;
double sig2=2.0*sigma*sigma;
double norm=std::sqrt(2.0*M_PI)*sigma;


for(int i=0; i< used.size(); i++){

	
		
		output*=std::exp(-xi_list[used[i]]*xi_list[used[i]]/sig2)/norm;
	
	
}
	
return output;	
	
}


void NewAmiCalc::randomize_xi(std::vector<double> &xi, int length, double max){

xi.resize(length);

for(int i=0; i<length; i++){
	xi[i]=random_real(-max,max);
}
	
}

std::complex<double> NewAmiCalc::evaluate_terms_simple_spectral(AmiBase::ami_parms &parms, AmiBase::terms ami_terms, AmiBase::ami_vars &external, AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref, AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list){
	

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Need factorized form for spectral evaluation - exiting");
}	

// get A product - easy part right now 
std::complex<double> A_prod=eval_spectral_product(external.energy_, xi_list, external.gamma_);


AmiBase::ami_vars this_external=external;

for(int i=0; i< this_external.energy_.size(); i++){
 
this_external.energy_[i]=-xi_list[i];	
		
}

std::complex<double> normal=amibase.evaluate(parms, ami_terms, this_external, unique_g, Rref, Eval_list);

	
return A_prod*normal;	
	
}
