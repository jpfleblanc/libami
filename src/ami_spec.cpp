#include "ami_spec.hpp"


AmiSpec::AmiSpec(double xc){
	
xi_cutoff=xc;	
	
}

std::complex<double> AmiSpec::get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps){
	if(ei.size()!= eps.size()){
	throw std::runtime_error("In A_t the epsilon_t and energy_t do not match in size - exiting");	
	}
	std::complex<double> output(0,0);
	for(int i=0; i< ei.size(); i++){
		
		output+=ei[i]*(double)eps[i];
		
	}
	
	return output;
}

std::complex<double> AmiSpec::get_X(X_t &Xsym, xi_t &xi){
	if(Xsym.size()!= xi.size()){
	throw std::runtime_error("In A_t the X_t and xi_t do not match in size - exiting");	
	}
	std::complex<double> output(0,0);
	for(int i=0; i< xi.size(); i++){
		
		output+=(double)Xsym[i]*xi[i];
		
	}
	
	return output;
}

std::complex<double> AmiSpec::eval_Aprod(A_prod_t &Ap, xi_t &xi, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){
	
	std::complex<double> output(1,0);
	
	// A(Sigma, X, E)  : Sigma: self-energy, X: is frequency, E: energy from k_vector 
	for(int i=0; i< Ap.size(); i++){
		
		std::complex<double> this_X=get_X( Ap[i].x_, xi);
		
		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
		std::complex<double> this_E=eval_tb(1.,0., this_k, mu);
		
		std::complex<double> this_sigma=get_sigma(this_k, this_X);
		
		
		output=output*A_eval(this_sigma, this_X, this_E);
		
		
		
	}
	
	
	return output;
	
}

//TODO: connor please do this part :)
std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X){
	
	
}



// todo: probably don't need to pass A if this function takes in X and E already 
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;	
	
}

void AmiSpec::generate_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0){
	
// first lets identify which G's have external frequencies - 
AmiBase::g_prod_t innert;
AmiBase::g_prod_t nert;

// std::vector<int> pole_indices;
// ami_term.g_list
for(int i=0; i<start_term.g_list.size(); i++){
	
	if(start_term.g_list[i].alpha_.back()!=0){
		nert.push_back(start_term.g_list[i]);
	}else{
		innert.push_back(start_term.g_list[i]);
	}
	
}	

// Next we want to create the product of ( a +d)*(b+d)*(c+d)... So we want to specify how many deltas are in each one
int max_deltas=nert.size();

std::vector<std::vector<int>> pp_v;
ami.get_pp_comb(max_deltas,pp_v);	

sp_terms new_terms;
new_terms.resize(pp_v.size()); // one new term for each combination in pp_v 

for(int i=0; i< pp_v.size(); i++){
	std::cout<<"Term:";
	for(int m=0; m< pp_v[i].size(); m++){
	std::cout<<pp_v[i][m]<<" ";
	
	}
	std::cout<<std::endl;
}

for(int i=0; i< pp_v.size(); i++){
	for(int m=0; m< pp_v[i].size(); m++){
	
	if(pp_v[i][m]==0){
				
		new_terms[i].dprod_.push_back(nert[m]);
		new_terms[i].delta_count++;
	}
	
	if(pp_v[i][m]==1){
	
	new_terms[i].ami_term_.g_list.push_back(nert[m]);
	}
	
	}
	
}

// Generate the A product from R0

A_prod_t this_Ap;
R0_to_Aprod(R0, this_Ap);

// std::cout<<"This ap has size "<<this_Ap<<std::endl;

// now attach the innert G's to each new term 

for(int i=0; i< new_terms.size(); i++){
	
new_terms[i].ami_term_.g_list.insert(new_terms[i].ami_term_.g_list.end(), innert.begin(), innert.end());
new_terms[i].ami_term_.p_list= start_term.p_list;
new_terms[i].aprod_=this_Ap;
		
}

// no update to the sign needed - the delta_count will tell about prefactors when evaluating
// each ami_term.p_list is unchanged from the original 
// What is missing is how to handle the A_prod_t - which at least at the start is defined by R0. and the same for every term. 



new_sp_terms=new_terms;	
	
	
	
}



// epsilon_i either by index or by value 
// int eps_index=-1;
// std::complex<double> eps_val;

// AmiBase::epsilon_t eps_;
// X_t x_;
// AmiBase::alpha_t alpha_;
// AmiBase::species_t species_;


// };
void AmiSpec::R0_to_Aprod(AmiBase::g_prod_t &R0, A_prod_t &Ap){
	
// R0 is g_prod_t 
// each g_prod_t has an epsilon and an alpha 	

Ap.clear();
Ap.resize(R0.size());

for(int i=0; i< R0.size(); i++){
	
X_t this_X;
this_X.resize(R0.size(),0);
this_X[i]=1;	

Ap[i].alpha_=R0[i].alpha_;
Ap[i].species_=R0[i].species_;
Ap[i].x_=this_X;
Ap[i].eps_index=i;
// Ap[i].eps_=R0[i].eps_; // not sure this is needed 

}	
	
	
}

std::complex<double> AmiSpec::eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu){
	
	std::complex<double> output(0,0);
	
	for(int i=0; i<k.size();i++){
		
	output+=-2.0*t*cos(k[i]);	
		
	}
	
	// uncomment this block if want t' in dispersion 
	// double term=-4.0*tp;
	// for(int i=0; i<k.size(); i++){
		
		// term=term*cos(k[i]);
	// }
	// output+=term;
	
	
	output -= mu;
	
	return output;
	
}


std::complex<double> AmiSpec::construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

std::complex<double> result=0;

NewAmiCalc::k_vector_t this_k=ami.construct_k(alpha, klist);
result=eval_tb(1.,0., this_k, mu);

	
return result;
	
}



// void generate_sp_terms(ami_term &start_term, sp_terms &new_sp_terms);
// void reduce_deltas(ami_sp_term &term);







