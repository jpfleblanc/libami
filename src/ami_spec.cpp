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
	
	for(int i=0; i< Ap.size(); i++){
		
		std::complex<double> this_X=get_X( Ap[i].x_, xi);
		
		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
		std::complex<double> this_E=eval_tb(1,0, this_k, mu);
		
		std::complex<double> this_sigma=get_sigma(this_k, this_X);
		
		
		output=output*A_eval(this_sigma, this_X, this_E);
		
		
		
	}
	
	
	return output;
	
}

// todo: probably don't need to pass A if this function takes in X and E already 
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;	
	
}

void generate_sp_terms(ami_term &start_term, sp_terms &new_sp_terms){
	
	
	
	
	
	
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


double AmiSpec::construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

double result=0;


result=eval_tb(1,0, ami.construct_k(alpha, klist), mu);

	
return result;
	
}



// void generate_sp_terms(ami_term &start_term, sp_terms &new_sp_terms);
// void reduce_deltas(ami_sp_term &term);







