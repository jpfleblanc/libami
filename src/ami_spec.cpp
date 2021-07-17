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


// todo: probably don't need to pass A if this function takes in X and E already 
std::complex<double> AmiSpec::A_eval(A_struct &A, std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;	
	
}