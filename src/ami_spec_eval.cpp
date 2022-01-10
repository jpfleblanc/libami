#include "ami_spec.hpp"

std::complex<double> AmiSpec::evaluate_sp_terms(AmiBase::ami_parms &parms, AmiSpec::sp_terms &sp_terms, AmiSpec::ami_spec_vars &vars){

std::complex<double> sum(0,0);

for(int i=0; i< sp_terms.size(); i++){

// std::cout<<"On term "<<i<<std::endl;
std::complex<double> sp_result=evaluate_sp_term(parms, sp_terms[i], vars);

// std::cout<<"gave "<<sp_result<<std::endl;

sum+=sp_result;
}

return sum;


}

std::complex<double> AmiSpec::evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, AmiSpec::ami_spec_vars &vars){

vars.frequency_.back()=vars.frequency_.back().real();

std::complex<double> output(0,0);

// std::cout<<"Entering eval"<<std::endl;
//TODO: Need some sort of warning or catch because this could be really dangerous....
vars.frequency_.back()=vars.frequency_.back().real();


AmiBase::ami_vars gprod_external(vars.xi_list_, vars.frequency_, vars.BETA_, vars.prefactor);

// AmiBase::ami_vars fprod_external(vars.xi_list_, vars.frequency_, vars.BETA_, vars.prefactor);

// TODO: is this sign swap necessary?
for(int i=0; i< gprod_external.energy_.size(); i++){

gprod_external.energy_[i]=-gprod_external.energy_[i];

}
 std::cout<<"Entering eval A"<<std::endl;

std::complex<double> A_prod=eval_Aprod(sp_term.aprod_, vars.xi_list_, vars.frequency_, vars.k_list_, vars.MU_);
// std::cout<< "A_prod:  "<<A_prod<<std::endl;

// std::cout<<"Entering eval G"<<std::endl;
std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, sp_term.ami_term_.g_list, gprod_external);
// std::cout<< "gprod:  "<<gprod<<std::endl;
	//std::cout<<"Entering eval F"<<std::endl;

std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, sp_term.ami_term_.p_list, gprod_external);
// std::cout<< "fprod:  "<<fprod<<std::endl;
std::complex<double> term_val(0,0);
std::complex<double> norm(0,0);

term_val=sp_term.ami_term_.sign*gprod*fprod;

std::complex<double> imag(0.,1.0);

norm=std::pow(-imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);
// std::cout<< "norm:  "<<norm<<std::endl;

// std::cout<<term_val<<" "<<A_prod<<" "<< norm<<std::endl;

term_val=term_val*A_prod*norm;

output+= term_val;

// std::cout<<"Exiting eval "<<output<<std::endl;

return output;



}
//


//TODO: Remove this...
std::complex<double> AmiSpec::evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, NewAmiCalc::ext_vars &ev,   AmiBase::ami_vars &external, NewAmiCalc::k_vect_list_t &klist,   xi_t &xi_list){

std::complex<double> output(0,0);

// std::cout<<"Entering eval"<<std::endl;

AmiBase::ami_vars gprod_external=external;

// TODO: is this sign swap necessary?
for(int i=0; i< gprod_external.energy_.size(); i++){

gprod_external.energy_[i]=-xi_list[i];

}

// std::cout<<"Entering eval A"<<std::endl;

std::complex<double> A_prod=eval_Aprod(sp_term.aprod_, xi_list, external.frequency_, klist, ev.MU_);

// std::cout<<"Entering eval G"<<std::endl;
std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, sp_term.ami_term_.g_list, gprod_external);

	// std::cout<<"Entering eval F"<<std::endl;
std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, sp_term.ami_term_.p_list, gprod_external);

std::complex<double> term_val(0,0);
std::complex<double> norm(0,0);

term_val=sp_term.ami_term_.sign*gprod*fprod;

std::complex<double> imag(0.,1.0);
norm=std::pow(-imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);

// std::cout<<term_val<<" "<<A_prod<<" "<< norm<<std::endl;

term_val=term_val*A_prod*norm;

output+= term_val;

// std::cout<<"Exiting eval"<<std::endl;

return output;


}





std::complex<double> AmiSpec::eval_Aprod(A_prod_t &Ap, xi_t &xi, AmiBase::frequency_t &freq, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

	std::complex<double> output(1,0);

	// A(Sigma, X, E)  : Sigma: self-energy, X: is frequency, E: energy from k_vector
	for(int i=0; i< Ap.size(); i++){


		std::complex<double> this_X=get_X( Ap[i].x_, xi, Ap[i].x_alpha_, freq);
		// std::cout<<"this_X Aprod:  "<<this_X<<std::endl;
		// for(int freq_count=0;freq_count<freq.size();freq_count++){
		// std::cout<<"freq:  "<<freq[freq_count]<<"  ";
	// }
	// std::cout<<"\n"<<std::endl;
		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
		std::complex<double> this_E=eval_tb(1.,0., this_k, mu);
		// std::cout<<"this X:  "<<this_X<<std::endl;
		// std::cout<<"this_k:  "<<this_k[0]<<","<<this_k[1]<<std::endl;
		// std::cout<<"this_E:  "<<this_E<<std::endl;
		//std::complex<double> this_sigma=get_sigma(this_k, this_X);
		//std::cout<<"got sigma"<<std::endl;
			std::complex<double> this_sigma(0,-0.5);// would replace this with a working sigma


		output=output*A_eval(this_sigma, this_X, this_E);



	}

	// std::cout<<"successsful return of eval_Aprod "<<output<<std::endl;
	return output;

}

// todo: probably don't need to pass A if this function takes in X and E already
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;

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
