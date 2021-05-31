#include "ami_base.hpp"


std::complex<double> AmiBase::evaluate(ami_parms &parms, terms ami_terms, ami_vars &external){

std::complex<double> output(0,0);

for(int i=0; i< ami_terms.size(); i++){
	// std::cout<<"On term="<<i<<" ";
	output+=evaluate_term(parms, ami_terms[i], external);
	
}




return output;
	
}

std::complex<double> AmiBase::evaluate_term(ami_parms &parms, term ami_term, ami_vars &external){

std::complex<double> gprod;

gprod=eval_gprod(parms, ami_term.g_list, external);

std::complex<double> fprod;

fprod=eval_fprod(parms, ami_term.p_list, external);

std::complex<double> output(0,0);

output=ami_term.sign*gprod*fprod;

// std::cout<<"Term gave "<< ami_term.sign*fprod<<" "<< gprod<<" "<<output<<std::endl;

return output;

}

std::complex<double> AmiBase::eval_fprod(ami_parms &parms,pole_array_t &p_list, ami_vars &external){
	
std::complex<double> output(1,0);

for(int i=0; i< p_list.size(); i++){
	
	output=output*fermi_pole(parms, p_list[i], external);
	
}
	
return output;	
	
}


