#include "ami_base.hpp"

void AmiBase::factorize_terms(terms &ami_terms, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list){
	
Ri_t Ri;

convert_terms_to_ri(ami_terms, Ri);

factorize_Rn(Ri, unique_g, Rref,Eval_list);	
	
	
return; 
	
}


