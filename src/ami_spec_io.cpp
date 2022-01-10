#include "ami_spec.hpp"


void AmiSpec::print_sp_terms(sp_terms &sp_terms){

for( int i=0; i< 	sp_terms.size(); i++){
std::cout<<"-------------------"<<std::endl;
std::cout<<"----  Term "<< i<<" ----"<<std::endl;	
std::cout<<"-------------------"<<std::endl;

print_sp_term(sp_terms[i]);

std::cout<<std::endl;
	
}
	
	
}

// ami_sp_term(){}	

// A_prod_t aprod_;
// AmiBase::term ami_term_;
// delta_prod_t dprod_;

// bool root=true;
// int delta_count=0;
	
// };

void AmiSpec::print_sp_term(ami_sp_term &term){
	
std::cout<<"delta_count="<< term.delta_count<<" and "<< term.dprod_.size()<<std::endl;

print_a_prod_t( term.aprod_);
std::cout<<std::endl;
print_delta_prod_t(term.dprod_);
std::cout<<std::endl;
amibase.print_term(term.ami_term_);

	
	
}

void AmiSpec::print_a_prod_t(A_prod_t &Aprod){
	std::cout<<std::endl;
std::cout<<"--- A PRODUCT - SIZE="<<Aprod.size()<<"---"<<std::endl;

for(int i=0; i< Aprod.size(); i++){

std::cout<<"Aterm "<<i<<std::endl;	

print_a_struct(Aprod[i]);
// std::cout<<std::endl;	
	
}

	
	
}

void AmiSpec::print_a_struct(A_struct &A){

std::cout<<"alpha - ";
print_int_vec(A.alpha_);	
std::cout<<std::endl;

std::cout<<"X -";
print_int_vec(A.x_);
std::cout<<std::endl;

std::cout<<"X_alpha -";
print_int_vec(A.x_alpha_);
std::cout<<std::endl;
	
	
}


void AmiSpec::print_delta_prod_t(delta_prod_t &delta_prod){
	
std::cout<<"--- Delta PRODUCT - SIZE="<< delta_prod.size()<<"---"<<std::endl;	
for(int i=0; i< delta_prod.size(); i++){

std::cout<<i<<" ";
print_delta_t(delta_prod[i]);

}	
std::cout<<std::endl;
	
	
}

void AmiSpec::print_delta_t(delta_t &delta){


std::cout<<"alpha - ";
print_int_vec(delta.alpha_);	
std::cout<<std::endl;

std::cout<<"eps -";
print_int_vec(delta.eps_);
std::cout<<std::endl;	
		
	
}

void AmiSpec::print_int_vec(std::vector<int> vec){

std::cout<<"(";
for (std::vector<int>::iterator it= vec.begin(); it != vec.end(); ++it){

std::cout<< *it << ' ';

}
std::cout<<")";
//std::cout << '\n';


}

