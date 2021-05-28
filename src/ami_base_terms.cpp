#include "ami_base.hpp"


void AmiBase::construct(int N_INT, g_prod_t R0, terms &terms_out){

terms_out.clear();
	
pole_array_t p;	
term start_term(1.0,p,R0);
terms these_terms;
these_terms.push_back(start_term);

terms next_terms;



for(int index=0; index<N_INT; index++){
// std::cout<<"On step "<<index<<std::endl;
// std::cout<<"Term count is "<<these_terms.size()<<std::endl;
// print_terms(these_terms);

integrate_step(index, these_terms, next_terms);
these_terms=next_terms;	
	
}
	
// put new terms in the list 
terms_out.insert(terms_out.end(), these_terms.begin(), these_terms.end());	
	
	
}

void AmiBase::integrate_step(int index, terms &in_terms, terms &out_terms){
	
	out_terms.clear();
	
for (int t_index=0; t_index< in_terms.size(); t_index++){

pole_array_t poles;
poles=find_poles(index, in_terms[t_index].g_list);

// std::cout<<"On term "<<t_index<<" Found n_poles="<<poles.size()<<std::endl;
// for(int i=0; i<poles.size(); i++){
	
	// std::cout<<"Pole "<<i<<" M="<<poles[i].multiplicity_<<std::endl;
	
// }

// for each pole 

for (int i=0; i < poles.size(); i++){


if(poles[i].multiplicity_==1){

term new_term;
new_term.g_list=simple_residue(in_terms[t_index].g_list, poles[i]);
// take sign from original term and multiply by new one
new_term.sign=in_terms[t_index].sign*get_simple_sign(index, in_terms[t_index].g_list,poles[i]);
// take poles from originating term 
new_term.p_list=in_terms[t_index].p_list;
new_term.p_list.push_back(poles[i]);

out_terms.push_back(new_term);	
	
}else{

terms new_terms;	
terms_general_residue(in_terms[t_index], poles[i], new_terms);

// print_terms(new_terms);
// put new terms in the list 
out_terms.insert(out_terms.end(), new_terms.begin(), new_terms.end());
	
}





}	

// std::cout<<"Out terms now contains "<<out_terms.size()<<std::endl;

}
	
	
	
}

void AmiBase::print_term(term &t){
std::cout<<"--- Term is ---"<<std::endl;
		for(int j=0; j< t.g_list.size(); j++){
		print_g_struct_info(t.g_list[j]);
		}
		std::cout<<"--poles--"<<std::endl;
		for(int j=0; j<t.p_list.size(); j++){
			print_pole_struct_info(t.p_list[j]);
		}
		
		std::cout<<"-------------"<<std::endl;	
	
}

void AmiBase::print_terms(terms &t){
	
	for(int i=0; i<t.size(); i++){
		
		std::cout<<"Term:"<<i<<std::endl;
		print_term(t[i]);
		std::cout<<"-------------"<<std::endl;
	}
	
	
}

void AmiBase::terms_general_residue(term &this_term, pole_struct this_pole, terms &out_terms){

out_terms.clear();

// split this_term into p_list that contains alpha[this_pole.index]!=0  and similarly for g_list. then make W_array from that 

// term innert_part, active_part;
// split_term(this_term, this_pole, innert_part,active_part);

// std::cout<<"This pole is attached to "<<this_pole.which_g_.size()<<" greens functions and has multiplicity "<< this_pole.multiplicity_<<std::endl;

// std::cout<<"Innert part has "<<innert_part.g_list.size()<<std::endl;
// print_term(innert_part);

// std::cout<<"Active part has "<<active_part.g_list.size()<<std::endl;
// print_term(active_part);



double starting_sign;
starting_sign=get_starting_sign(this_term.g_list,this_pole);

term W;
W.g_list=reduce_gprod(this_term.g_list,this_pole);

// std::cout<<"W part has "<<std::endl;
// print_term(W);

// add the pole to the p-list here 
// this is not necessary - the variables external do not enter into derivatives. 
// TODO: test if enabling this changes anything
// W.p_list=active_part.p_list;

W.p_list.push_back(this_pole);

W.sign=starting_sign;
// the W 

terms int_terms;
int_terms.push_back(W);

// std::cout<<"Int terms before derivative "<<std::endl;
// print_terms(int_terms);


for (int m=0; m< this_pole.multiplicity_-1; m++){
	terms temp_terms;
	for(int i=0; i< int_terms.size(); i++){
		// std::cout<<"On derivative loop "<<m<<" "<<i<<std::endl;
	
	take_term_derivative(int_terms[i], this_pole, temp_terms);
	
	
	}
	// std::cout<<"On term "<<std::endl;
	// print_term(int_t
	
	int_terms=temp_terms;
	
}

// now we have all of the terms and their derivatives. so now it is safe to sub in the poles 

// std::cout<<"Int terms contains "<<int_terms.size()<<std::endl;
// print_terms(int_terms);

for( int i=0 ;i< int_terms.size(); i++){
	
	for(int j=0; j< int_terms[i].g_list.size(); j++){
		
		int_terms[i].g_list[j]=update_G_pole(int_terms[i].g_list[j],this_pole);
		
		
		
	}
	
	// for every term have to put the pole list back 
	
	int_terms[i].p_list.insert(int_terms[i].p_list.end(), this_term.p_list.begin(), this_term.p_list.end());
	
	
	//TODO: Test if enabling this changes anything 
	// for(int j=0; j< int_terms[i].p_list.size(); j++){
	
	// int_terms[i].p_list[j]=update_Z_pole(int_terms[i].p_list[j], this_pole);
	
	// }
	
	// for every term need to put back the innert parts 
	
	// int_terms[i].g_list.insert(int_terms[i].g_list.end(), innert_part.g_list.begin(), innert_part.g_list.end());
	// int_terms[i].p_list.insert(int_terms[i].p_list.end(), innert_part.p_list.begin(), innert_part.p_list.end());
	
	
}




// if this gets slow should swap instead
out_terms=int_terms;



}

void AmiBase::take_term_derivative(term &in_term, pole_struct &pole, terms &out_terms){
	// out_terms.clear();
	
terms fd_terms;
terms gd_terms;

// std::cout<<"Taking derivative of term "<<std::endl;
// std::cout<<"Current term is "<<std::endl;

// print_term(in_term);

// We don't have a list. we just have one always. 
// first we take the derivative of the product of p_list f's 

/*  // This was a function to take derivatives of chains of f functions. we just have one always 
for( int one=0; one< in_term.p_list.size(); one++){
	term temp;
	temp.g_list=in_term.g_list;
	temp.p_list=in_term.p_list;
	
	
	pole_struct der=in_term.p_list[one];
	der.der_++;
	
	temp.p_list.push_back(der);
	
	for(int m=0; m< in_term.p_list.size(); m++){
	
	if(m==one){continue;}
	
	temp.p_list.push_back(in_term.p_list[m]);
	
	}
	
	fd_terms.push_back(temp);
	
	
} */

for( int one=0; one< in_term.p_list.size(); one++){
	term temp;
	temp.g_list=in_term.g_list;
	temp.p_list=in_term.p_list;
	
	
	temp.p_list[one].der_++;
		
	
	fd_terms.push_back(temp);
	
	
}




// now do the gprod terms 


for(int i=0; i< in_term.g_list.size(); i++){
	
	int alpha=in_term.g_list[i].alpha_[pole.index_];
	
	term temp_gd_term;
	
	temp_gd_term.p_list=in_term.p_list;
		
	if(alpha!=0){
		
		temp_gd_term.sign=-(double)alpha*in_term.sign;
		for(int m=0; m< in_term.g_list.size(); m++){
			if(i==m){
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
			}else{
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
			}
			
			
		}
		
	gd_terms.push_back(temp_gd_term);	
		
	}
	
	
	
	
}

// std::cout<<"At end of derivative gd and fd terms are sizes "<<gd_terms.size()<<" "<<fd_terms.size()<<std::endl;
// std::cout<<"GD terms "<<std::endl;
// print_terms(gd_terms);
// std::cout<<"FD terms"<< std::endl;
// print_terms(fd_terms);
	
out_terms.insert(out_terms.end(), fd_terms.begin(), fd_terms.end());

out_terms.insert(out_terms.end(), gd_terms.begin(), gd_terms.end());

// std::cout<<"Out terms has size "<<out_terms.size()<<std::endl;	
	
	
}

/* 
void AmiBase::take_gprod_term_derivative(term &in_term, pole_struct &pole, terms &out_terms){
	
term fd_term;
terms gd_terms;

// fd term is easy because there is just one 
fd_term=in_term;
pole_struct pole_der=pole;
pole_der.der_++;
fd_term.p_list.push_back(pole_der);

for(int i=0; i< in_term.g_list.size(); i++){
	
	int alpha=in_term.g_list[i].alpha_[pole.index_];
	
	term temp_gd_term;
	
	temp_gd_term.p_list.push_back(
	
	if(alpha!=0){
		
		temp_gd_term.sign=-(double)alpha*in_term.sign;
		for(int m=0; m< in_term.g_list.size(); m++){
			if(i==m){
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
			}else{
				temp_gd_term.g_list.push_back(in_term.g_list[m]);
			}
			
			
		}
		
	gd_terms.push_back(temp_gd_term);	
		
	}
	
	
	
	
}
	
	
}

 */



void AmiBase::split_term(term &this_term, pole_struct this_pole, term &innert_part, term &active_part){

innert_part.g_list.clear();
innert_part.p_list.clear();
active_part.g_list.clear();
active_part.p_list.clear();

// keep the old sign on the innert part 
innert_part.sign=this_term.sign;
active_part.sign=1.0;

for (int i=0; i< this_term.g_list.size(); i++){
		
	if(this_term.g_list[i].alpha_[this_pole.index_]==0){
		innert_part.g_list.push_back(this_term.g_list[i]);
	}else{
		active_part.g_list.push_back(this_term.g_list[i]);
	}
		
}

for( int i=0; i< this_term.p_list.size(); i++){
	
	if(this_term.p_list[i].alpha_[this_pole.index_]==0){
		innert_part.p_list.push_back(this_term.p_list[i]);
	}else{
		active_part.p_list.push_back(this_term.p_list[i]);
	}
	
	
}


	
	
	
}



