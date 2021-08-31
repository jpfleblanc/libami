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


void AmiSpec::resolve_deltas( sp_terms &sp_terms){

for(int i=0; i< sp_terms.size(); i++){
resolve_deltas(sp_terms[i]);
}	
	
}


// TODO: if the prefactor is not +1 or -1, then need to modify the sign for the term to account for the delta rule delta(2*x)=delta(x)/|2|

// TODO: if there is a delta function that is squared, then it is the same as the single delta function but still gives the volume element due to normalization in the larger space. 
void AmiSpec::resolve_deltas(ami_sp_term &sp_term){
	
std::cout<<"Number of deltas is "<<sp_term.dprod_.size()<<std::endl;

if(	sp_term.dprod_.size()==0){ return;}
	
// look at each delta. create a pole for each, and assign an index_ to it that represents which xi will be replaced.  do this for each and make sure each delta gets a unique xi to replace 
// for each pole with respect to its xi values generate the actual pole 
// for each pole replace the xi in every Aterm, g_prod, fermi_pole AND other remaining deltas. once done mark the delta for removal
// for removal define an empty integer vec of size delta.size(). initialize to zero. set to 1 for each resolved delta.  check at the end that all the entries are 1.  and then remove all of the deltas. otherwise throw an error.  OR remove only the entries that are resolved. 

AmiBase::pole_array_t pv;
pv.resize(sp_term.dprod_.size());

int xi_size=sp_term.dprod_[0].eps_.size();
int delta_size=sp_term.dprod_.size();

std::vector<int> used, assigned;
used.resize(xi_size,0);
assigned.resize(delta_size,0);

for(int i=0; i< xi_size; i++){

// search through deltas to find a non-zero and then break 

	for(int m=0; m< sp_term.dprod_.size(); m++){
		if( assigned[m]==0){
		if(sp_term.dprod_[m].eps_[i]!=0 ){
			pv[m].index_=i;
			used[i]=1;
			assigned[m]=1;
			break;
			
		}
		}
		
	}

if( count(used.begin(), used.end(),1)==delta_size){
	break;
}
	
	
}

if(count(used.begin(), used.end(),1)!= delta_size){
	throw std::runtime_error("Didn't resolve deltas correctly - exiting");
	exit(0);
}

// at this point we should have a list of xi we are going to replace. so generate the poles 

 for(int i=0; i< sp_term.dprod_.size(); i++){
	
	pv[i].eps_.resize( sp_term.dprod_[i].eps_.size());
	
	pv[i].alpha_.resize( sp_term.dprod_[i].alpha_.size());
	
	int this_index=pv[i].index_;
	int this_prefactor=sp_term.dprod_[i].eps_[this_index];
	
	
	for( int m=0; m< pv[i].alpha_.size(); m++){
		pv[i].alpha_[m]=sp_term.dprod_[i].alpha_[m]*(-this_prefactor);
	}
	
	for(int m=0; m< pv[i].eps_.size(); m++){
		if(m==this_index){pv[i].eps_[m]=0;}
		else{
			
			pv[i].eps_[m]= sp_term.dprod_[i].eps_[m]*(-this_prefactor);
		}
		
	}
	
	
}

// now at this point we should have all of the poles accounted for 
// std::cout<<std::count(used.begin(), used.end(),1)<<std::endl;

std::cout<<"BEFORE: Printing poles "<<std::endl;

for(int i=0; i< pv.size(); i++){

std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
print_int_vec(pv[i].alpha_);
std::cout<<" | eps-";
print_int_vec(pv[i].eps_);
std::cout<<std::endl;	
	
	
}


// so now replace poles 
for(int i=0; i< pv.size(); i++){

replace_xi(i,pv, sp_term);	
	
	
}

std::cout<<"AFTER: Printing poles "<<std::endl;

for(int i=0; i< pv.size(); i++){

std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
print_int_vec(pv[i].alpha_);
std::cout<<" | eps-";
print_int_vec(pv[i].eps_);
std::cout<<std::endl;	
	
	
}


// we should have a catch in case something goes wrong. 
sp_term.dprod_.clear();



// for(int i=0; i< sp_term.dprod_.size(); i++){
	// pv[i].eps_=sp_term.dprod_[i].eps_;
	// pv[i].alpha_=sp_term.dprod_[i].alpha_;
// }



	
	
	
	
	
}

void AmiSpec::replace_xi(int i, AmiBase::pole_array_t &pv, ami_sp_term &sp_term){


AmiBase::pole_struct this_pole=pv[i];

// update other poles j>i 
for(int j=0; j< pv.size(); j++){
	if(j==i){ continue;}
	update_spec_pole(pv[i], pv[j].alpha_, pv[j].eps_);
	
}

// update the fermi poles in the sp_term

for(int j=0; j< sp_term.ami_term_.p_list.size(); j++){
	
	update_spec_pole(pv[i],sp_term.ami_term_.p_list[j].alpha_,sp_term.ami_term_.p_list[j].eps_);
	
}

// update the g_prod 
for(int j=0; j< sp_term.ami_term_.g_list.size(); j++){
	
	update_spec_pole(pv[i], sp_term.ami_term_.g_list[j].alpha_, sp_term.ami_term_.g_list[j].eps_);
	
}

// update the Aprod
// note that Aprod is defined by the X_t x_ values not eps_ 
for(int j=0; j< sp_term.aprod_.size(); j++){
	
	update_spec_pole(pv[i], sp_term.aprod_[j].alpha_, sp_term.aprod_[j].x_);
	
}


	
	
	
}

void AmiSpec::update_spec_pole(AmiBase::pole_struct &source_pole, AmiBase::alpha_t &target_alpha, AmiBase::epsilon_t &target_eps){
	
// std::cout<<"Updating : starting with: alpha - ";
// print_int_vec(target_alpha);
// std::cout<<" | and eps=";
// print_int_vec(target_eps);
// std::cout<<std::endl;

// std::cout<<"USING pole: alpha - ";
// print_int_vec(source_pole.alpha_);
// std::cout<<" | and eps=";
// print_int_vec(source_pole.eps_);
// std::cout<<std::endl;
	
	
int index=source_pole.index_;	
if(target_eps[index]==0){ return;} // nothing to do 

for(int m=0; m< target_alpha.size(); m++){
	target_alpha[m]+=source_pole.alpha_[m];
}

for(int m=0; m< target_eps.size(); m++){
	if(m==index){ target_eps[m]=0;}
	else{
	target_eps[m]+=source_pole.eps_[m];
	}
}

	
	

// std::cout<<"Updating : resulting in: alpha - ";
// print_int_vec(target_alpha);
// std::cout<<" | and eps=";
// print_int_vec(target_eps);
// std::cout<<std::endl;
		
	
	return;
	
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







