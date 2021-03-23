#include "ami_calc.hpp"



	std::vector<int> NewAmiCalc::toBinary(int n,int length)
{
	  std::vector<int> r;
    while(n!=0) {
			
			if(n%2==0){ 
			
			r.push_back(0);}else{
				
				r.push_back(1);}
			
			n/=2;
			
			}
		
std::vector<int> temp(length-r.size(),0);
// if(r.size()!=length){
std::reverse(r.begin(), r.end());
temp.insert(temp.end(), r.begin(), r.end());

// }	
			
			
    return temp;
}

void NewAmiCalc::get_pp_comb(int length, std::vector< std::vector<int> > &ppv){
	
int max = std::pow(2,length);	
	
	
}


std::complex<double> NewAmiCalc::get_xb_from_pole( AmiBase::pole_struct pole, AmiBase::ami_vars external){

std::complex<double> output(0,0);

// std::cout<<"Evaluating energies"<<std::endl;
for (int i=0; i< pole.eps_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.eps_[i])*external.energy_[i];

}

for (int i=0; i< pole.alpha_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.alpha_[i])*external.frequency_[i];

}


// std::cout<<"Output is "<<output<<std::endl;

return output;

}

void NewAmiCalc::remap_external(AmiBase::ami_vars &external,AmiBase::ami_vars &this_external,AmiBase::Pi_t &Unique_poles, std::vector<int> &pp){

// given pp - for each zero in pp, need to select from Unique_poles one of the possible poles. it doesn't matter so long as they are different indices 

AmiBase::pole_array_t pole_list;
std::vector<int> used;

for(int i=0; i< pp.size(); i++){
	
if(pp[i]==0){

bool found=false;

for(int pole=0; pole<Unique_poles[i].size(); pole++){

int this_xb_index=Unique_poles[i][pole].index_;

if(std::find(used.begin(), used.end(), this_xb_index) == used.end()){
	// was not found in used 
	used.push_back(this_xb_index);
	found=true;
	pole_list.push_back(Unique_poles[i][pole]);
	break; // found a pole so exit the 'pole' for loop 
}



}

// if made it here without finding an appropriate pole then throw an error
// This likely needs to be generalized in some way to find a pole combination when this fails. 
if(!found){ throw std::runtime_error("Failed to find a pole combination: in remap_external for spectral representation");}

}	
	
}

// at this point should have a list of poles in pole_list.  now we want to recursively reassign values to this_external based on replacing the values in the pole.index_ (xb) with the appropriate value 

this_external=external;

for(int i=0; i< pole_list.size(); i++){
	
	std::complex<double> this_xb_val=get_xb_from_pole(pole_list[i],this_external);
	this_external.energy_[pole_list[i].index_]=this_xb_val;
	
}
	
		
	return;
	
}



std::complex<double> NewAmiCalc::evaluate_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g, AmiBase::Pi_t &Unique_poles, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list){
	
	

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	throw std::runtime_error("Called the evaluate spectral function without optimization");
	// return evaluate(parms,R_array, P_array, S_array, external);
}


std::complex<double> final_result;


 // std::cout<<"Evaluating Result for construction: ";

int dim=parms.N_INT_;


for(int opt_term=0; opt_term<Rref.size(); opt_term++){

// first look at what is in the opt_term and now iterate over the different combinations of delta functions 
int ndelta=0;

std::vector<int> delta_indices;

for(int i=0; i< Rref[opt_term].size(); i++){

if( Unique_poles[Rref[opt_term][i].first].size()!=0){
ndelta++;
delta_indices.push_back(i);
}	
	
	
}
// these are the delta combinations to be applied to delta_indices 
std::vector<std::vector<int>> pp_v;
get_pp_comb(ndelta,pp_v);

for(int delta_term=0; delta_term< pp_v.size(); delta_term++){
	
	
	
	std::vector<int> pp;
	pp.resize(Rref[opt_term].size(),1);
	// expand the pp_comb to include the extra Rref entries 
	// lets define 1 to be normal, and zero to be a delta 
	for(int i=0; i< delta_indices.size(); i++){
		
		pp[delta_indices[i]]=pp_v[delta_term][i];
		
	}
	// now pp should be same length as Rref and tell which are deltas and which are not 

	int ndelta = count(pp.begin(), pp.end(), 0);// need this later for prefactor 

// Next we need to manipulate this_external to address the mapping - which has not been figured out yet!


// std::cout<<"dim="<<dim<< std::endl;
AmiBase::SorF_t SorF_result;

AmiBase::ami_vars this_external; //=external;

remap_external( external, this_external, Unique_poles, pp);
// now we have the external value for this pp combination.  so generate everything as done before 
// EXCEPT - exclude the delta G's 
// add a prefactor (-i pi)^(ndelta)




if(dim==0){
std::complex<double> gprod;

gprod=amibase.eval_gprod(parms, R_array[0][0], this_external);

return gprod;

}

if (dim==1){
AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;	
	
SF_left=amibase.dot(S_array[0], amibase.fermi(parms,P_array[0], this_external));	

SorF_result=SF_left;
	
}

for (int i=0; i< dim-1; i++){

AmiBase::SorF_t SF_left, SF_right;
AmiBase::SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=amibase.dot(S_array[i], amibase.fermi(parms,P_array[i], this_external));
 // std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot
//std::cout<<"Before i="<<i<<std::endl;
SF_right=amibase.dot(S_array[i+1], amibase.fermi(parms,P_array[i+1], this_external));

// std::cout<<i<<std::endl;

 SorF_result=amibase.cross(SF_left,SF_right);
 // std::cout<<"xS["<<i+1<<"].f(P["<<i+1<<"])";
 

// std::cout<<"After i "<<i<<"steps, K contains "<<std::endl;
// for( int x=0; x< SorF_result[0].size(); x++){
// std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
// }



}
 



// std::cout<<"*R for dim="<<dim<< std::endl;
// final_result=star(parms, SorF_result, R_array[dim], external);

// final_result=amibase.optimized_star(parms, SorF_result, unique_g, Rref,Eval_list, this_external);
std::complex<double> this_result=optimized_spectral_star(parms, SorF_result, unique_g, Rref,Eval_list, this_external,pp);
std::complex<double> imag(0.,1.0);
final_result=final_result+this_result*std::pow(-imag*M_PI, ndelta);


// if(double_check!=final_result){
	// std::cout<<"Double checking got"<<std::endl;
	// std::cout<<double_check<<" "<<final_result<<std::endl;
	
// }


} // end delta_term 

}// end opt_term loop 




return final_result;


	
	
}


std::complex<double> NewAmiCalc::optimized_spectral_star(AmiBase::ami_parms &parms, AmiBase::SorF_t K, AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, AmiBase::ami_vars external, std::vector<int> &pp){


std::complex<double> output=0;
std::complex<double> term;
// std::complex<double> gprod;


std::vector< std::complex<double>> unique_vals;

for(int i=0; i< unique_g.size(); i++){
// we pretend each G is a g_prod and use the same function as the normal star operation for evaluating products of G's 
AmiBase::g_prod_t this_term;
this_term.push_back(unique_g[i]);

unique_vals.push_back(amibase.eval_gprod(parms,this_term,external)*external.prefactor); // This removes the overall prefactor for each g. we then add that back later 
		
}

// now I have the numerical values of our unique G's 

bool verbose=0;

for(int i=0; i< Eval_list.size(); i++){
	
	std::complex<double> ksum(0,0);
	for(int j=0; j< Eval_list[i].size(); j++){
		
		ksum+=K[0][Eval_list[i][j].first]*double(Eval_list[i][j].second);
		
	}
	// in principle, every entry in eval_list has the same Rref terms 
	AmiBase::ref_v_t pair_vec= Rref[i];  
	std::complex<double> this_gprod(1,0);
	
	// std::cout<<"Term comprised of unique G indexes ";
	for(int j=0; j< pair_vec.size(); j++){
	
	// we exclude the delta function cases from being multiplied here 
	if(pp[j]!=0){
	this_gprod=this_gprod*unique_vals[pair_vec[j].first];
	}
	
	}
	

	
	term=ksum*this_gprod*external.prefactor; // add back the overall prefactor for this term 
	

	output+=term;
	


	
}





return output;

}





// We would feed this the UniqueG vector 
void NewAmiCalc::collect_spectral_poles(AmiBase::g_prod_t &gprod, AmiBase::Pi_t &pa){
	// clear the array 
pa.clear();
pa.resize(gprod.size());	
	
int size=gprod[0].eps_.size();

// first lets identify which G's have external frequencies - 
std::vector<int> pole_indices;

for(int i=0; i<gprod.size(); i++){
	
	if(gprod[i].alpha_.back()!=0){
		pole_indices.push_back(i);
	}
	
}

// once we have a set of G's that need to be treated we can decide what to do

// Challenge - in order to optimally use delta's we need to decide how to associate each pole with an delta(xi) and make sure we can find a pole for each of these G's.  
// Try 1: Lets assume that we can find one set that works for the all delta case, and use that for the intermediates. 

// first assign poles to anything with a single epsilon 
std::vector<int> used(size,0);

for(int i=0; i< pole_indices.size(); i++){
	
	if(std::count(gprod[pole_indices[i]].eps_.begin(),gprod[pole_indices[i]].eps_.end(),0) == size-1){
		
		for(int xb=0; xb<size; xb++){
			AmiBase::pole_struct this_pole;
			bool result=false;
			
			// if(used[xb]==0){
			result=find_spectral_pole(xb, gprod[pole_indices[i]], this_pole);
			if(result){
				// used[xb]=1;
				pa[pole_indices[i]].push_back(this_pole);
				
			}
			
		}
		
	}
	
	
}

	
	
}

// This function looks at g and obtains the pole relative to the epsilon index xb that must be non-zero. If it is zero then returns false.  If true then the pole is returned 
bool NewAmiCalc::find_spectral_pole(int xb, AmiBase::g_struct &g, AmiBase::pole_struct &pole){

//
AmiBase::pole_struct pole_out;

if(g.alpha_.back()==0){ return false;}	
if(g.eps_[xb]==0){ return false;}

	
pole_out.eps_=g.eps_;
pole_out.alpha_=g.alpha_;
pole_out.index_=xb;

// rather than check for xb just flip all the signs then set the xb one to zero 
for(int i=0; i< pole_out.eps_.size(); i++){
	pole_out.eps_[i]=-pole_out.eps_[i];
}
pole_out.eps_[xb]=0;
for(int i=0; i< pole_out.alpha_.size(); i++){
	pole_out.alpha_[i]=-pole_out.alpha_[i];
}

pole=pole_out;

return true; 	
	
}
/* 
// TODO do this for the OPT notation 
void AmiCalc::produce_spectral_solutions(AmiCalc::solution_set &ss, AmiCalc::solution_set_vec_t &ssv){

// clear the vector just in case 
ssv.clear();	

// presume that the ss should have a spectral notation for every G in the R0.
if(ss.INT_XI.size==0){
	ss.INT_XI.resize(ss.R0_.size(),1);
}

// Now, for R[last].  for each term R[last][i] get multiplied by K[i] in the end.  Here we are going to do the non-optimized notation first.  Where we generate multiple solution set vectors.  Later we can used the optimized notation to substantially simplify this - but it means that different terms have different INT_XI's so this is simpler as a starting point. 





	
	
	
	
} */

