#include "ami_base.hpp"

 /**
 * This is the primary AMI symbolic integration function.  It takes a starting integrand defined by `g_prod_t` R0 and `ami_parms` object, and returns the S, P and R arrays necessary for symbolic evaluation. 
 * @param[in] parms : `ami_parms` object, basic parameters for AMI. 
 * @param[in] R0 : `g_prod_t` integrand to be processed.
 * @param[out] R_array: Resultant `R_t` 
 * @param[out] P_array : Resultant `P_t`
 * @param[out] S_array : Resultant `S_t`. 
*/
void AmiBase::construct(ami_parms &parms, g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array){

R_array.clear();
P_array.clear();
S_array.clear();
int dim=parms.N_INT_;

g_prod_array_t R0_array;
R0_array.push_back(R0);
R_array.push_back(R0_array);


pole_array_t poles;


//Explanation: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.


for (int index=0; index<dim; index++){

update_gprod_general(index, index, R_array, P_array,S_array);

}


}



/**
*
* Primary loop.  Takes R[n] as input and outputs R[n+1], P[n+1] and S[n+1].
* @param[in] int_index: indicator for which is the integration variable
* @param[in] array_index: indicator for which element of `alpha_t` corresponds to that integration variable. Typically the same as int_index, unless you want to perform out of order. 
* @param[in] R_array: R[n] input 
* @param[out] R_array: updated with R[n+1] on output. 
* @param[out] S_array: updated with S[n+1] on output.
* @param[out] P_array: updated with P[n+1] on output.
*
*/
void AmiBase::update_gprod_general(int int_index, int array_index, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array){

Ri_t g_in;


g_prod_array_t temp_g_array;
Pi_t temp_pole_array;
Si_t temp_sign_array;


for (int j=0; j<R_array[array_index].size(); j++)
{

sign_t s_out, temp_s;
pole_array_t p_out, temp_p;

pole_array_t poles;
poles=find_poles(int_index, R_array[array_index][j]);  


pole_array_t cor_poles;
sign_t col_signs;


for (int i=0; i < poles.size(); i++){


if(poles[i].multiplicity_==1){ 

g_in.push_back(simple_residue(R_array[array_index][j],poles[i]));

s_out.push_back(get_simple_sign(int_index,R_array[array_index][j], poles[i]));

p_out.push_back(poles[i]);

}else
{

Ri_t temp_g;

evaluate_general_residue(R_array[array_index][j], poles[i], temp_g, temp_p, temp_s);


g_in.insert(std::end(g_in), std::begin(temp_g), std::end(temp_g));
s_out.insert(std::end(s_out),std::begin(temp_s), std::end(temp_s));
p_out.insert(std::end(p_out),std::begin(temp_p), std::end(temp_p));


}


}


temp_pole_array.push_back(p_out);
temp_sign_array.push_back(s_out);

p_out.clear();
s_out.clear();

}


R_array.push_back(g_in);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);

}


/**
* Used in pole multiplicity == 1 - aka, simple pole.  Performs residue theorem on G_in and produces a new `g_prod_t` as output. 
*
*@param[in] G_in: `g_prod_t` to evaluate residues with respect to pole.
*@param[in] pole: `pole_struct` defining the pole. 
*
*@return is of type `g_prod_t`. 
*/
AmiBase::g_prod_t AmiBase::simple_residue(AmiBase::g_prod_t G_in, AmiBase::pole_struct pole){

g_prod_t residue;

if (pole.multiplicity_!=1){std::cerr<<"Simple residue called for pole with M!=1"<<std::endl;}


for (int i=0; i< G_in.size(); i++){


if(i != pole.which_g_[0]){


residue.push_back(update_G_pole(G_in[i],pole) );


}


}

return residue;
}


/**
* Manipulates a Green's function to replace residue variable 'z' with complex pole. 
*
*/
AmiBase::g_struct AmiBase::update_G_pole(AmiBase::g_struct g_in, AmiBase::pole_struct pole){

AmiBase::g_struct g_new;

for (int i=0; i< g_in.alpha_.size(); i++){

if (i != pole.index_){
g_new.alpha_.push_back( g_in.alpha_[i] + g_in.alpha_[pole.index_]*pole.alpha_[i]   );
}
else{
g_new.alpha_.push_back(0);
}

}


for (int i=0; i< g_in.eps_.size(); i++){
g_new.eps_.push_back( g_in.eps_[i] + g_in.alpha_[pole.index_]*pole.eps_[i]   );
}


if ( g_new.alpha_.size() != g_in.alpha_.size()){std::cerr<<"Error: Something may be wrong. Alphas of G and pole not the same size" <<std::endl; }


return g_new;

}



/**
* Manipulates a Pole_struct to replace residue variable 'z' with complex pole. 
*
*/
AmiBase::pole_struct AmiBase::update_Z_pole(AmiBase::pole_struct p_in, AmiBase::pole_struct pole){

AmiBase::pole_struct p_new;

for (int i=0; i< p_in.alpha_.size(); i++){

if (i != pole.index_){
p_new.alpha_.push_back( p_in.alpha_[i] + p_in.alpha_[pole.index_]*pole.alpha_[i]   );
}
else{
p_new.alpha_.push_back(0);
}

}


for (int i=0; i< p_in.eps_.size(); i++){
p_new.eps_.push_back( p_in.eps_[i] + p_in.alpha_[pole.index_]*pole.eps_[i]   );
}


if ( p_new.alpha_.size() != p_in.alpha_.size()){std::cerr<<"Error: Something may be wrong. Alphas of G and pole not the same size" <<std::endl; }


return p_new;

}



/**

* Simple function to scan through `g_prod_t` and collect an std::vector of poles with respect to index.

*/
AmiBase::pole_array_t AmiBase::find_poles(int index, AmiBase::g_prod_t &R){

pole_array_t pole_array;
pole_struct pole;


pole.alpha_.reserve(R[0].alpha_.size());
pole.eps_.reserve(R[0].eps_.size());

pole.index_=index;

if(pole.index_>R[0].alpha_.size()){

std::cerr<<"WARNING: Pole index exceeds length of g_prod_t: This may result in errors.  Exiting find_poles for safety."<<std::endl;	
return pole_array;	
}

for (int i=0; i< R.size(); i++){


  if (R[i].alpha_[index] != 0){

  pole.which_g_.push_back(i);

  for (int j=0; j<R[i].alpha_.size(); j++){
  if (j!=index){pole.alpha_.push_back(-R[i].alpha_[index]*R[i].alpha_[j]); }
  else {pole.alpha_.push_back(0); }
  }

  for (int j=0; j<R[i].eps_.size(); j++){
  pole.eps_.push_back(-R[i].alpha_[index]*R[i].eps_[j]) ;
  }

//Context:  So now you have a pole... need to check if it is already in the array, and then if it isn't add it.  if it is, modify the multiplicity of the pole in the array it is equal to. and add it's which_g_ index to the pole_array with multiplicity

bool duplicate=false;
for(int ploop=0; ploop<pole_array.size();ploop++)
{

if(pole_equiv(pole_array[ploop], pole) ){ // std::cerr<<"Duplicate pole detected!"<<std::endl;
duplicate=true;

pole_array[ploop].multiplicity_+=1;
pole_array[ploop].which_g_.push_back(pole.which_g_[0]);

break;
}

}

// added may 26 2021
// extra check if the pole is only a fermionic frequency
 
bool non_zero=false; 
if(drop_matsubara_poles){ 
int zcount=std::count(pole.eps_.begin(), pole.eps_.end(),0);
if( zcount<pole.eps_.size()){
	
	non_zero=true;
}
}else{
	non_zero=true;
}


if(non_zero){
	
// if (!duplicate && non_zero ){
if (!duplicate){
pole_array.push_back(pole);
}
}

// }


pole.alpha_.clear();
pole.eps_.clear();
pole.which_g_.clear();

}




}

return pole_array;

}



bool AmiBase::pole_equiv (pole_struct pole1, pole_struct pole2){

bool result=true;

if(pole1.eps_.size()!= pole2.eps_.size()){ return false;}
if(pole1.alpha_.size()!= pole2.alpha_.size()){ return false;}

for (int i=0; i< pole1.eps_.size(); i++)
{
if (pole1.eps_[i] != pole2.eps_[i]){ return false; break;}

}

for (int i=0; i< pole1.alpha_.size(); i++)
{
if (pole1.alpha_[i] != pole2.alpha_[i]){ return false; break;}

}



return result;
}


bool AmiBase::g_equiv (g_struct g1, g_struct g2){

bool result=true;

if(g1.eps_.size()!= g2.eps_.size()){ return false;}
if(g1.alpha_.size()!= g2.alpha_.size()){ return false;}

for (int i=0; i< g1.eps_.size(); i++)
{
if (g1.eps_[i] != g2.eps_[i]){ return false; break;}

}

for (int i=0; i< g1.alpha_.size(); i++)
{
if (g1.alpha_[i] != g2.alpha_[i]){ return false; break;}

}



return result;
}



AmiBase::sign_t AmiBase::find_signs(int index, g_prod_t &R){

sign_t sign;

  for (int i=0; i< R.size(); i++){


	if (R[i].alpha_[index] != 0){

	sign.push_back((double)R[i].alpha_[index]);
		}

 
  }

return sign;
}


double AmiBase::get_simple_sign(int index,g_prod_t &R, pole_struct pole){
	
	double sign;

 sign=R[pole.which_g_[0]].alpha_[index];

  

return sign;
	
	
	
}

/**
* Depricated version of update_gprod for only simple poles. Replaced by update_gprod_general.
*/
void AmiBase::update_gprod_simple(int index, AmiBase::R_t &R_array, AmiBase::P_t &P_array,AmiBase::S_t &S_array){

int next=index+1;

g_prod_array_t temp_g_array;
Pi_t temp_pole_array;
Si_t temp_sign_array;

	for (int j=0; j<R_array[index].size(); j++){



	pole_array_t poles;
	poles=find_poles(index, R_array[index][j]);  

	sign_t signs;
	signs=find_signs(index,R_array[index][j]);

	pole_array_t cor_poles;




		for (int i=0; i < poles.size(); i++){

		
		// append to the next R_array
			for(int m=0; m< poles[i].multiplicity_ ; m++){
			temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));
			cor_poles.push_back(poles[i]);
			}
		}
	
	temp_pole_array.push_back(cor_poles);
	temp_sign_array.push_back(signs);
	
	}


R_array.push_back(temp_g_array);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


}






/**
 *
 * Void function takes a product of Green's functions - an element of an `Ri_t` object - and obtains the residue for a specific pole. Output is an array of such elements, a full `Ri_T` object, along with the respective poles and signs for `Pi_t` and  `Si_t`. Includes derivative. 
 * @param[in] G_in : Product of green's functions
 * @param[in] pole : Pole to evaluate residue
 * @param[out] Ri_out : Resultant `g_prod_t` 
 * @param[out] poles : Resultant `pole_array_t`
 * @param[out] signs : Resultant `sign_t`. 
*/
void AmiBase::evaluate_general_residue(AmiBase::g_prod_t G_in, AmiBase::pole_struct pole, AmiBase::Ri_t &Ri_out, AmiBase::pole_array_t &poles, AmiBase::sign_t &signs){

// std::cout<<"------Starting a residue---------"<<std::endl;

double starting_sign;
starting_sign=get_starting_sign(G_in,pole);

g_prod_t W_array;
W_array=reduce_gprod(G_in,pole);  // this is h(z) in f(z)h(z)


Ri_t temp_ri, temp_ri_out, collect_ri;
pole_array_t temp_poles_out, collect_poles;
sign_t temp_signs_out, collect_signs;

temp_ri.push_back(W_array);
sign_t temp_signs;
pole_array_t temp_pole_in;
temp_signs.push_back(starting_sign);
temp_pole_in.push_back(pole);
// take derivatives
for (int m=0; m< pole.multiplicity_-1; m++){
// std::cout<<"On derivative "<<m+1<<std::endl;
// std::cout<<"temp ri has size "<<temp_ri.size()<<std::endl;
for(int i=0; i< temp_ri.size(); i++){
// First store the fermi derivative part 	
	// std::cout<<"On i="<<i<<std::endl;
	
take_derivative_gprod(temp_ri[i], temp_pole_in[i], temp_signs[i], temp_ri_out, temp_poles_out,temp_signs_out);	

// std::cout<<"Temp ri out has size "<<temp_ri_out.size()<<std::endl;
// std::cout<<"Temp poles out has size "<<temp_poles_out.size()<<std::endl;
// std::cout<<"Temp signs out has size "<<temp_signs_out.size()<<std::endl;
// for(int k=0; k< temp_poles_out.size(); k++){
// std::cout<<"Pole "<<k<<" with der="<< temp_poles_out[k].der_<<std::endl;	
	
// }

for(int add=0; add< temp_ri_out.size(); add++){
collect_ri.push_back(temp_ri_out[add]);
collect_poles.push_back(temp_poles_out[add]);
collect_signs.push_back(temp_signs_out[add]);

}

}

// std::cout<<"Temp_ri_out has size "<<temp_ri_out.size()<<std::endl;
// std::cout<<"Collect ri has size "<<collect_ri.size()<<std::endl;

temp_ri=collect_ri;// temp_ri_out;
temp_signs=collect_signs;//temp_signs_out;
temp_pole_in=collect_poles;//temp_poles_out;

collect_ri.clear();
collect_poles.clear();
collect_signs.clear();

}

signs=temp_signs;//signs_out;
poles=temp_pole_in;//s_out;
Ri_out=temp_ri;//_out;
		
	// sub in pole at the end 
	for(int i=0; i< Ri_out.size();i++){
	for(int j=0; j<Ri_out[i].size();j++){
	//std::cout<<i<<" "<<j<<std::endl;
	Ri_out[i][j]=update_G_pole(Ri_out[i][j],pole);
	}
	}
	
}

/**
* Take derivative of `g_prod_t` with respect to index stored in `pole_struct`. This is accomplished via multiple applications of chain rule.  This is described somewhat in PRB 99 035120 - though the notation changed somewhat.  The prefactor of -1 is no longer absorbed into a green's function, but instead is stored in the `S_t` array.  What is omitted there is discussion of derivatives of the fermi functions.  Here they are simply tracked by defining an integer `der_` of the `pole-struct` that gets incremented. 
*
*/
void AmiBase::take_derivative_gprod(AmiBase::g_prod_t &g_prod, AmiBase::pole_struct pole, double start_sign, AmiBase::Ri_t &r_out, AmiBase::pole_array_t &poles, AmiBase::sign_t &signs){
poles.clear();
signs.clear();
r_out.clear();	


pole_struct fermi_pole=pole;
fermi_pole.der_++;

r_out.push_back(g_prod);
poles.push_back(fermi_pole);	
signs.push_back(start_sign);

g_prod_t temp_gprod;



for(int term=0; term< g_prod.size(); term++){
int alpha=g_prod[term].alpha_[pole.index_];

if(alpha!=0){

 signs.push_back(-(double)alpha*start_sign);

	for(int m=0; m< g_prod.size(); m++){
	if(term==m){
	
	temp_gprod.push_back(g_prod[m]);
	temp_gprod.push_back(g_prod[m]);
	}
	else{
		
	temp_gprod.push_back(g_prod[m]);
	}
    		
	}
r_out.push_back(temp_gprod);
temp_gprod.clear();

	
}



}	


for(int m=0; m<signs.size()-1; m++){
poles.push_back(pole);	
}
	
}


/**
* Obtains the elements of `S_t` arrays, modified for poles with multiplicity!=1. 
*/
double AmiBase::get_starting_sign(AmiBase::g_prod_t G_in, AmiBase::pole_struct pole){

double result=1.0;

for (int i=0; i< pole.which_g_.size(); i++){

result=result*(double)G_in[pole.which_g_[i]].alpha_[pole.index_];

}

// std::cout<<"On pole with multiplicity "<<pole.multiplicity_<<std::endl;

result=result/(double)factorial(pole.multiplicity_-1 );


return result;
}

/**
* Removes green's functions from `g_prod_t` to account for numerator of residue (z-z0)^m
*/
AmiBase::g_prod_t AmiBase::reduce_gprod(AmiBase::g_prod_t G_in, AmiBase::pole_struct pole){

g_prod_t reduced;

for (int i=0; i< G_in.size(); i++){

bool add=true;

for(int j=0; j< pole.which_g_.size();j++){

if(pole.which_g_[j]==i){add=false; break;}

}

if(add){ reduced.push_back(G_in[i]);}

}

return reduced;
}







int AmiBase::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}