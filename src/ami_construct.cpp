//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami.hpp"
#include <iomanip>


AmiCalc::AmiCalc(ami_parms &parms)
{
 
std::cout<<"Ami Constructor Called"<<std::endl;

}

AmiCalc::AmiCalc()
{
 
std::cout<<"Ami Constructor Called"<<std::endl;

}

// TODO do this actually need alps::params?

void AmiCalc::construct(ami_parms &parms, g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array){

// std::cout<<"Here the code constructs the analytic solution and stores it in the R, P and S arrays"<<std::endl;

int dim=parms.N_INT_;

g_prod_array_t R0_array;
R0_array.push_back(R0);
R_array.push_back(R0_array);


pole_array_t poles;

// this will be iterated from 0 to dim exclusive


//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.


for (int index=0; index<dim; index++){
// std::cout<<"Working on integration number "<< index <<std::endl;

//update_gprod_general(index, R_array, P_array,S_array);
update_gprod_simple(index, R_array, P_array,S_array);

}






}




  //everything from here on down is evaluation.
void AmiCalc::evaluate(ami_parms& parms){
  

  // std::cout<<"Here the code will evaluate a given solution"<<std::endl;


}



AmiCalc::pole_array_t AmiCalc::find_poles(int index, AmiCalc::g_prod_t &R){

pole_array_t pole_array;
pole_struct pole;

// TODO: Optimize by predefining sizes.
//int size=R[0].alpha_.size();
//int epssize=R[0].eps_.size();

pole.alpha_.reserve(R[0].alpha_.size());
pole.eps_.reserve(R[0].eps_.size());

pole.index_=index;

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

/// So now you have a pole... need to check if it is already in the array, and then if it isn't add it.  if it is, modify the multiplicity of the pole in the array it is equal to. and add it's which_g_ index to the pole_array with multiplicity

bool duplicate=false;
for(int ploop=0; ploop<pole_array.size();ploop++)
{

if(pole_equiv(pole_array[ploop], pole) ){
duplicate=true;

pole_array[ploop].multiplicity_+=1;
pole_array[ploop].which_g_.push_back(pole.which_g_[0]);
std::cerr<<"Duplicate pole detected!"<<std::endl;
break;
}

}

if (duplicate==false){
pole_array.push_back(pole);
}



///


pole.alpha_.clear();
pole.eps_.clear();
pole.which_g_.clear();

}

//std::cout<<"Found x poles "<< pole_array.size()<<std::endl;


}



return pole_array;

}


bool AmiCalc::pole_equiv(pole_struct pole1, pole_struct pole2){

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



AmiCalc::sign_t AmiCalc::find_signs(int index, AmiCalc::g_prod_t &R){

sign_t sign;



  for (int i=0; i< R.size(); i++){


	if (R[i].alpha_[index] != 0){

	sign.push_back((double)R[i].alpha_[index]);
		}

  
 
  }



return sign;
}


AmiCalc::g_prod_t AmiCalc::simple_residue(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole){

g_prod_t residue;

if (pole.multiplicity_!=1){std::cerr<<"Simple residue called for pole with M!=1"<<std::endl;}
for (int i=0; i< G_in.size(); i++){


if(i != pole.which_g_[0]){

//update_G_pole(G_in[i],pole);
residue.push_back(update_G_pole(G_in[i],pole) );

}


}


return residue;
}
///









////
AmiCalc::g_struct AmiCalc::update_G_pole(AmiCalc::g_struct g_in, AmiCalc::pole_struct pole){

struct g_struct g_new;

for (int i=0; i< g_in.alpha_.size(); i++){

if (i != pole.index_){
g_new.alpha_.push_back( g_in.alpha_[i] + g_in.alpha_[pole.index_]*pole.alpha_[i]   );
}
else{
g_new.alpha_.push_back(0);
}

}

// TODO: Can these two loops be combined into a single one?

for (int i=0; i< g_in.eps_.size(); i++){
g_new.eps_.push_back( g_in.eps_[i] + g_in.alpha_[pole.index_]*pole.eps_[i]   );
}


if ( g_new.alpha_.size() != g_in.alpha_.size()){std::cerr<<"Maybe something wrong? Alphas not the same size" <<std::endl; }


return g_new;

}


AmiCalc::g_prod_t AmiCalc::reduce_gprod(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole){

g_prod_t reduced;

for (int i=0; i< G_in.size(); i++){

bool add=true;

for(int j=0; j< pole.which_g_.size();j++){

if(pole.which_g_[j]==i){add=false; break;}

}

if(add){ reduced.push_back(G_in[i]);}

}

// std::cout<<"Lengths in reduction are "<< reduced.size()<<" "<< G_in.size()<<std::endl;
// std::cout<<"For pole with "<< pole.which_g_.size()<<std::endl;

return reduced;
}

double AmiCalc::get_starting_sign(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole){

double result=1.0;

for (int i=0; i< pole.which_g_.size(); i++){

result=result*(double)G_in[pole.which_g_[i]].alpha_[pole.index_]/(double)factorial(pole.multiplicity_-1 );

}

return result;
}


void AmiCalc::take_derivatives(AmiCalc::Ri_t &Wi, AmiCalc::pole_struct pole, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){

g_prod_array_t temp_g_array;
sign_t temp_signs;
pole_array_t temp_poles;
for(int i=0; i< Wi.size(); i++){

print_g_prod_array(Wi);

  for(int dindex=0; dindex< Wi[i].size(); dindex++){
  g_prod_t temp_gprod;

  int alpha=Wi[i][dindex].alpha_[pole.index_];
  // std::cout<<"Alpha was "<< alpha<<std::endl;
   if(alpha !=0 ){
// populate signs and poles
	temp_signs.push_back(-(double)alpha*signs[i]);
        temp_poles.push_back(pole);
//
	for(int term=0; term< Wi[i].size(); term++){
	if(term==dindex){
	temp_gprod.push_back(Wi[i][term]);
	temp_gprod.push_back(Wi[i][term]);
	}
	else{
	temp_gprod.push_back(Wi[i][term]);
	}


	}//close for
    temp_g_array.push_back(temp_gprod);
    // std::cout<<"temp array has length "<< temp_g_array.size()<<std::endl;
    
    

    }//close if


  }// close dindex



} // close i

// update input Wi with output.
Wi=temp_g_array;
signs=temp_signs;
poles=temp_poles;

}

// in the end, Ri will hold the array of updated G's
void AmiCalc::evaluate_general_residue(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole, AmiCalc::Ri_t &Wi, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){


g_prod_t W_array;

// remove extra terms from Mth order residue
W_array=reduce_gprod(G_in, pole);
// Get overall starting sign and prefactor

//std::cout<<"Lengths are "<< W_array.size()<<" "<< G_in.size()<<std::endl;

double starting_sign;
starting_sign=get_starting_sign(G_in,pole);
signs.push_back(starting_sign);

// put W_array into Wi to start
Wi.push_back(W_array);

// take M-1 derivatives

for(int i=0; i< pole.multiplicity_-1; i++){
// std::cout<<"Evaluating multipole for M= "<< pole.multiplicity_ <<std::endl;
// std::cout<<"Length before is "<<Wi.size()<<std::endl;
take_derivatives( Wi, pole, poles, signs); 
// std::cout<<"Length after is "<<Wi.size()<<std::endl;
}

// update each term in Ri with the pole value


for(int i=0; i< Wi.size();i++){
for(int j=0; j<Wi[i].size();j++){
//std::cout<<i<<" "<<j<<std::endl;
Wi[i][j]=update_G_pole(Wi[i][j],pole);
}
}



}


void AmiCalc::update_gprod_general(int index, AmiCalc::R_t &R_array, AmiCalc::P_t &P_array, AmiCalc::S_t &S_array){

//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.

int next=index+1;

// std::cout<<"Size of R["<< index<<"]="<< R_array[index].size() <<std::endl;


Ri_t g_in;


g_prod_array_t temp_g_array;
Pi_t temp_pole_array;
Si_t temp_sign_array;


for (int j=0; j<R_array[index].size(); j++)
{

/*


pole_array_t poles;
poles=find_poles(index, R_array[index][j]);  

sign_t signs;
signs=find_signs(index,R_array[index][j]);

temp_pole_array.push_back(poles);
temp_sign_array.push_back(signs);
*/

sign_t s_out;
pole_array_t p_out;

pole_array_t poles;
poles=find_poles(index, R_array[index][j]);  
p_out=poles;

//std::cout<<"Testing pole equivalence "<< pole_equiv(poles[0],poles[0]) <<std::endl;



for (int i=0; i < poles.size(); i++){

// g new g_prod

//std::cout<<"working on pole number "<< i<<std::endl;
if(poles[i].multiplicity_==1){ 
g_in.push_back(simple_residue(R_array[index][j],poles[i]));
s_out=find_signs(index,R_array[index][j]);
}else
{

evaluate_general_residue(R_array[index][j], poles[i], g_in, p_out, s_out);

}
// append to the next R_array
//temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));

}
//poles.clear();

temp_pole_array.push_back(p_out);
temp_sign_array.push_back(s_out);
//temp_g_array.push_back(g_in);

}




R_array.push_back(g_in);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


// std::cout<<"Size of R["<< next<<"]="<< R_array[next].size() <<std::endl;

}


void AmiCalc::update_gprod_simple(int index, AmiCalc::R_t &R_array, AmiCalc::P_t &P_array,AmiCalc::S_t &S_array){

//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.

int next=index+1;

// std::cout<<"Size of R["<< index<<"]="<< R_array[index].size() <<std::endl;


g_prod_array_t temp_g_array;
Pi_t temp_pole_array;
Si_t temp_sign_array;


for (int j=0; j<R_array[index].size(); j++)
{



pole_array_t poles;
poles=find_poles(index, R_array[index][j]);  

sign_t signs;
signs=find_signs(index,R_array[index][j]);

temp_pole_array.push_back(poles);
temp_sign_array.push_back(signs);



for (int i=0; i < poles.size(); i++){

// g new g_prod


// append to the next R_array
temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));

}
//poles.clear();

}


R_array.push_back(temp_g_array);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


// std::cout<<"Size of R["<< next<<"]="<< R_array[next].size() <<std::endl;

}
