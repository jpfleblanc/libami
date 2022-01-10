//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami.hpp"
#include <iomanip>


AmiCalc::AmiCalc(ami_parms &parms)
{
 
// std::cout<<"Ami Constructor Called"<<std::endl;

}

AmiCalc::AmiCalc()
{
 
// std::cout<<"Ami Constructor Called"<<std::endl;

}

// TODO do this actually need alps::params?

void AmiCalc::construct(ami_parms &parms, g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array){

// std::cout<<"Here the code constructs the analytic solution and stores it in the R, P and S arrays"<<std::endl;
R_array.clear();
P_array.clear();
S_array.clear();
int dim=parms.N_INT_;

g_prod_array_t R0_array;
R0_array.push_back(R0);
R_array.push_back(R0_array);


pole_array_t poles;

// this will be iterated from 0 to dim exclusive


//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.


for (int index=0; index<dim; index++){
// std::cout<<"Working on integration number "<< index <<std::endl;

update_gprod_general(index, index, R_array, P_array,S_array);
// update_gprod_simple(index, R_array, P_array,S_array);

/* std::cout<<"------------------"<<std::endl;
std::cout<<"After integral "<<index<<" array sizes are"<<std::endl;
std::cout<<"------------------"<<std::endl;

std::cout<<"R["<<index+1<<"] is length"<<R_array[index+1].size()<<std::endl;

int count=0;
for(int i=0; i< S_array.size();i++){
	count+= S_array[i].size();
}
std::cout<<"S total size is "<< count<<std::endl;
 */

}






}

// attempt to get different construction (more stable) by skipping the multipole indexes 
void AmiCalc::minimal_construct(ami_parms &parms, g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array){


R_array.clear();
P_array.clear();
S_array.clear();
int dim=parms.N_INT_;

g_prod_array_t R0_array;
R0_array.push_back(R0);
R_array.push_back(R0_array);


pole_array_t poles;

// this will be iterated from 0 to dim exclusive


//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.
int array_index=0;

std::vector<int> skipped;

for (int int_index=0; int_index<dim; int_index++){

update_gprod_general_minimal(int_index, array_index, R_array, P_array,S_array);

if(R_array.size()==array_index+2){array_index++;}
else{skipped.push_back(int_index);
std::cout<<"Skipped integral "<< int_index<<" due to multiplicity "<<std::endl;
}

}

// do the skipped indexes regardless of multiplicity 
for(int index=0; index< skipped.size(); index++){
	
update_gprod_general(skipped[index], array_index, R_array, P_array,S_array);
array_index++;	
	
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
// std::cerr<<"Duplicate pole detected!"<<std::endl;
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




}

// std::cout<<"Found x poles "<< pole_array.size()<<std::endl;
// std::cout<<"With multiplicities ";
// for(int p=0; p< pole_array.size();p++){
	// std::cout<< pole_array[p].multiplicity_<<std::endl;
// }

return pole_array;

}


bool AmiCalc::pole_equiv (pole_struct pole1, pole_struct pole2){

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

// std::cout<<"Updating simple residue "<<std::endl;
// print_pole_struct_info(pole);
// print_g_prod_info(G_in);


for (int i=0; i< G_in.size(); i++){


if(i != pole.which_g_[0]){

//update_G_pole(G_in[i],pole);
residue.push_back(update_G_pole(G_in[i],pole) );


}


}
// std::cout<<"Updated Gprod is now "<<std::endl;
// print_g_prod_info(residue);


return residue;
}
///









////
AmiCalc::g_struct AmiCalc::update_G_pole(AmiCalc::g_struct g_in, AmiCalc::pole_struct pole){

AmiCalc::g_struct g_new;

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

// std::cout<<"Reduction starts with size "<< G_in.size()<<std::endl;

for (int i=0; i< G_in.size(); i++){

bool add=true;

for(int j=0; j< pole.which_g_.size();j++){

if(pole.which_g_[j]==i){add=false; break;}

}

if(add){ reduced.push_back(G_in[i]);}

}

// std::cout<<"Lengths in reduction are "<< reduced.size()<<" "<< G_in.size()<<std::endl;
// std::cout<<"For pole with "<< pole.which_g_.size()<<std::endl;
// std::cout<<"Reduction ends with size "<< reduced.size()<<std::endl;
return reduced;
}

double AmiCalc::get_starting_sign(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole){

double result=1.0;

for (int i=0; i< pole.which_g_.size(); i++){

result=result*(double)G_in[pole.which_g_[i]].alpha_[pole.index_];

}

// std::cout<<"Starting sign for m="<<pole.multiplicity_<<" is ";

result=result/(double)factorial(pole.multiplicity_-1 );

// std::cout<<result<<std::endl;

return result;
}

AmiCalc::g_struct AmiCalc::der_fix(g_struct &g_in, double alpha){
g_struct out;
out.alpha_.resize(g_in.alpha_.size());
out.eps_.resize(g_in.eps_.size());
	
for(int i=0; i< g_in.alpha_.size(); i++){
	out.alpha_[i]=-g_in.alpha_[i]*alpha;
}

for(int i=0; i< g_in.eps_.size(); i++){
	out.eps_[i]=-g_in.eps_[i]*alpha;
}

	
return out;	
}

// TODO: This is depricated 
void AmiCalc::take_derivatives(AmiCalc::Ri_t &Wi, AmiCalc::pole_struct pole, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){

// std::cout<<"Before derivatives, sign list is of size "<<signs.size()<<std::endl;

// std::cout<<"Wi has "<< Wi.size()<<std::endl;
g_prod_array_t temp_g_array;
sign_t temp_signs;
pole_array_t temp_poles;
for(int i=0; i< Wi.size(); i++){
// std::cout<<"On i loop "<< i<<std::endl;
// print_g_prod_array(Wi);

  for(int dindex=0; dindex< Wi[i].size(); dindex++){
  g_prod_t temp_gprod;

  int alpha=Wi[i][dindex].alpha_[pole.index_];
  // std::cout<<"Alpha was "<< alpha<<std::endl;
   if(alpha !=0 ){
// populate signs and poles
// std::cout<<"Adding sign into temp "<<-(double)alpha*signs[i]<<std::endl;
// std::cout<<(double)alpha<<" "<< signs[i]<<std::endl;
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

signs=temp_signs;

} // close i

// std::cout<<"Temp signs of size "<<temp_signs.size()<<std::endl;
// for(int m=0; m< temp_signs.size(); m++){
	// std::cout<<temp_signs[m]<<std::endl;
// }
// update input Wi with output.
Wi=temp_g_array;
// std::cout<<"Signs and temp are "<< signs.size()<<" "<< temp_signs.size()<<std::endl;
// std::cout<<"Poles and temp are "<< poles.size()<<" "<< temp_poles.size()<<std::endl;

poles=temp_poles;
// print_pole_array(poles);
// std::cout<<"Signs and temp are "<< signs.size()<<" "<< temp_signs.size()<<std::endl;
// std::cout<<"Poles and temp are "<< poles.size()<<" "<< temp_poles.size()<<std::endl;

}

void AmiCalc::evaluate_general_residue(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole, AmiCalc::Ri_t &Ri_out, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){

// poles.clear();
// signs.clear();

double starting_sign;
starting_sign=get_starting_sign(G_in,pole);

g_prod_t W_array;
W_array=reduce_gprod(G_in,pole);  // this is h(z) in f(z)h(z)



Ri_t temp_ri, temp_ri_out;
pole_array_t temp_poles_out;
sign_t temp_signs_out;

temp_ri.push_back(W_array);
sign_t temp_signs;
pole_array_t temp_pole_in;
temp_signs.push_back(starting_sign);
temp_pole_in.push_back(pole);
// take derivatives
for (int m=0; m< pole.multiplicity_-1; m++){

for(int i=0; i< temp_ri.size(); i++){
// First store the fermi derivative part 	
	
	
// temp_ri and temp_signs have same length
take_derivative_gprod(temp_ri[i], temp_pole_in[i], temp_signs[i], temp_ri_out, temp_poles_out,temp_signs_out);	
// std::cout<<"temp_ri_out.size() is "<<temp_ri_out.size()<<std::endl;

}

// if(temp_ri_out.size()==0){break;}

temp_ri=temp_ri_out;
temp_signs=temp_signs_out;
temp_pole_in=temp_poles_out;

}

signs=temp_signs_out;
poles=temp_poles_out;
Ri_out=temp_ri_out;

// std::cout<<"Signs poles and Ri are size "<<signs.size()<<" "<<poles.size()<<" "<<Ri_out.size()<<std::endl;

// sub in pole at the end 
for(int i=0; i< Ri_out.size();i++){
for(int j=0; j<Ri_out[i].size();j++){
//std::cout<<i<<" "<<j<<std::endl;
Ri_out[i][j]=update_G_pole(Ri_out[i][j],pole);
}
}
	
}


void AmiCalc::take_derivative_gprod(AmiCalc::g_prod_t &g_prod, AmiCalc::pole_struct pole, double start_sign, AmiCalc::Ri_t &r_out, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){
poles.clear();
signs.clear();
r_out.clear();	


pole_struct fermi_pole=pole;
fermi_pole.der_++;

// std::cout<<"Added fermi part to derivative"<<std::endl;
// fermi derivative part 
r_out.push_back(g_prod);
poles.push_back(fermi_pole);	
signs.push_back(start_sign);
////
 

g_prod_t temp_gprod;



for(int term=0; term< g_prod.size(); term++){
int alpha=g_prod[term].alpha_[pole.index_];
/* std::cout<<"Alpha is "<<alpha<<" for term "<< term<< " with start sign "<< start_sign <<std::endl;

std::cout<<"Starting G_prod for term  "<< term<<" is"<<std::endl;
print_g_prod_info(g_prod);

std::cout<<"Pole is "<<std::endl;
print_pole_struct_info(pole); */

if(alpha!=0){
	
	// std::cout<<"g_prod.size() is "<< g_prod.size()<<std::endl;

 // signs.push_back(start_sign);
 signs.push_back(-(double)alpha*start_sign);
 // std::cout<<"Pushing back sign "<< start_sign<<std::endl;
	for(int m=0; m< g_prod.size(); m++){
	if(term==m){
		// std::cout<<"Pushing back twice"<<std::endl;
	temp_gprod.push_back(g_prod[m]);
	temp_gprod.push_back(g_prod[m]);
	// g_struct temp_g=der_fix(g_prod[m], alpha);
	// temp_gprod.push_back(temp_g);
	}
	else{
		// std::cout<<"Pushing back once"<<std::endl;
	temp_gprod.push_back(g_prod[m]);
	}
    // std::cout<<"After m= "<< m<<" the temp_g_prod is "<<std::endl;
    // print_g_prod_info(temp_gprod);	
		
		
	}
r_out.push_back(temp_gprod);
temp_gprod.clear();

	
}



}	


// TODO : maybe this should be signs.size()-C, where C is defined to be the starting poles.size() ?
for(int m=0; m<signs.size()-1; m++){
poles.push_back(pole);	
}
// for(int m=0; m<signs.size(); m++){
// poles.push_back(pole);	
// }

	
 // std::cout<<"Sizes of r_out and poles and signs must match "<< r_out.size()<<" "<< poles.size()<<" "<<signs.size()	<<std::endl;
	
}


/* // in the end, Ri will hold the array of updated G's
void AmiCalc::evaluate_general_residue(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole, AmiCalc::Ri_t &Wi, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs){

double starting_sign;
starting_sign=get_starting_sign(G_in,pole);
AmiCalc::sign_t temp_signs;
temp_signs.push_back(starting_sign);

g_prod_t W_array;

// remove extra terms from Mth order residue
W_array=reduce_gprod(G_in, pole);
// Get overall starting sign and prefactor

std::cout<<"Lengths are "<< W_array.size()<<" "<< G_in.size()<<std::endl;
std::cout<<"Working with pole mult="<<pole.multiplicity_<<std::endl;
std::cout<<"With starting sign="<<temp_signs[0]<<std::endl;




// put W_array into Wi to start
Wi.push_back(W_array);

// take M-1 derivatives

for(int i=0; i< pole.multiplicity_-1; i++){
// std::cout<<"Evaluating multipole for M= "<< pole.multiplicity_ <<std::endl;
// std::cout<<"Length before is "<<Wi.size()<<std::endl;
take_derivatives( Wi, pole, poles, temp_signs); 

 std::cout<<"Number of signs after derivative is "<<signs.size()<<std::endl;
}

// update each term in Ri with the pole value
print_pole_array(poles);

for(int i=0; i< Wi.size();i++){
for(int j=0; j<Wi[i].size();j++){
//std::cout<<i<<" "<<j<<std::endl;
Wi[i][j]=update_G_pole(Wi[i][j],pole);
}
}

signs=temp_signs;

} */


void AmiCalc::update_gprod_general_minimal(int int_index, int array_index,  AmiCalc::R_t &R_array, AmiCalc::P_t &P_array, AmiCalc::S_t &S_array){

//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.

//int next=index+1;

// std::cout<<"Size of R["<< index<<"]="<< R_array[index].size() <<std::endl;


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
//p_out=poles;

pole_array_t cor_poles;
sign_t col_signs;
//std::cout<<"Testing pole equivalence "<< pole_equiv(poles[0],poles[0]) <<std::endl;



for (int i=0; i < poles.size(); i++){

if(poles[i].multiplicity_==1){ 

g_in.push_back(simple_residue(R_array[array_index][j],poles[i]));

s_out.push_back(get_simple_sign(int_index,R_array[array_index][j], poles[i]));

p_out.push_back(poles[i]);

}else
{
// if any multiplicity is found, then do nothing 
return;

}
// append to the next R_array
//temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));


}
//poles.clear();

// for(int p=0; p< p_out.size(); p++){
// for(int m=0; m< p_out[p].multiplicity_; m++){	
// cor_poles.push_back(p_out[p]);	
	
// }	
// }
// print_pole_array(p_out);

temp_pole_array.push_back(p_out);
temp_sign_array.push_back(s_out);

p_out.clear();
s_out.clear();
// p_out.clear();
// cor_poles.clear();
// s_out.clear();
//temp_g_array.push_back(g_in);

}




R_array.push_back(g_in);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


// std::cout<<"Size of R["<< next<<"]="<< R_array[next].size() <<std::endl;

}


void AmiCalc::update_gprod_general(int int_index, int array_index, AmiCalc::R_t &R_array, AmiCalc::P_t &P_array, AmiCalc::S_t &S_array){

//TODO: for each GProd in R, for each pole in that Gprod, find residue for each and collect those in an array.

//int next=index+1;

// std::cout<<"Size of R["<< index<<"]="<< R_array[index].size() <<std::endl;


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
//p_out=poles;

pole_array_t cor_poles;
sign_t col_signs;
//std::cout<<"Testing pole equivalence "<< pole_equiv(poles[0],poles[0]) <<std::endl;



for (int i=0; i < poles.size(); i++){

// g new g_prod
// std::cout<<"Printing pole "<<i<<std::endl;
// print_pole_struct_info(poles[i]);

// std::cout<<"working on pole number "<< i<<std::endl;
if(poles[i].multiplicity_==1){ 
// std::cout<<"Pole "<<i<<" has multiplicity 1"<<std::endl;
// print_pole_struct_info(poles[i]);

// print_g_prod_array(R_array[index]);

g_in.push_back(simple_residue(R_array[array_index][j],poles[i]));

s_out.push_back(get_simple_sign(int_index,R_array[array_index][j], poles[i]));
// std::cout<<"Pushed back sign of "<< get_simple_sign(index,R_array[index][j], poles[i])<<std::endl;
// s_out=find_signs(index,R_array[index][j]);
p_out.push_back(poles[i]);

}else
{
	// std::cout<<"Pole "<<i<<" has multiplicity "<< poles[i].multiplicity_<<std::endl;
Ri_t temp_g;

evaluate_general_residue(R_array[array_index][j], poles[i], temp_g, temp_p, temp_s);

// std::cout<<"Lengths of these should be the same "<< temp_g.size()<<" "<< temp_p.size()<<" "<< temp_s.size()<<std::endl;
//a.insert(std::end(a), std::begin(b), std::end(b));

g_in.insert(std::end(g_in), std::begin(temp_g), std::end(temp_g));
s_out.insert(std::end(s_out),std::begin(temp_s), std::end(temp_s));
p_out.insert(std::end(p_out),std::begin(temp_p), std::end(temp_p));


// std::cout<<"temp_g returned with size "<< temp_g.size()<<std::endl;



// std::cout<<"Finished general residue "<<std::endl;
// print_pole_array(p_out);


}
// append to the next R_array
//temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));


}
//poles.clear();

// for(int p=0; p< p_out.size(); p++){
// for(int m=0; m< p_out[p].multiplicity_; m++){	
// cor_poles.push_back(p_out[p]);	
	
// }	
// }
// print_pole_array(p_out);

temp_pole_array.push_back(p_out);
temp_sign_array.push_back(s_out);

p_out.clear();
s_out.clear();
// p_out.clear();
// cor_poles.clear();
// s_out.clear();
//temp_g_array.push_back(g_in);

}




R_array.push_back(g_in);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


// std::cout<<"Size of R["<< next<<"]="<< R_array[next].size() <<std::endl;

}

double AmiCalc::get_simple_sign(int index,AmiCalc::g_prod_t &R,AmiCalc::pole_struct pole){
	
	double sign;

 sign=R[pole.which_g_[0]].alpha_[index];

  

return sign;
	
	
	
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

pole_array_t cor_poles;




for (int i=0; i < poles.size(); i++){

// g new g_prod


// append to the next R_array
for(int m=0; m< poles[i].multiplicity_ ; m++){
temp_g_array.push_back(simple_residue(R_array[index][j],poles[i]));
cor_poles.push_back(poles[i]);
}
}
//poles.clear();



temp_pole_array.push_back(cor_poles);
temp_sign_array.push_back(signs);
// cor_poles.clear();
// signs.clear();
}


R_array.push_back(temp_g_array);
P_array.push_back(temp_pole_array);
S_array.push_back(temp_sign_array);


// std::cout<<"Size of R["<< next<<"]="<< R_array[next].size() <<std::endl;

}

