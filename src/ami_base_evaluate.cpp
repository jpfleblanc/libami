#include "ami_base.hpp"


std::complex<double> AmiBase::evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external){


 // std::cout<<"Evaluating Result for construction: ";

int dim=parms.N_INT_;

// std::cout<<"dim="<<dim<< std::endl;
SorF_t SorF_result;


if(dim==0){
std::complex<double> gprod;

gprod=eval_gprod(parms, R_array[0][0], external);

return gprod;

}

if (dim==1){
SorF_t SF_left, SF_right;
SorF_t S_double_left, S_double_right;	
	
SF_left=dot(S_array[0], fermi(parms,P_array[0], external));	

SorF_result=SF_left;
	
}

for (int i=0; i< dim-1; i++){

SorF_t SF_left, SF_right;
SorF_t S_double_left, S_double_right;

// need line that converts the Si_t from integers to doubles? So that the dot operator has doubles both left and right entries.



if(i==0){
// do dot operation
SF_left=dot(S_array[i], fermi(parms,P_array[i], external));
 // std::cout<<"S["<<i<<"].f(P["<<i<<"])";
}
else {SF_left=SorF_result;}

// do dot
//std::cout<<"Before i="<<i<<std::endl;
SF_right=dot(S_array[i+1], fermi(parms,P_array[i+1], external));

// std::cout<<i<<std::endl;

 SorF_result=cross(SF_left,SF_right);
 // std::cout<<"xS["<<i+1<<"].f(P["<<i+1<<"])";
 

// std::cout<<"After i "<<i<<"steps, K contains "<<std::endl;
// for( int x=0; x< SorF_result[0].size(); x++){
// std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
// }



}
 


std::complex<double> final_result;

// std::cout<<"*R for dim="<<dim<< std::endl;
final_result=star(parms, SorF_result, R_array[dim], external);


return final_result;


}

std::complex<double> AmiBase::star(ami_parms &parms, SorF_t K, Ri_t R, ami_vars external){

// std::cout<<"Size of left is "<< K.size() <<std::endl;
// std::cout<<"Size of left array is "<< K[0].size() <<std::endl;
// std::cout<<"Size of right is "<< R.size() <<std::endl;

std::complex<double> output=0;
std::complex<double> term;
std::complex<double> gprod;

std::complex<double> t0,t10;

// std::ofstream file;
// file.open("outfile.dat",  std::ofstream::out | std::ofstream::app);
bool print_output=false;

for( int i=0; i< K[0].size(); i++)
{



gprod=eval_gprod(parms, R[i], external);
term=K[0][i]*gprod;

// if(i==0){t0=term;}
// if(i==10){t10=term;}

//std::isnan(std::real(term))

// if(abs(term)<1){
output+= term;
// }
// if( true ){output+= term;

// std::cout<<"Term apparently is finite? "<< term << std::endl;
// }
// if(std::real(term)>10){
// print_g_prod_info(R[i]);
// }
// if(term.real() > 100 || term.imag()>100 ){
// std::cout<<"In star K[]*R"<<std::endl;
// std::cout<< std::setprecision(20)<< i<<" "<< K[0][i] <<" "<< std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;
// print_output=true;
// }
 // file <<i<<" "<< K[0][i] <<" "<< std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;

}

// file.close();

// std::cout<<"t0 plus t10 gives "<< t0<<" "<<t10<<" "<< t0+t10<<std::endl;

// if(print_output){
// std::cout<<"Returning value of "<< output <<std::endl;
// }

return output;

}


/// All the special function structures.


AmiBase::SorF_t AmiBase::dot(Si_t Si, SorF_t fermi){

SorF_t output;
output.reserve(Si.size());

for( int i=0; i< Si.size(); i++){

std::vector<std::complex<double> > line;
line.reserve(Si[i].size());
for(int j=0; j< Si[i].size(); j++){

// std::cout<<"Dot term "<<std::endl;
// std::cout<<Si[i][j]<<" "<<fermi[i][j]<<std::endl;
line.push_back(Si[i][j]*fermi[i][j]);

}

output.push_back(line);

}


return output;
}

AmiBase::SorF_t AmiBase::cross(SorF_t left, SorF_t right){
// should check that size() of left==1 or else this won't work. Throw an error.
// also, left[0].size()==right.size()

// std::cout<<"Size of left is "<< left.size() <<std::endl;
// std::cout<<"left[0] and right is "<< left[0].size()<< " "<< right.size() <<std::endl;

SorF_t output;


std::vector<std::complex<double> > line;
line.reserve(left[0].size()*right.back().size());// just guess that the last one is the biggest for reservation

for( int i=0; i< left[0].size(); i++){
for( int rj=0; rj< right[i].size(); rj++){
//std::cout<<"i and rj are "<<i<<" "<<rj<<std::endl ;
line.push_back(left[0][i]*right[i][rj]);

}
}


// std::cout<<"Lengths are "<< line.size()<<std::endl;

output.push_back(line);


return output;


}

AmiBase::SorF_t AmiBase::fermi(ami_parms &parms, Pi_t Pi, ami_vars external){

SorF_t output;
output.reserve(Pi.size());


for (int i=0; i< Pi.size(); i++){

std::vector<std::complex<double> > group;
group.reserve(Pi[i].size());

  for (int j=0; j< Pi[i].size(); j++){
// std::cout<<"Fermi term "<< fermi_pole(parms, Pi[i][j], external)<<std::endl;
  group.push_back( fermi_pole(parms, Pi[i][j], external));

}

output.push_back(group);


}



return output;

}


std::complex<double>  AmiBase::fermi_pole(ami_parms &parms, pole_struct pole, ami_vars external){

// std::cout<<"Working on pole"<<std::endl;
// print_pole_struct_info(pole);

std::complex<double>  output,term;
int eta=0;

double beta=external.BETA_;
double E_REG=parms.E_REG_;
// std::cout<<"Evaluating with energy regulator of "<< E_REG<<std::endl;


// In order to generalize to have fermi and bose lines, here to 'sigma' needs to be considered. 

// example fixed
//
// 1) create a stat map for the frequencies std::vector<int> stat_map: populate 1 for fermi and 0 for bose. length is equal to alpha.size() 
// 2) simply replace eta=eta + 1*stat_map[i] 
//

// alternate fix.  parms.TYPE_ is 0 for sigma, 1 for Pi etc.  So if parms.TYPE_==1 and pole.alpha_.back()==1, don't add one. else add one to eta.

// TODO: should this be improved somehow? Should frequencies know their stat type?
// handle all but one external
for (int i=0; i< pole.alpha_.size()-1; i++){
//eta+= pole.alpha_[i];
if(pole.alpha_[i]!=0){
	eta++;
}
}

// handle external based on graph type 
if(pole.alpha_.back()!=0 && parms.TYPE_!=AmiBase::Pi_phuu && parms.TYPE_!=AmiBase::Pi_phud  && parms.TYPE_ !=AmiBase::doubleocc && parms.TYPE_!=AmiBase::Pi_ppuu && parms.TYPE_!=Pi_ppud){


eta++;	

// std::cout<<"External triggered"<<std::endl;
}

// if this is a double occupancy graph then the external line is bosonic. so it is a bosonic matsubara integral. so eta needs to be incremented IF the pole is for the last integration index
if(parms.TYPE_==AmiBase::doubleocc){

if(pole.index_==pole.alpha_.size()-1){
eta++;
// std::cout<<"Incrementing eta since this is a bosonic integral for "<< pole.index_<<std::endl;
}	
	
}

// END TODO


/* 
for (int i=0; i< pole.alpha_.size(); i++){
//eta+= pole.alpha_[i];
if(pole.alpha_[i]!=0){
	eta++;
}
} */

// could put infor into ami_vars external as to what the state type of the external variables is.
std::complex<double>  E= get_energy_from_pole(pole,external);
// std::cout<<"Energy is "<<E<<std::endl;

// if(eta%2==0){
// std::cout<<"eta was even"<<std::endl;
// }else{
	// std::cout<<"eta was odd"<<std::endl;
// }


double sigma= pow(-1.0, double(eta));
// if(eta==0){sigma=-1;}

std::complex<double> zero(0,0);
std::complex<double> im(0,1);

// TODO: This regulator may not be sufficient for derivative functions.
// if(eta%2==1 && abs(E.real())<abs(E_REG)){
// E+=E_REG*sgn(E.real());
// }	

// TEST - always add the regulator
// TODO: changed on may 18 2020 for molecules - should instead reconstruct ami solutions 
// std::cout<<"DB is "<<drop_bosonic_diverge<<std::endl;

if(E==zero && sigma==-1  ){  // && pole.der_==0    Not sure if pole derivative plays any role
	// std::cout<<"Bosonic function at zero energy - must vanish, setting to zero"<<std::endl;
if(drop_bosonic_diverge){	
// std::cout<<"Bosonic function at zero energy - must vanish, setting to zero"<<std::endl;
return zero;	// TODO: need to test this might be an approximation 
}
E+=E_REG;	
}else{
	if(sgn(E.real())!=0){
	E+=E_REG*sgn(E.real());}
	else{
		E+=E_REG;
	}
}

// july 23rd adding back in a different regulator 

// if(sgn(E.real())!=0){
	// E+=E_REG*sgn(E.real());}
	// else{
		// E+=E_REG;
	// }

// if(sgn(E.real())!=0){
	// E+=E_REG*sgn(E.real());}
	// else{
		// E+=E_REG;
	// }

if(drop_der && pole.der_!=0){

// if(pole.der_!=0){
	
	return zero;
}


// if(pole.der_==0){
// output=1.0/(sigma*std::exp(beta*(E))+1.0);
// return output;
// }

int m=pole.der_;

// compute m'th derivative

// std::cout<<m<<" "<<sigma<<" "<<beta<<" "<<E<<std::endl;
output=fermi_bose(m,sigma,beta,E);


/*
if(pole.der_==0){
output=1.0/(sigma*std::exp(beta*(E))+1.0);
// return output;
}else{  // compute m'th derivative

for( int k=0; k<m+1; k++){
	// std::cout<<"On k'th derivative "<<k<<std::endl;
	term= frk(m,k)*std::exp(k*beta*(E))*std::pow(sigma,k) *std::pow(-1.0, k+1)/std::pow(sigma*std::exp(beta*(E))+1.0, k+1) ;
	output+= term;
	// std::cout<<"Fermi Pole Term evaluated to "<< term << " at energy "<< E<<" with sigma "<<sigma<< " betaE is "<< beta*E<<" in exponent "<< std::exp(beta*(E))<< std::endl;
}

// TODO: double check that this multiplication is general 
output=output*std::pow(beta,m)*(-1.0);

// output=0.0;
}
*/

// TODO I'm not sure this is correct for doubleocc 
// if external was integrated over and bosonic, then the above (-1.0) should not be there. ... I think, and the starting fermi turns into a bosonic function  should generalize this
if(parms.TYPE_==AmiBase::doubleocc){
output=-1.0*output;	
}


// std::cout<<"Evaluated Fermi derivative function and got "<<output<< " at energy "<< E<<std::endl;

return output;
}

// this computes the mth order derivative of the fermi or bose distribution functions at beta, for energy E. sigma=1.0 for fermi and -1.0 for bose. 
std::complex<double> AmiBase::fermi_bose(int m, double sigma, double beta, std::complex<double> E){

std::complex<double> output,term;
output=0.0;

if(m==0){
output=1.0/(sigma*std::exp(beta*(E))+1.0);
// return output;
}else{  // compute m'th derivative

for( int k=0; k<m+1; k++){
	// std::cout<<"On k'th derivative "<<k<<std::endl;
	// term= frk(m,k)*std::exp(k*beta*(E))*std::pow(sigma,k) *std::pow(-1.0, k+1)/std::pow(sigma*std::exp(beta*(E))+1.0, k+1) ;
	term= frk(m,k)*std::pow(sigma,k) *std::pow(-1.0, k+1)*(1.0/(sigma*std::exp(beta*(E))+1.0)/std::pow(sigma +std::exp(-beta*(E)), k)) ;
	output+= term;
	// std::cout<<"Fermi Pole Term evaluated to "<< term << " at energy "<< E<<" with sigma "<<sigma<< " betaE is "<< beta*E<<" in exponent "<< std::exp(beta*(E))<< std::endl;
}

// TODO: double check that this multiplication is general 
output=output*std::pow(beta,m)*(-1.0);
}
	
	// std::cout<<"Fermi Pole output evaluated to "<< output << " at energy "<< E<<" with sigma "<<sigma<< " betaE is "<< beta*E<<" in exponent "<< std::exp(beta*(E))<< std::endl;

return output;	
	
}




std::complex<double> AmiBase::get_energy_from_pole( pole_struct pole, ami_vars external){

std::complex<double> output(0,0);

// std::cout<<"Evaluating energies"<<std::endl;
for (int i=0; i< pole.eps_.size(); i++){
// std::cout<<"Pole "<<double(pole.eps_[i])<<" external e is "<< external.energy_[i]<<" mult is "<<double(pole.eps_[i])*external.energy_[i];
output+= double(pole.eps_[i])*external.energy_[i];

}


// std::cout<<"Output is "<<output<<std::endl;

return output;

}


std::complex<double>  AmiBase::get_energy_from_g( g_struct g, ami_vars external){

std::complex<double>  output=0;

for (int i=0; i< g.eps_.size(); i++){

output+= double(g.eps_[i])*external.energy_[i];

}



return output;

}


// TODO: Comment completely 
std::complex<double> AmiBase::eval_gprod(ami_parms &parms, g_prod_t g_prod, ami_vars external){
std::complex<double> output(0,0);

std::complex<double> denom_prod(1,0);
double prefactor=external.prefactor;
// std::cout<<"Evaluating with prefactor "<< prefactor<<std::endl;

double E_REG=parms.E_REG_;

bool verbose=false;
// int N_EXT=parms.N_EXT_;

// std::cout<<"Eval Gprod"<<std::endl;

for(int i=0; i< g_prod.size(); i++){
std::complex<double> alphadenom(0,0);
std::complex<double> epsdenom(0,0);


// TODO: This should multiply from the back and stop when it hits the end of num_ext freq. 

// precompute a list of non-zero alpha_ and eps_ vector indices, so these only loop over non-zero entries.
// TODO: this change has almost no practical speedup
/* for(int a=g_prod[i].alpha_.size()-1; a>g_prod[i].alpha_.size()-1-N_EXT; a--){
alphadenom+=double(g_prod[i].alpha_[a])*external.frequency_[a];	
} */

 for(int a=0; a< g_prod[i].alpha_.size(); a++){
alphadenom+=double(g_prod[i].alpha_[a])*external.frequency_[a];

// std::cout<<"Alpha's and frequencies " << g_prod[i].alpha_[a] <<" " << external.frequency_[a] << std::endl;

} 
// TODO: Right here, if the denom==0 still, then the R entry was empty, so regulate the next section, eps -> eps+i0+

std::complex<double> zero(0,0);
std::complex<double> im(0,1);
/* 
if(alphadenom==zero){
alphadenom+=E_REG*im+E_REG;
// std::cout<<"Added ereg in gprod_eval"<<std::endl;
// verbose=true;
// alphadenom+=E_REG;	
} */

// std::cout<<"Energies"<<std::endl;
// Unsure. should this be -=? given that my epsilon is the positive quantity?
for(int a=0; a< g_prod[i].eps_.size(); a++){
epsdenom+=double(g_prod[i].eps_[a])*external.energy_[a];
//denom-=double(g_prod[i].eps_[a])*external.energy_[a];
// if(verbose){
// std::cout<<a<<" "<<g_prod[i].eps_[a]<<" ext "<< external.energy_[a]<<" "<<std::endl;
// }
}

if(alphadenom==zero){
	// return zero;
	double val=E_REG*sgn(epsdenom.real());
alphadenom+=val+val*im;
// std::cout<<"Added ereg in gprod_eval "<<E_REG<<std::endl;
// verbose=true;
// alphadenom+=E_REG;	
}


/* 
if(verbose){
std::cout<<"Epsdenom gave "<<epsdenom<<std::endl;
}
// std::cout<<"End"<<std::endl;

if(epsdenom==zero){
	std::cout<<"Epsilon returned zero"<<std::endl;
}

if(std::abs(epsdenom)< E_REG){
	std::cout<<"Small epsilon denominator detected"<<std::endl;
}
 */
//denom+=get_energy_from_g(parms, g_prod[i], external);

//if (denom==std::complex<double>((0,0))){denom=E_REG;}//  prefactor=1.0;}//0.0;}
// std::cout<<"This denom gives term "<<alphadenom+epsdenom<<std::endl;

// if(std::abs(epsdenom)>E_REG){
denom_prod=denom_prod*(alphadenom+epsdenom);
// }

/* 
if((1.0/denom_prod).real()> 100){
	std::cout<<"Denominator product is "<< denom_prod<<" with prefactor "<< prefactor<<std::endl;
std::cout<<"returning value "<< 1.0/denom_prod*prefactor<<std::endl;

for(int a=0; a< g_prod[i].eps_.size(); a++){
// epsdenom+=double(g_prod[i].eps_[a])*external.energy_[a];
//denom-=double(g_prod[i].eps_[a])*external.energy_[a];
// if(verbose){
std::cout<<a<<" "<<g_prod[i].eps_[a]<<" ext "<< external.energy_[a]<<" "<<std::endl;
// }
}

	
}  */



}

// std::cout<<"Denominator product is "<< denom_prod<<" with prefactor "<< prefactor<<std::endl;
// std::cout<<"returning value "<< 1.0/denom_prod*prefactor<<std::endl;

output=1.0/denom_prod*prefactor;



return output;
}


// Using notation to match 10.1103/PhysRevB.101.075113 
// They produced coefficients to the fermi functions and put them in a table. 
// We derivated a general expression for those coefficients - we believe this to be general but have only checked up to 6th order I think 
// TODO: Precomputing these might be a smart thing to do 
double AmiBase::frk(int r, int k){
double output=0.0;



for(int m=0; m< k+1; m++){
	
output+= binomialCoeff(k,m)*std::pow(m,r)*(std::pow(-1,k-m));	
	
	
}

// std::cout<<"Evaluating Frk function "<<r<<" "<<k<<" = "<<output<<std::endl;
	
	return output;
	
}

// Returns value of Binomial Coefficient C(n, k)  
int AmiBase::binomialCoeff(int n, int k){  
    int res = 1;  
  
    // Since C(n, k) = C(n, n-k)  
    if ( k > n - k )  
        k = n - k;  
  
    // Calculate value of  
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]  
    for (int i = 0; i < k; ++i)  
    {  
        res *= (n - i);  
        res /= (i + 1);  
    }  
  
    return res;  
}  



