#include "ami_base.hpp"

// TODO: since there is some duplication of code in these two evaluate commands they should be separated
std::complex<double> AmiBase::evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external,g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list){

	// std::cout<<"Evaluating with external"<<std::endl;

	// for(int i=0; i< external.energy_.size(); i++){
		// std::cout<<external.energy_[i]<<" ";
	// }
	// std::cout<<std::endl;


	// for(int i=0; i< external.frequency_.size(); i++){
		// std::cout<<external.frequency_[i]<<" ";
	// }
	// std::cout<<std::endl;

if(Rref.size()==0|| unique_g.size()==0||Eval_list.size()==0){
	return evaluate(parms,R_array, P_array, S_array, external);
}

// std::complex<double> double_check=evaluate(parms,R_array, P_array, S_array, external);

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

// std::cout<<"For dim=1 K contains "<<std::endl;
// for( int x=0; x< SorF_result[0].size(); x++){
// std::cout<<x<<" "<< SorF_result[0][x]<<std::endl;
// }

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
// final_result=star(parms, SorF_result, R_array[dim], external);

final_result=optimized_star(parms, SorF_result, unique_g, Rref,Eval_list, external);


// std::cout<<"In AMI result got "<<final_result<<std::endl;

// if(double_check!=final_result){
	// std::cout<<"Double checking got"<<std::endl;
	// std::cout<<double_check<<" "<<final_result<<std::endl;

// }

return final_result;




}


std::complex<double> AmiBase::optimized_star(ami_parms &parms, SorF_t K, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list, ami_vars external){


std::complex<double> output=0;
std::complex<double> term;
// std::complex<double> gprod;


std::vector< std::complex<double>> unique_vals;
// std::cout<<"External prefactor is "<<external.prefactor<<std::endl;
for(int i=0; i< unique_g.size(); i++){
// we pretend each G is a g_prod and use the same function as the normal star operation for evaluating products of G's
g_prod_t this_term;
this_term.push_back(unique_g[i]);

unique_vals.push_back(eval_gprod(parms,this_term,external)*external.prefactor); // This removes the overall prefactor for each g. we then add that back later

}


// std::cout<<"Debugging Energy values "<<std::endl;

// for(int i=0; i< external.energy_.size(); i++){

	// std::cout<<external.energy_[i]<<std::endl;

// }


// std::cout<<"Debugging Unique G's "<<std::endl;
// for(int i=0; i<unique_vals.size(); i++){
	// std::cout<<i<<std::endl;
	// print_g_struct_info(unique_g[i]);

// }

// std::cout<<"Debugging unique values "<<std::endl;
// for(int i=0; i<unique_vals.size(); i++){
	// std::cout<<i<<" "<<unique_vals[i]<<std::endl;

// }

// now I have the numerical values of our unique G's

bool verbose=0;

for(int i=0; i< Eval_list.size(); i++){

	std::complex<double> ksum(0,0);
	for(int j=0; j< Eval_list[i].size(); j++){
		// std::cout<<j<<" "<<K[0][Eval_list[i][j].first]<<" "<<double(Eval_list[i][j].second)<<std::endl;
		ksum+=K[0][Eval_list[i][j].first]*double(Eval_list[i][j].second);

	}
	// in principle, every entry in eval_list has the same Rref terms
	ref_v_t pair_vec= Rref[i]; // just grab the first one
	std::complex<double> this_gprod(1,0);

	// std::cout<<"Term comprised of unique G indexes ";
	for(int j=0; j< pair_vec.size(); j++){
	// std::cout<<pair_vec[j].first<<" ";
	this_gprod=this_gprod*unique_vals[pair_vec[j].first];

	}

	// std::cout<<std::endl;

	term=ksum*this_gprod*external.prefactor; // add back the overall prefactor for this term


	// std::cout<<"In optimized star K[]*R"<<std::endl;
// std::cout<< std::setprecision(20)<< i<<" "<< ksum <<" "<< std::real(this_gprod)<<" "<<std::imag(this_gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;

	output+=term;


}





return output;

}



 /**
 * This is a primary AMI symbolic evaluation function.  It takes a solution defined by S, P, and R arrays and the external parameters and returns a `complex<double>` restult.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] R_array: Input `R_t`
 * @param[in] P_array : Input `P_t`
 * @param[in] S_array : Input `S_t`
 * @param[in] external : Input `ami_vars` containing point to evaluate.
*/
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

output+= term;

// std::cout<<"In star K[]*R"<<std::endl;
// std::cout<< std::setprecision(20)<< i<<" "<< K[0][i] <<" "<< std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<< std::imag(term) <<std::endl;
// print_output=true;
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


line.push_back(Si[i][j]*fermi[i][j]);

}

output.push_back(line);

}


return output;
}

AmiBase::SorF_t AmiBase::cross(SorF_t left, SorF_t right){

SorF_t output;


std::vector<std::complex<double> > line;
line.reserve(left[0].size()*right.back().size());// just guess that the last one is the biggest for reservation

for( int i=0; i< left[0].size(); i++){
for( int rj=0; rj< right[i].size(); rj++){
//std::cout<<"i and rj are "<<i<<" "<<rj<<std::endl ;
line.push_back(left[0][i]*right[i][rj]);

}
}


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

// Spectral evaluation only 
std::complex<double> freq_shift(0,0);
if(pole.x_alpha_.size()!=0){
	for(int i=0; i< pole.x_alpha_.size(); i++){
		
		freq_shift-=external.frequency_[i]*(double)pole.x_alpha_[i];
		
	}
	
}
// std::cout<<"Freq shift is "<<freq_shift<<std::endl;


// std::complex<double> freq_shift(0,0);
// if(is_real_external){

	// if(pole.alpha_.back()!=0){

		// freq_shift=external.frequency_.back()*(double)pole.alpha_.back();
		// pole.alpha_.back()=0;

	// }


// }


// In order to generalize to have fermi and bose lines, here to 'sigma' needs to be considered.

// example fixed
//
// 1) create a stat map for the frequencies std::vector<int> stat_map: populate 1 for fermi and 0 for bose. length is equal to alpha.size()
// 2) simply replace eta=eta + 1*stat_map[i]
//

// alternate fix.  parms.TYPE_ is 0 for sigma, 1 for Pi etc.  So if parms.TYPE_==1 and pole.alpha_.back()==1 (or -1), don't add one. else add one to eta.

// TODO: should this be improved somehow? Should frequencies know their stat type?
// handle all but one external
for (int i=0; i< pole.alpha_.size()-1; i++){
//eta+= pole.alpha_[i];
if(pole.alpha_[i]!=0){
	eta++;
}
}



// handle external based on graph type
if(pole.alpha_.back()!=0 && parms.TYPE_!=AmiBase::Pi_phuu && parms.TYPE_!=AmiBase::Pi_phud  && parms.TYPE_ !=AmiBase::doubleocc && parms.TYPE_!=AmiBase::Pi_ppuu && parms.TYPE_!=AmiBase::Pi_ppud && parms.TYPE_!=AmiBase::FORCE){


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


// could put infor into ami_vars external as to what the state type of the external variables is.
std::complex<double>  E= get_energy_from_pole(pole,external);

// std::cout<<"Energy and shift "<< E<<" "<<freq_shift<<std::endl;

// In the case of spectral poles the freq_shift might not be zero 
E=E+freq_shift;



// if(std::abs(E.real())<E_REG){return std::complex<double>(0,0);}

// std::cout<<"Energy is "<<E<<std::endl;

// if(eta%2==0){
// std::cout<<"eta was even"<<std::endl;
// }else{
	// std::cout<<"eta was odd"<<std::endl;
// }


double sigma= pow(-1.0, double(eta));
// if(eta==0){sigma=-1;}

// std::cout<<"eta and sigma "<< eta<<" "<<sigma<<std::endl;

std::complex<double> zero(0,0);
std::complex<double> im(0,1);

// TODO: This regulator may not be sufficient for derivative functions.
// if(eta%2==1 && abs(E.real())<abs(E_REG)){
// E+=E_REG*sgn(E.real());
// }


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

// std::cout<<"Fermi Pole Term evaluated to "<< term << " at energy "<< E<<" with sigma "<<sigma<< " betaE is "<< beta*E<<" in exponent "<< std::exp(beta*(E))<< std::endl;


// TODO I'm not sure this is correct for doubleocc
// if external was integrated over and bosonic, then the above (-1.0) should not be there. ... I think, and the starting fermi turns into a bosonic function  should generalize this
if(parms.TYPE_==AmiBase::doubleocc){
output=-1.0*output;
}


// std::cout<<"Evaluated Fermi derivative function and got "<<output<< " at energy "<< E <<" with der="<<pole.der_<<" and sigma="<<sigma<<std::endl;

// std::cout<<"Fermi Pole Term evaluated to "<< output << " at energy "<< E<<" with der="<<pole.der_<<" with sigma "<<sigma<< " betaE is "<< beta*E<<" in exponent "<< std::exp(beta*(E))<< std::endl;

return output;
}

/// This computes the mth order derivative of the fermi or bose distribution functions at beta, for energy E. sigma=1.0 for fermi and -1.0 for bose.
std::complex<double> AmiBase::fermi_bose(int m, double sigma, double beta, std::complex<double> E){

std::complex<double> output,term;
output=0.0;



 // if(sigma==-1){
	// E=std::complex<double>(std::abs(E.real()), E.imag());
	// if(E.real()<0 ){
		// std::complex<double> zero(0,0);
		// return zero;
	// }

// }

if(m==0){
	
	double arg=std::real(beta*E);
	double arg_amp=std::abs(arg);
	
	if(arg>exp_max_arg){
		double arg_sign=(double)sgn(arg);
			if(arg_sign>0){
				output=0;
			}
			else{
				output=1;
			}
			
		
				
		
	}else{
output=1.0/(sigma*std::exp(beta*(E))+1.0);
	}
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

// Evaluating energies for pole
for (int i=0; i< pole.eps_.size(); i++){

output+= double(pole.eps_[i])*external.energy_[i];

}

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


 for(int a=0; a< g_prod[i].alpha_.size(); a++){
alphadenom+=double(g_prod[i].alpha_[a])*external.frequency_[a];

}
// TODO: Right here, if the denom==0 still, then the R entry was empty, so regulate the next section, eps -> eps+i0+

std::complex<double> zero(0,0);
std::complex<double> im(0,1);


// std::cout<<"Energies"<<std::endl;
// Unsure. should this be -=? given that my epsilon is the positive quantity?
for(int a=0; a< g_prod[i].eps_.size(); a++){
epsdenom+=double(g_prod[i].eps_[a])*external.energy_[a];
//denom-=double(g_prod[i].eps_[a])*external.energy_[a];
// if(verbose){
// std::cout<<a<<" "<<g_prod[i].eps_[a]<<" ext "<< external.energy_[a]<<" "<<std::endl;
// }
}

// For safety this is disabled.  No regulator whatsoever.
// if(alphadenom==zero){
	// return zero;
	// double val=E_REG*sgn(epsdenom.real());
	// alphadenom+=val*im;
// alphadenom+=E_REG*sgn(epsdenom.real())+E_REG*sgn(epsdenom.imag())*im;
// alphadenom+=E_REG*sgn(epsdenom.real())+E_REG*im;
// std::cout<<"Added ereg in gprod_eval "<<alphadenom<<std::endl;
// verbose=true;
// alphadenom+=E_REG;
// }


denom_prod=denom_prod*(alphadenom+epsdenom);
// std::cout<<"alphadenom:  "<<alphadenom<<"  epsdenom:  "<<epsdenom<<"  prefactor:  "<<prefactor<<std::endl;


}


output=1.0/denom_prod*prefactor;



return output;
}


/// Using notation to match 10.1103/PhysRevB.101.075113
/// They produced coefficients to the fermi functions and put them in a table.
/// We derive a general expression for those coefficients - we believe this to be general but have only checked up to 6th order I think
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
