
#pragma once
#include "ami_base.hpp"
 
#include <string> 
#include <sstream> 
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <experimental/filesystem>


/**
 * @class AmiSpec
 *
 *  
 *
 * @brief  A spectral application of the AmiBase class.  Experimental.
 *
 * @note N/A
 *
 * @author JPF LeBlanc 
 *
 * @version Revision: 0.4 
 *
 * @date Date: 2020/11/03  
 *
 * 
 * Contact: jleblanc@mun.ca
 *
 *
 *
 *
 */

class AmiSpec 
{

public:

AmiBase amibase;

double xi_cutoff=10;

// Define structures that are not in the ami_base class 

typedef AmiBase::epsilon_t X_t ;


// the A_struct stores the symbolic information about A.  you will need your own function to take the A_struct and some self energy and return a value. 
struct A_struct{
A_struct(AmiBase::epsilon_t eps, X_t x){
eps_=eps;
x_=x;
species_=0;
}

/// uninitialized variant
A_struct(){
} 

// epsilon_i either by index or by value 
int eps_index=-1;
std::complex<double> eps_val;

AmiBase::epsilon_t eps_;
X_t x_;
AmiBase::alpha_t alpha_;
AmiBase::species_t species_;


};

typedef std::vector<A_struct> A_prod_t;


typedef AmiBase::pole_struct delta_t;
typedef std::vector< delta_t > delta_prod_t;
typedef AmiBase::energy_t xi_t;

struct ami_sp_term{

ami_sp_term(AmiBase::term this_term, A_prod_t aprod, delta_prod_t dprod){
aprod_=aprod;
ami_term_=this_term;
	
}
	
//uninitialized 	
ami_sp_term(){}	

A_prod_t aprod_;
AmiBase::term ami_term_;
delta_prod_t dprod_;
	
};

typedef std::vector<ami_sp_term> sp_terms;




std::complex<double> A_eval(A_struct &A, std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E);
std::complex<double> get_X(X_t &Xsym, xi_t &xi);
std::complex<double> get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps);


AmiSpec(double xc);



private:


};


