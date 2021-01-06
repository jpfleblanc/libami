/*
 * Copyright (C) 2018 JPF LeBlanc jleblanc@mun.ca  See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */


#pragma once

 
#include <string> 
#include <sstream> 
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>





/**
 * @class AmiBase
 *
 *  
 *
 * @brief  The primary class of libami. 
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

class AmiBase
{
	
public:

/// Not to be underappreciated - returns the sign of a value - or zero if it is uniquely zero. 
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


/// Simple factorial function - Nothing special. 
int factorial(int n);

// AmiVars

bool drop_bosonic_diverge=false; // if set to true, then E_REG not needed because bosonic divergences will simply be set to zero.  Rigorously this may not be correct. 


// External list of energies and frequencies 
/// The energy of each denominator will always appear as a linear combination of these initial (pre integration) energies, \f$\epsilon_1, \epsilon_2\f$ etc  
typedef std::vector<std::complex<double>> energy_t;

/// This is the list of independent and external frequencies.  Typically only the last elements for external frequencies are non-zero - but one can evaluate intermediate steps where multiple external frequencies are necessary. 
///
typedef std::vector< std::complex<double>  > frequency_t;


// Fundamental objects

// the symbolic epsilon

/// Vector of type `int`.  This is the symbolic \f$-\epsilon\f$ described in Ref[AMI].  We use the convention \f$G=\frac{1}{\omega+\epsilon}\f$.  
typedef std::vector<int> epsilon_t;


/// Vector of type `int`.  This is the symbolic representation of the frequency, as a linear combination of possible entries.  Typically contains only values of 0, -1 and +1. Other values at intermediate steps typically represent an error.  \f$X=\sum\limits_{i} i\nu_i \alpha_i\f$. 
typedef std::vector<int> alpha_t;

/// Indicator for multi-species green's function or energy dispersions (spin-up vs spin-dn, multiband, etc)
typedef int species_t;

/// Indicator for statistics. A future version might use this more frequently.  Current version presumes all integration frequencies are Fermionic. 
typedef enum {Bose,Fermi} stat_type ;

// Ideally these types will not appear in the ami_base class 
// Graph types are an unfortunate oversight on my part - and should not strictly be in AMI 
typedef enum {Sigma,Pi_phuu, Pi_phud,Hartree, Bare, Greens, density, doubleocc, Pi_ppuu, Pi_ppud, DOS,ENERGY} graph_type ;

typedef enum {hubbard,coulomb} int_type;
typedef enum {tb, fp, hf} disp_type;


/// the ami_vars struct if the basic information required for the evaluation stage of AMI result.  Specifically it is a list of numerical values for energies of each line and values for each frequency.  Also stored is the possibility of an overall prefactor. Also required is a value of \f$\beta=\frac{1}{k_B T}\f$ needed for evaluation of Fermi/Bose distributions. 

struct ami_vars{

ami_vars(energy_t eps, frequency_t freq){
energy_= eps;
frequency_= freq;
prefactor=1.0;

}

ami_vars(){prefactor=1.0;}

energy_t energy_;
frequency_t frequency_;
double prefactor;
double BETA_;

};


// TODO: does ami_parms need TYPE_, int_type_, dispersion_ etc? 
// TODO: Define new objects that satisfy the requirements of TYPE_ etc 
// TODO: is E_REG used during construction?  I think not. 

/// Parameters for the construction phase of AMI.
struct ami_parms{
ami_parms(int N_INT,  double E_REG){
N_INT_=N_INT;
E_REG_=E_REG;
N_EXT_=1;
TYPE_=static_cast<AmiBase::graph_type>(0); /// by default sigma
int_type_=static_cast<AmiBase::int_type>(0); /// by default is hubbard 
dispersion_=static_cast<AmiBase::disp_type>(0);// by default is tight-binding.
}

ami_parms(int N_INT, double E_REG, graph_type TYPE){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=TYPE;
int_type_=static_cast<AmiBase::int_type>(0);
N_EXT_=1;
}

ami_parms(int N_INT, double E_REG, graph_type TYPE, int_type inter, disp_type disp){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=TYPE;
N_EXT_=1;
int_type_=static_cast<AmiBase::int_type>(inter);
dispersion_=static_cast<AmiBase::disp_type>(disp);
}

ami_parms(){}

int N_INT_;
int N_EXT_;
double E_REG_;// does E_Reg belong here? 


// these should be replaced with something more specific to AMI 
graph_type TYPE_;
int_type int_type_;
disp_type dispersion_;

};



/// A Green's function structure.  This is a symbolic vector of `epsilon_t` and vector of `alpha_t`.  Also needed is what the statistics of the line are.  While this could be determined from alpha - it is better to store it.  For multistate systems the species_ index might be useful.

struct g_struct{
g_struct(epsilon_t eps, alpha_t alpha, stat_type stat){
eps_=eps;
alpha_=alpha;
stat_=stat;
species_=0;

}

/// Constructor assumes fermi statistics if not specified for partially initialized structure
g_struct(epsilon_t eps, alpha_t alpha){
eps_=eps;
alpha_=alpha;
stat_=Fermi;
species_=0;

}

/// uninitialized variant
g_struct(){} 

epsilon_t eps_;
// std::vector<int> eps_indices_;
alpha_t alpha_;
stat_type stat_;
species_t species_;

};

/// Pole structure. Equivalent to `g_struct`, but kept separate. Tracks multiplicity, and which green's function it is attached to. Also it tracks how many derivatives to take when evaluated at a fermi function.

struct pole_struct{

epsilon_t eps_; 
alpha_t alpha_;
int index_; 
int multiplicity_=1;
int der_=0; /**< Count of derivatives */
std::vector<int> which_g_; /**< G-structure that is attached to it */

};



/// AMI output datatypes

typedef std::vector<double> sign_t;
typedef std::vector<g_struct> g_prod_t;
typedef std::vector<pole_struct> pole_array_t;

typedef std::vector<sign_t> sign_array_t;
typedef std::vector<g_prod_t> g_prod_array_t;



// S, P , R at each step

typedef std::vector<g_prod_t> Ri_t;
typedef std::vector<pole_array_t> Pi_t;
typedef std::vector<sign_t> Si_t;

// final output arrays

typedef std::vector<Ri_t> R_t; 
typedef std::vector<Pi_t> P_t;
typedef std::vector<Si_t> S_t;

// typedefs for evaluation

/*
Though perhaps strangely named - represents the multiplication of sign and pole arrays defined in Equation (19)
*/
typedef std::vector< std::vector<std::complex<double> > > SorF_t;



// Functions for poles
//Testing Priority: 1 - this is an essential function - if it fails nothing will work 
/*
Given an array of Green's functions, finds all poles with respect to frequency index, and checks for multiplicities. Stores multiplicities of the poles as well as `which_g_`, an identifier that specifies which green's function it was attached to.  
*/
pole_array_t find_poles(int index, g_prod_t &R);

// Residue Functions
// TODO: Simple_residue might be depricated - check usage
g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);

// Testing Priority: 1 - updating G with a given pole is an essential function - if it fails nothing will work 
g_struct update_G_pole(g_struct g_in,pole_struct pole);

// Functions for signs
sign_t find_signs(int index, g_prod_t &R);
double get_simple_sign(int index,g_prod_t &R,pole_struct pole);

// Functions for R's P's and S's
// TODO: Of these three only 'update_gprod_general' is actually used.  'Simple' is the case with no multi-poles.  
void update_gprod_simple(int index, R_t &R_array, P_t &P_array, S_t &S_array);

// Testing Priority: 1 - This is the central loop of the code. However it is a very complicated function so again test the internal functions and not the function itself 
void update_gprod_general(int int_index, int array_index,  R_t &R_array, P_t &P_array, S_t &S_array);

// Functions for Evaluation
// Testing Priority: 1 - these should all be testable - and form the backbone of the evaluation 
// This is the star function from AMI paper below equation 20: Define the function here and reference the paper: prb 99 035120
std::complex<double> star(ami_parms &parms, SorF_t H, Ri_t R, ami_vars external);

// This function evaluates equation (20) - it takes an array of poles and evaluates the fermi function for each pole and keeps the array structure 
SorF_t fermi(ami_parms &parms, Pi_t pi, ami_vars external);

// This is the central evaluation of the fermi and bose functions.  It also includes evaluating arbitrary derivatives of the functions.  See frk function that is rather complicated .  This function is also the MOST challenging function for numerical evaluation.  It is the most likely source of issue or errors and some thought should go into testing this carefully
//Testing Priority: 1
std::complex<double>  fermi_pole(ami_parms &parms, pole_struct pole, ami_vars external);

// Cross and dot are other operations like 'star' above that are defined after equation (20) in the AMI paper 
SorF_t cross(SorF_t left, SorF_t right);
SorF_t dot(Si_t Si, SorF_t fermi);

// Testing Priority: 2 should be easy to test 
std::complex<double> get_energy_from_pole(pole_struct pole, ami_vars external);
std::complex<double>  get_energy_from_g(g_struct g, ami_vars external);

std::complex<double> eval_gprod(ami_parms &parms, g_prod_t g_prod, ami_vars external);

/* math example 

 * The distance between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$.
 * \f[
    |I_2|=\left| \int_{0}^T \psi(t) 
             \left\{ 
                u(a,t)-
                \int_{\gamma(t)}^a 
                \frac{d\theta}{k(\theta,t)}
                \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
             \right\} dt
          \right|
 * \f]
 */


/**
 *
 * Checks if two Poles have similar characteristics.
 * These characteristics include `epsilon_t` size and `alpha_t` size.
 * Additionally, it checks each value of `epsilon_t` vector  and `alpha_t` vector are the same.
 * @param[in] pole1 First pole you want to check and of type `pole_struct`
 * @param[in] pole2 Second pole you want to check and of type `pole_struct`
 * 
*/

bool pole_equiv (pole_struct pole1, pole_struct pole2);

/// Testing Priority: 1
// evaluate_general_residue is the primary function called in the main loop by 'update_gprod_general'

void evaluate_general_residue(g_prod_t G_in, pole_struct pole, Ri_t &Ri_out, pole_array_t &poles, sign_t &signs);


void take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

// This is actually a pretty important function. probably needs a more clear name and documentation as to what it does 
// TODO: This became depricated - unsure how 
// the der_fix function absorbed a minus sign into the alpha and epsilon values of one of the two green's functions before pushing it into the new array 
// I think i bailed on amir's definition and instead put a minus sign into the sign array - much cleaner 
// g_struct der_fix(g_struct &g_in, double alpha);

/**
* As part of the construction of the `Si_t` arrays, the initial sign is extracted prior to taking derivatives.  Also is multiplied by \f$1/(M-1)!\f$ where M is the `multiplicity_` of the respective pole. 
*
*/
double get_starting_sign(g_prod_t G_in, pole_struct pole);

//This function removes the inert parts of the gprod in the context of taking derivatives 
g_prod_t reduce_gprod(g_prod_t G_in, pole_struct pole);


// derivatives of fermi functions
/**
In order to take arbitrary derivatives of fermi functions - we need to know the prefactor.  It is given in general by this function.  TODO: This could be replaced with an algorithmic differentiation tool. 


\f[

\sum\limits_{k=0}^{m} frk(m,k)*std::exp(k*beta*(E))*std::pow(sigma,k) *std::pow(-1.0, k+1)/std::pow(sigma*std::exp(beta*(E))+1.0, k+1)


frk function itself returns 
for(int m=0; m< k+1; m++){
	
output+= binomialCoeff(k,m)*std::pow(m,r)*(std::pow(-1,k-m));	

}
\f]

*/
double frk(int r, int k);

/*
Recursive binomial coefficient function 
*/
int binomialCoeff(int n, int k);



// Functions 

  ///Default Constructor
  AmiBase();
  /// Constructor with ami_parms
  AmiBase(ami_parms &parms);


  // The construction
void construct(ami_parms &parms,  g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array);


/* 
This is the primary evaluation which takes again `ami_parms`, the outputs from `construct` as well as the `ami_vars` external values that enter into the expression 
*/

std::complex<double> evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external);




private:

};