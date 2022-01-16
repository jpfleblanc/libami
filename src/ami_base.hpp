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
#include <iomanip>
#include <algorithm>
#include <functional> 





/**
 * @class AmiBase
 *
 *  
 *
 * @brief  The primary class of libami. 
 *
 * @note See https://github.com/jpfleblanc/libami
 *
 * @author JPF LeBlanc 
 *
 * @version Revision: 0.61 
 *
 * @date Date: 2022/01/11  
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

/// Returns the sign of a value - or zero if it is uniquely zero. 
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Maximum argument allowed for fermi and bose functions - to prevent inf due to double precision numbers
double exp_max_arg=500.0;


/// Simple factorial function - Nothing special. 
int factorial(int n);

// AmiVars

bool drop_bosonic_diverge=false; // if set to true, then E_REG not needed because bosonic divergences will simply be set to zero.  Rigorously this may not be correct. 
bool drop_der=false; // if set to true, then E_REG not needed because bosonic divergences will simply be set to zero.  Rigorously this may not be correct.
bool drop_matsubara_poles=false; // if set to true, ignores Matsubara poles with zero energy 
// bool is_real_external=false;

// External list of energies and frequencies 
/// The energy of each denominator will always appear as a linear combination of these initial (pre integration) energies, \f$\epsilon_1, \epsilon_2\f$ ..etc  
/// By convention, the energy_t contains the NEGATIVE of the energy of a given Green's function line, \f$ 1/(X+E) \f$ where \f$ E=-\epsilon \f$.
typedef std::vector<std::complex<double>> energy_t;

/// This is the list of internal and external frequencies values.  Typically only the last elements for external frequencies are non-zero - but one can evaluate intermediate steps where multiple external frequencies are non-zero. 
///
typedef std::vector< std::complex<double>  > frequency_t;


// Fundamental objects

// the symbolic epsilon

/// Vector of type `int` with elements \f$ a_i\f$.  This is the symbolic representation of the energy \f$E=-\sum\limits_{i}a_i\epsilon_i\f$ described in AMI paper (https://doi.org/10.1103/PhysRevB.99.035120).  We use the convention \f$G=\frac{1}{X+\epsilon}\f$.  It is the coefficients for a linear combination of a set of possible values. 
typedef std::vector<int> epsilon_t;


/// Vector of type `int` with elements \f$ \alpha_i\f$.  This is the symbolic representation of the frequency, as a linear combination of possible entries.  Typically contains only values of 0, -1 and +1. Other values at intermediate steps typically represent an error.  \f$X=\sum\limits_{i} i\nu_i \alpha_i\f$. 
typedef std::vector<int> alpha_t;

/// Indicator for multi-species Green's function or energy dispersions (spin-up vs spin-dn, multiband, etc).  Technically this is not relevant for libami, but may be useful. 
typedef int species_t;

/// Indicator for statistics. A future version might use this more frequently.  Current version presumes all integration frequencies are Fermionic. 
typedef enum {Bose,Fermi} stat_type ;

// Ideally these types will not appear in the ami_base class 
/// Graph types will be removed in a future release.  Current support is limited to Sigma and Pi_phuu graph types. 
typedef enum {Sigma,Pi_phuu, Pi_phud,Hartree, Bare, Greens, density, doubleocc, Pi_ppuu, Pi_ppud, DOS,ENERGY, FORCE} graph_type ;

/// To be removed in a future release
typedef enum {hubbard,coulomb} int_type;
/// To be removed in a future release
typedef enum {tb, fp, hf} disp_type;
/// To be removed in a future release
typedef enum {matsubara, real} ext_type;

ext_type ext_freq_type=matsubara;



/// The `ami_vars` struct is the basic information required for the evaluation stage of AMI result.  It contains the variable internal/external quantities.  Specifically it is a list of numerical values for energies of each line and values for each frequency.  Also stored is the possibility of an overall prefactor. Also required is a value of \f$\beta=\frac{1}{k_B T}\f$ needed for evaluation of Fermi/Bose distributions. 
struct ami_vars{

ami_vars(energy_t eps, frequency_t freq){
energy_= eps;
frequency_= freq;
prefactor=1.0;

}

ami_vars(energy_t eps, frequency_t freq, double Bta){
energy_= eps;
frequency_= freq;
BETA_=Bta;

}

ami_vars(energy_t eps, frequency_t freq, double Bta, double pf){
energy_= eps;
frequency_= freq;
prefactor=pf;
BETA_=Bta;

}


ami_vars(){prefactor=1.0;}

/// Numerical values of energies. 
energy_t energy_;
/// Numerical Values of frequencies.
frequency_t frequency_;
/// Overall prefactor - default(1).
double prefactor;
/// Required value of inverse temperature, \f$\beta\f$.
double BETA_;

///Experimental parameter for spectral representation. 
double gamma_=0;

};


/// Parameters for AMI construction/evaluation.
struct ami_parms{
ami_parms(int N_INT,  double E_REG){
N_INT_=N_INT;
E_REG_=E_REG;
N_EXT_=1;
TYPE_=static_cast<AmiBase::graph_type>(0); /// by default sigma.
int_type_=static_cast<AmiBase::int_type>(0); /// by default is hubbard. 
dispersion_=static_cast<AmiBase::disp_type>(0);/// by default is tight-binding.
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

/// Number of integrations to perform.
int N_INT_;
/// Hardcoded as (1) in this version - represents number of external variables. 
int N_EXT_=1;
/// Possible energy regulator for evaluation of Fermi/Bose functions to control divergences.  Should be zero by default and switched on if needed.
double E_REG_=0;

graph_type TYPE_;
int_type int_type_;
disp_type dispersion_;

};



/// A Green's function structure.  This is a symbolic vector of `epsilon_t` and vector of `alpha_t`.  Also needed is what the statistics of the line are.  While this could be determined from alpha - it is better to store it.  For multistate systems the species_ index might be useful.

struct g_struct{
	
/// Constructor with state_type specified.	
g_struct(epsilon_t eps, alpha_t alpha, stat_type stat){
eps_=eps;
alpha_=alpha;
stat_=stat;
species_=0;

}

/// Constructor assumes fermi statistics if not specified for partially initialized structure.
g_struct(epsilon_t eps, alpha_t alpha){
eps_=eps;
alpha_=alpha;
stat_=Fermi;
species_=0;

}

/// Uninitialized variant.
g_struct(){
	
stat_=Fermi;
species_=0;
} 

/// Symbolic coefficients for a linear combination of epsilons.
epsilon_t eps_;

/// Symbolic coefficients for a linear combination of frequencies.
alpha_t alpha_;
/// Mark Fermi/Bose stats - Fermi is required in current version.
stat_type stat_;
species_t species_;

// Experimental Spectral representation.  Implementation incomplete.
int pp=-1;  // pp=0 means this G represents a principle part integral. pp=1 it is a delta function. else it is inert 

};

/// Pole structure. Equivalent to `g_struct`, but kept separate. Tracks multiplicity, and which Green's function it is attached to. Also it tracks how many derivatives to take when evaluated at a fermi function.
struct pole_struct{

pole_struct(epsilon_t eps, alpha_t alpha){
eps_=eps;
alpha_=alpha;
}

pole_struct(){}

epsilon_t eps_; 
alpha_t alpha_;
/// Index that specifies which frequency it is a pole with respect to.
int index_; 
/// The multiplicity of the pole, starts at 1 and increments as needed.
int multiplicity_=1;
int der_=0; /**< Counter for derivatives. */
std::vector<int> which_g_; /**< Index to identify which `g_struct` a pole originated from.*/

/// Experimental component of Spectral evaluation.
alpha_t x_alpha_;

};




/// Basic element of the S array in SPR notation.
typedef std::vector<double> sign_t;
/// Basic element of the R array in SPR notation.
typedef std::vector<g_struct> g_prod_t;
/// Basic element of the P array in SPR notation.
typedef std::vector<pole_struct> pole_array_t;

typedef std::vector<sign_t> sign_array_t;
typedef std::vector<g_prod_t> g_prod_array_t;



/// R result after each integration step.
typedef std::vector<g_prod_t> Ri_t;
/// P result after each integration step.
typedef std::vector<pole_array_t> Pi_t;
/// S result after each integration step.
typedef std::vector<sign_t> Si_t;

/// final output R arrays. 
typedef std::vector<Ri_t> R_t; 
/// final output P arrays.
typedef std::vector<Pi_t> P_t;
/// final output S arrays.
typedef std::vector<Si_t> S_t;

// typedefs for evaluation

/** Term Structure for term-by-term evaluation.  Conceptually simpler than SPR construction.
* Storage translates to \f$ \prod{f(p_i)}\prod{G_j}\times sign \f$.
*
*/
struct term{
	
	term(){
		
	}
	
term(double s, pole_array_t p, g_prod_t g){
sign=s;
p_list=p;
g_list=g;
}	

/// Sign prefactor 
double sign=1;
/// List of poles, \f$ \prod{f(p_i)}\f$. 
pole_array_t p_list;
/// List of Green's functions, \f$ \prod{G_j}\f$.
g_prod_t g_list;	
	
};

/// The storage for the term-by-term construction.  Each `term` struct is an element of the `terms` vector.
typedef std::vector< term > terms;




/// Construction function for term-by-term construction. 
void construct(int N_INT, g_prod_t R0, terms &terms_out);
/// For user simplicity this is a wrapper function to make ami_terms and SPR calls similar.
void construct(ami_parms &parms,  g_prod_t R0, terms &terms_out);
/// Evaluate Terms. 
std::complex<double> evaluate(ami_parms &parms, terms &ami_terms, ami_vars &external);
/// Evaluate a single Term.  Usage is identical to `evaluate` function. 
std::complex<double> evaluate_term(ami_parms &parms, term &ami_term, ami_vars &external);
/// Evaluate the numerator of a term - product of fermi/bose functions. 
std::complex<double> eval_fprod(ami_parms &parms,pole_array_t &p_list, ami_vars &external);

/// Integrates a single Matsubara index. 
void integrate_step(int index, terms &in_terms, terms &out_terms);

void split_term(term &this_term, pole_struct this_pole, term &innert_part, term &active_part);

/// Primary residue function for term construction. 
void terms_general_residue(term &this_term, pole_struct this_pole, terms &out_terms);
/// Derivative for term construction 
void take_term_derivative(term &in_term, pole_struct &pole, terms &out_terms);

/// Screen IO for debugging. 
void print_terms(terms &t);
/// Screen IO for debugging.
void print_term(term &t);

/// Convert terms to Ri structure for optimization.  This does not create a usable `terms` object. It is purely an intermediate step for the optimization functions.
void convert_terms_to_ri(terms &ami_terms, Ri_t &Ri);




/**
Although perhaps strangely named - represents the multiplication of sign and pole arrays defined in Equation (19) (https://doi.org/10.1103/PhysRevB.99.035120).
*/
typedef std::vector< std::vector<std::complex<double> > > SorF_t;



// Functions for poles
/*
Given an array of Green's functions, finds all poles with respect to frequency index, and checks for multiplicities. Stores multiplicities of the poles as well as `which_g_`, an identifier that specifies which Green's function it was attached to.  
*/
pole_array_t find_poles(int index, g_prod_t &R);

// Residue Functions
g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);


g_struct update_G_pole(g_struct g_in,pole_struct pole);
pole_struct update_Z_pole(AmiBase::pole_struct p_in, AmiBase::pole_struct pole);

// Functions for signs
sign_t find_signs(int index, g_prod_t &R);
double get_simple_sign(int index,g_prod_t &R,pole_struct pole);

// Functions for R's P's and S's
// TODO: Of these three only 'update_gprod_general' is actually used.  'Simple' is the case with no multi-poles.  
void update_gprod_simple(int index, R_t &R_array, P_t &P_array, S_t &S_array);

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

/// Cross and dot are other operations like 'star' above that are defined after equation (20) in the AMI paper  (https://doi.org/10.1103/PhysRevB.99.035120).
SorF_t cross(SorF_t left, SorF_t right);
SorF_t dot(Si_t Si, SorF_t fermi);

// Testing Priority: 2 should be easy to test 
/// Given a set of external energies, beta, and frequencies, will evaluate the energy of a pole_struct. 
std::complex<double> get_energy_from_pole(pole_struct pole, ami_vars external);

/// Given a set of external energies, beta, and frequencies, will evaluate the energy of a g_struct. 
std::complex<double>  get_energy_from_g(g_struct g, ami_vars external);

/// Evaluates a product of Green's functions. 
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
 * @param[in] pole1 First pole you want to compare of type `pole_struct`.
 * @param[in] pole2 Second pole you want to compare of type `pole_struct`.
 * 
*/
bool pole_equiv (pole_struct pole1, pole_struct pole2);

/**
 *
 * Checks if two `g_struct` have similar characteristics.
 * These characteristics include `epsilon_t` size and `alpha_t` size.
 * Additionally, it checks each value of `epsilon_t` vector  and `alpha_t` vector are the same.
 * @param[in] g1 First Green's function you want to compare of type `g_struct`.
 * @param[in] g2 Second Green's function you want to compare of type `g_struct`.
 * 
*/
bool g_equiv (g_struct g1, g_struct g2);

// Testing Priority: 1
// evaluate_general_residue is the primary function called in the main loop by 'update_gprod_general'

void evaluate_general_residue(g_prod_t G_in, pole_struct pole, Ri_t &Ri_out, pole_array_t &poles, sign_t &signs);


void take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

// This is actually a pretty important function. probably needs a more clear name and documentation as to what it does 
// TODO: This became depricated - unsure how 
// the der_fix function absorbed a minus sign into the alpha and epsilon values of one of the two Green's functions before pushing it into the new array 
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
In order to take arbitrary derivatives of fermi functions - we need to know the prefactor.  The fermi/bose function is given in general by the function: 
\f[

\sum\limits_{k=0}^{m} f_{rk}(m,k) \sigma^k (-1.0)^{k+1} \frac{1}{\sigma e^{\beta E}+1.0} \frac{1}{[\sigma + e^{-\beta E}]^k}
\f]

The frk function itself returns 
\f[
f_{rk}(m,k)=\sum\limits_{m=0}^{k+1} binomialCoeff(k,m) m^r (-1)^{k-m}.
\f]

*/
double frk(int r, int k);



///Recursive construction of fermi_bose derivatives. 
std::complex<double> fermi_bose(int m, double sigma, double beta, std::complex<double> E);



/*
Recursive binomial coefficient function 
*/
int binomialCoeff(int n, int k);



// Functions 

  ///Default Constructor.  Constructor is empty.  Currently no initialization is required in most cases.
  AmiBase();
  //Possibly deprecated constructor. 
  /// Constructor with ami_parms. 
  AmiBase(ami_parms &parms);


  // The construction
void construct(ami_parms &parms,  g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array);


/* 
This is the primary evaluation which takes again `ami_parms`, the outputs from `construct` as well as the `ami_vars` external values that enter into the expression 
*/

std::complex<double> evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external);


// Optimization functions - not strictly necessary but might implement automatically 


typedef std::pair<int,int> ref_t;
typedef std::vector<ref_t> ref_v_t;
typedef std::vector<ref_v_t> R_ref_t;
typedef R_ref_t ref_eval_t;

/// This is an optimized version of the evaluate function. For simplicity if the new objects are empty the evaluate function is called directly. 
std::complex<double> evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external,g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list);

void factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list);
void reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list);

std::complex<double> optimized_star(ami_parms &parms, SorF_t K, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list, ami_vars external);

bool pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign);
bool g_struct_equiv(g_struct &g1, g_struct &g2, int &sign);
void print_g_struct_info(g_struct g);
void print_epsilon_info(epsilon_t eps);
void print_alpha_info(alpha_t alpha);
void print_pole_struct_info(pole_struct g);
////

// experimental
void derivative_opt(g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list);


// term optimization

void factorize_terms(terms &ami_terms, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list);
std::complex<double> evaluate(ami_parms &parms, terms &ami_terms, ami_vars &external, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list);


private:

};