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
#include <experimental/filesystem>
//namespace std::filesystem = std::experimental::filesystem;

int factorial(int n);
// this is a template for the signum function to produce +1, -1 or 0 as the sign 
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


/**
 * @class AmiCalc
 *
 *  LIBAMI
 * (Note, this needs exactly one somewhere)
 *
 * @brief Provide an example
 *
 * This class is meant as an example.  It is not useful by itself
 * rather its usefulness is only a function of how much it helps
 * the reader.  It is in a sense defined by the person who reads it
 * and otherwise does not exist in any real form. 
 *
 * @note Attempts at zen rarely work.
 *
 * @author JPF LeBlanc 
 *
 * @version Revision: 1.5 
 *
 * @date Date: 2020/09/01  
 *
 * 
 * Contact: jleblanc@mun.ca
 *
 * Created on: testdate later
 *
 *
 */


class AmiCalc 
{

public:


///----- Hartree fock functions are part of epsilon - and these should wind up external and be removed eventually 
// HARTREE FOCK ENERGY OBJECTS and functions

void read_hf(std::string pgrid_filename, std::string ptoi_filename, std::string sigma_filename);
std::vector<int> ptoi;
std::vector<double> pgrid, sigma_hf;
double hf_mu, hf_kstep;
double rs;

double get_hf_sigma(double kk);
double hf_energy(double kk);


// Similar to Hartree fock - these are not strictly part of AMI and should be removed 
// molecule functions 

void read_hii(std::string filename, int maxval);
std::vector<std::complex<double>> global_hii;
bool not_molecule=1; // When removed check the instance of this variable 

bool density_warning=true;



// Type definitions for AMI

// Fundamental objects

// the symbolic epsilon

/// Vector of type `int`
typedef std::vector<int> epsilon_t;
/// Vector of type `int`
typedef std::vector<int> alpha_t;

typedef double hopping_t;
typedef std::vector<hopping_t> hopping_list_t;
typedef int species_t;

typedef enum {Bose,Fermi} stat_type ;

// Graph types are an unfortunate oversight on my part - and should not strictly be in AMI 
typedef enum {Sigma,Pi_phuu, Pi_phud,Hartree, Bare, Greens, density, doubleocc, Pi_ppuu, Pi_ppud} graph_type ;

typedef enum {hubbard,coulomb} int_type;
typedef enum {tb, fp, hf} disp_type;

// TODO: Check what this does 
bool flatten_warning=true;

// AMI Parameter structure
// This is a little weird.  The ami_parms structure is defined here. But there is no internal variable set for this... meaning it has to be passed to all of the functions... probably this should be changed. Explicit passing is 'ok'.

struct ami_parms{
ami_parms(int N_INT,  double E_REG){
N_INT_=N_INT;
E_REG_=E_REG;
N_EXT_=1;
TYPE_=static_cast<AmiCalc::graph_type>(0); /// by default sigma
int_type_=static_cast<AmiCalc::int_type>(0); /// by default is hubbard 
dispersion_=static_cast<AmiCalc::disp_type>(0);// by default is tight-binding.
}

ami_parms(int N_INT, double E_REG, graph_type TYPE){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=TYPE;
int_type_=static_cast<AmiCalc::int_type>(0);
N_EXT_=1;
}

ami_parms(int N_INT, double E_REG, graph_type TYPE, int_type inter, disp_type disp){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=TYPE;
N_EXT_=1;
int_type_=static_cast<AmiCalc::int_type>(inter);
dispersion_=static_cast<AmiCalc::disp_type>(disp);
}

ami_parms(){}

int N_INT_;
int N_EXT_;
double E_REG_;
graph_type TYPE_;

int_type int_type_;
disp_type dispersion_;

};


// Green's function structure which has a symbolic epsilon and set of alpha values
struct g_struct{
g_struct(epsilon_t eps, alpha_t alpha, stat_type stat){
eps_=eps;
alpha_=alpha;
stat_=stat;
species_=0;

}

// Assume fermi statistics if not specified for partially initialized structure
g_struct(epsilon_t eps, alpha_t alpha){
eps_=eps;
alpha_=alpha;
stat_=Fermi;
species_=0;

}

g_struct(){} // uninitialized variant

epsilon_t eps_;
// std::vector<int> eps_indices_;
alpha_t alpha_;
stat_type stat_;
species_t species_;

};

/// Pole structure. Equivalent to G structure, but kept separate. Tracks multiplicity, and which green's function it is attached to. Also it tracks how many derivatives to take when evaluated at a fermi function.

struct pole_struct{

epsilon_t eps_; 
alpha_t alpha_;
int index_; 
int multiplicity_=1;
int der_=0; /**< Count of derivatives */
std::vector<int> which_g_; /**< G-structure that is attached to it */

};


// External Variable structure and typedefs

/// external list of energies used with epsilon_t to construct energies.

typedef std::vector<std::complex<double>> energy_t;
typedef std::vector< std::complex<double>  > frequency_t;

typedef std::vector< double> k_vector_t;
typedef std::vector< k_vector_t> k_vect_list_t;



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

typedef std::vector<g_prod_array_t> R_t;  // TODO: is there a reason this is not a vector of Ri_t ???
typedef std::vector<Pi_t> P_t;
typedef std::vector<Si_t> S_t;

// typedefs for evaluation

typedef std::vector< std::vector<std::complex<double> > > SorF_t;


// various structures

// internal state - list of internal k'sb_type
// ext_vars and external_variable_list : loaded from a file and populated with external parameters to be evaluated
// ami_vars - actual input to ami routines for constructing solutions

struct internal_state{

internal_state(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
if(2*order_-1>0){
t_list_.resize(2*order_-1,1);}
disp_=static_cast<AmiCalc::disp_type>(0); // by default tight binding unless necessary to change 

mink_=0;
maxk_=2*M_PI;

}

internal_state(){}

void initialize(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
if(2*order_-1>0){
t_list_.resize(2*order_-1,1);}

disp_=static_cast<AmiCalc::disp_type>(0); // by default tight binding unless necessary to change 
mink_=0;
maxk_=2*M_PI;
	
}

k_vect_list_t internal_k_list_;
int order_;
int dim_;

hopping_list_t t_list_;
disp_type disp_;

double mink_,maxk_;

};


struct ext_vars{

ext_vars(int dim){
	
KDIM_=dim;
// external_k_vector_.assign(dim,0.0);	
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);
external_freq_.resize(1);


	
}
ext_vars(int dim, double beta, std::complex<double> mu){
	
KDIM_=dim;
// external_k_vector_.assign(dim,0.0);
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);	
external_freq_.resize(1);

BETA_=beta;
MU_=mu;

H_=0;
	
}

ext_vars(int dim, double beta, std::complex<double> mu, double H){
	
KDIM_=dim;
// TODO: make this a vector of external momentum vectors 
// external_k_vector_.assign(dim,0.0);	
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);
external_freq_.resize(1);

BETA_=beta;
MU_=mu;

H_=H;
	
}

ext_vars(){
	external_freq_.resize(1);
	
}


// k_vector_t external_k_vector_;
k_vector_t dummy_k_;
k_vect_list_t external_k_list_;
	
int KDIM_;
frequency_t external_freq_;

double BETA_;
double H_;
std::complex<double> MU_;
};

typedef std::vector< ext_vars> external_variable_list;



// ami_vars should be constructed from ext_vars and internal_state

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

typedef std::vector<ami_vars> ami_vars_list;


// The solution set structure is the complete package that is passed for evaluation 

struct solution_set{
	
solution_set(g_prod_t test_R0,S_t  test_S, P_t test_P, R_t test_R,ami_parms test_amiparms, double prefactor){

R0_=test_R0;
S_=test_S;
P_=test_P;
R_=test_R;
ami_parms_=test_amiparms;
prefactor_=prefactor;


}

solution_set(g_prod_t test_R0,S_t  test_S, P_t test_P, R_t test_R,ami_parms test_amiparms, double prefactor, std::vector<alpha_t> bose){

R0_=test_R0;
S_=test_S;
P_=test_P;
R_=test_R;
ami_parms_=test_amiparms;
prefactor_=prefactor;

bose_alphas_=bose;

}

solution_set(){}
	
g_prod_t R0_;
S_t S_;
P_t P_;
R_t R_;
ami_parms ami_parms_;
double prefactor_;
int loops_;

std::vector<alpha_t> bose_alphas_;

	
};


// I'm not sure if this is actually used. 
// TODO: this is likely depricated 
// TODO: this is used by amicalc codes 
struct evaluation_set{
evaluation_set(){}

evaluation_set(solution_set s,	ext_vars ext){
ext_vars_=ext;
sol_=s;	
}
	
	
ext_vars ext_vars_;
solution_set sol_;
	
};


/// These don't explicitly need to exist in ami library - but not harmful 
typedef std::vector<solution_set> solution_set_vec_t;
typedef std::vector< std::vector<solution_set> > solution_set_matrix_t;
typedef std::vector< solution_set_matrix_t > gg_solution_set_matrix_t;




// IO function for external variables

//TODO: The external variable list is hardcoded for what I need.  A user might want something very different. So their external_variable_list might contain completely different things.  Some thought into how to restructure this with epsilon is warranted

void read_external(std::string filename, external_variable_list &extern_list);

// At one point we thought we would write S,P,R to files along with prefactors etc. Largely this was abbandoned but technically these work.
// Test Priority: 3
void read_text_S_solutions(std::string filename, S_t &s_array);
void read_text_P_solutions(std::string eps_filename,std::string alpha_filename, P_t &p_array);
void read_text_R_solutions(std::string eps_filename,std::string alpha_filename, R_t &r_array, int size);
void read_text_R0(std::string alpha_filename, std::string eps_filename, g_prod_t &R0);
double load_prefactor(std::string filename, std::string mul_file, int order);
double load_mul(std::string filename);
void write_S_readable(S_t &s_array);
void write_P_readable(P_t &p_array);
void write_R_readable(R_t &r_array);

void load_solutions(std::string top_directory, solution_set_matrix_t &AMI_MATRIX, int MAX_ORDER, double EREG);

/// Evaluation
// evaluate functions 
// TODO: Can't recall why there are two. Found it cleaner for some reason to separate the real and imaginary measurements. 
/// Test Priority: 1 - however, this is a very complex function So rather than testing it directly, go into the function and find the functions it depends on and test those. I dont' think we can test the whole thing. 
void evaluate_solutions(std::vector<std::complex<double>> &results, solution_set &AMI, ami_vars_list &ami_eval_vars_list);
void evaluate_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars);


void construct_ami_vars_list(AmiCalc::g_prod_t &R0, double prefactor, AmiCalc::internal_state &state, AmiCalc::external_variable_list &external,AmiCalc::ami_vars_list &vars_list);
ami_vars construct_ami_vars(AmiCalc::g_prod_t &R0, double prefactor, AmiCalc::internal_state &state, AmiCalc::ext_vars &external);

// Energy and epsilon functions remain problematic
// testing priority: 4
energy_t construct_energy(AmiCalc::g_prod_t &R0, AmiCalc::internal_state &state, AmiCalc::ext_vars &external);
std::complex<double> eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, std::complex<double> mu  , disp_type disp );
std::complex<double> eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, disp_type disp);
k_vector_t construct_k(AmiCalc::alpha_t alpha, AmiCalc::k_vect_list_t &k);



/// Various functions

g_prod_t construct_multipole_example();
g_prod_t construct_example();
g_prod_t construct_example_Y();
g_prod_t construct_example_J();


// These are depricated
ami_vars construct_ext_example_Y();
ami_vars construct_ext_example_J();
ami_vars construct_ext_example_Sigma();
ami_vars construct_ext_multipole_example();
ami_vars construct_6ord_ext_multipole_example();
ami_vars construct_4ord_ext_multipole_example();

ami_vars construct_random_example_J(std::mt19937 &rng);

// helper functions - could move to helper.hpp
// testing priority: 4 - these are mostly debugging functions to see what is happening 
void print_g_prod_info(g_prod_t g);
void print_g_struct_info(g_struct g);
void print_epsilon_info(epsilon_t eps);
void print_alpha_info(alpha_t alpha);
void print_complex_array(std::vector<std::complex<double>> &array);
void print_array(std::vector<std::vector<double>> &array);
void print_int_vector(std::vector<int> &array);

void print_pole_struct_info(pole_struct g);

void print_sign_array(sign_array_t signs);
void print_signs(sign_t signs);

void print_g_prod_array(g_prod_array_t g_array);
void print_pole_array(pole_array_t g);

void print_S(int dim, S_t &s_array);
void print_P( int dim, P_t &p_array);
void print_R( int dim, R_t &r_array);
void print_final( int dim, R_t &r_array, P_t &p_array, S_t &s_array);

void print_Pi( int dim, Pi_t &Pi_array);

void write_P(P_t &p_array);
void write_S(S_t &s_array);



// Functions for poles
//Testing Priority: 1 - this is an essential function - if it fails nothing will work 

pole_array_t find_poles(int index, g_prod_t &R);

// Residue Functions
// TODO: Simple_residue might be depricated - check usage
g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);

// Testing Priority: 1 - updating G with a given pole is an essential function - if it fails nothing will work 
g_struct update_G_pole(g_struct g_in,pole_struct pole);


// Functions for signs
// TODO: check usage of 'simple' for depricated 
sign_t find_signs(int index, g_prod_t &R);
double get_simple_sign(int index,g_prod_t &R,pole_struct pole);

// Functions for R's P's and S's
// TODO: Of these three only 'update_gprod_general' is actually used.  'Simple' is the case with no multi-poles.  
void update_gprod_simple(int index, R_t &R_array, P_t &P_array, S_t &S_array);

// Testing Priority: 1 - This is the central loop of the code. However it is a very complicated function so again test the internal functions and not the function itself 
void update_gprod_general(int int_index, int array_index,  R_t &R_array, P_t &P_array, S_t &S_array);

// This was an attempt to reorder the integration to try to minimize the number of resulting terms. It is not strictly necessary and its status is unknown. TODO: likely remove 
void update_gprod_general_minimal(int int_index, int array_index, R_t &R_array, P_t &P_array, S_t &S_array);

// Functions for Evaluation
// Testing Priority: 1 - these should all be testable - and form the backbone of the evaluation 
// This is the star function from AMI paper below equation 20: Define the function here and reference the paper: prb 99 035120
std::complex<double> star(ami_parms &parms, SorF_t H, Ri_t R, ami_vars external);
//TODO: This is depricated star function 
std::complex<double> star(SorF_t H, Ri_t R);

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


// Not sure what this function was used for. Most likely depricated 
energy_t random_energy(int N, std::mt19937 &rng);


// Multipole functions

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

// TODO: Check which of these is used or if both are used 
// Testing Priority: 2 should be easy 
void take_derivatives(AmiCalc::Ri_t &Wi, AmiCalc::pole_struct pole, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs);
void take_derivative_gprod(AmiCalc::g_prod_t &g_prod, AmiCalc::pole_struct pole, double start_sign, AmiCalc::Ri_t &r_out, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs);

// This is actually a pretty important function. probably needs a more clear name and documentation as to what it does 
g_struct der_fix(g_struct &g_in, double alpha);

double get_starting_sign(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole);
g_prod_t reduce_gprod(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole);


// derivatives of fermi functions
double frk(int r, int k);
int binomialCoeff(int n, int k);


// most below here is the AmiCalc class stuff 
public:

  ///Default Constructor
  AmiCalc();
  /// Constructor with ami_parms
  AmiCalc(ami_parms &parms);
  ///define parameter defaults
//  static void define_parameters(alps::params &p);  // this is needed for the evaluate stage
  ///the main calculation
//  void run(alps::params& parms);
  /// The construction

  void construct(ami_parms &parms);  // needs to know how many integrals to do
  
  /**
 * Construct AMI solution set .
 *
 * @param[out] mean the mean of `data`, or `NaN` if `data` is empty
 * @param[out] stdDev the unbiased (sample) standard deviation, or `NaN`
 *     if `data` contains fewer than 2 elements
 * @param[in] data the data to analyze
 */
 // Testing Priority: 1 - should be testable by picking a few specific inputs and knowing what the output should be 
 // construct is the MAIN construction function - so really important 
  void construct(ami_parms &parms,  g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array);
  // flip order of integration to try to reduce multipole issues
  // TODO: probably depricated and should be removed 
  void minimal_construct(ami_parms &parms, g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array);

  /// The evaluation - depricated
  void evaluate(ami_parms &parms);
  //void evaluate(alps::params &parms, R_t &R_array, P_t &P_array, S_t &S_array);


/** \brief Evaluation Function .
      * \param level an integer setting how useful to be
      * \return Output that is extra useful
      * 
      * This method does unbelievably useful things.  
      * And returns exceptionally useful results.
      * Use it everyday with good health.
      */
  std::complex<double> evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external);
  
// TODO: Check if depricated   
  void evaluate_multi_random(int NDAT, ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, std::mt19937 &rng);
  

private:




}; 




