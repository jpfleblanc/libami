/*
 * Copyright (C) 2018 JPF LeBlanc jleblanc@mun.ca  See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

#pragma once

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>

int factorial(int n);


class AmiCalc 
{

public:

/// Type definitions for AMI

//Fundamental objects

 
typedef std::vector<int> epsilon_t;  /// the symbolic epsilon
typedef std::vector<int> alpha_t;
typedef enum {Bose,Fermi} stat_type ;
typedef enum {Sigma,Pi,Hartree} graph_type ;

// AMI Parameter structure
// This is a little weird.  The ami_parms structure is defined here. But there is no internal variable set for this... meaning it has to be passed to all of the functions... probably this should be changed. Explicit passing is 'ok'.

struct ami_parms{
ami_parms(int N_INT,  double E_REG){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=static_cast<AmiCalc::graph_type>(0); /// by default sigma
}

ami_parms(int N_INT, double E_REG, graph_type TYPE){
N_INT_=N_INT;
E_REG_=E_REG;
TYPE_=TYPE;
}

ami_parms(){}

int N_INT_;
double E_REG_;
graph_type TYPE_;

};


// Green's function structure which has a symbolic epsilon and set of alpha values
struct g_struct{
g_struct(epsilon_t eps, alpha_t alpha, stat_type stat){
eps_=eps;
alpha_=alpha;
stat_=stat;

}

// Assume fermi statistics if not specified for partially initialized structure
g_struct(epsilon_t eps, alpha_t alpha){
eps_=eps;
alpha_=alpha;
stat_=Fermi;

}

g_struct(){} // uninitialized variant

epsilon_t eps_;
alpha_t alpha_;
stat_type stat_;

};

// Pole structure. Equivalent to G structure, but kept separate. Tracks multiplicity, and which green's function it is attached to.
struct pole_struct{


epsilon_t eps_;
alpha_t alpha_;
int index_;
int multiplicity_=1;
std::vector<int> which_g_;

};


// External Variable structure and typedefs

/// external list of energies used with epsilon_t to construct energies

typedef std::vector<std::complex<double>> energy_t;
typedef std::vector< std::complex<double>  > frequency_t;

typedef std::vector< double> k_vector_t;
typedef std::vector< k_vector_t> k_vect_list_t;





typedef std::vector<double> sign_t;
typedef std::vector<g_struct> g_prod_t;
typedef std::vector<pole_struct> pole_array_t;

typedef std::vector<sign_t> sign_array_t;
typedef std::vector<g_prod_t> g_prod_array_t;



// arrays of arrays OMG this is annoying

typedef std::vector<g_prod_t> Ri_t;
typedef std::vector<pole_array_t> Pi_t;
typedef std::vector<sign_t> Si_t;

// final output arrays

typedef std::vector<g_prod_array_t> R_t;  // TODO: is there a reason this is not a vector of Ri_t ???
typedef std::vector<Pi_t> P_t;
typedef std::vector<Si_t> S_t;

// typedefs for evaluation

typedef std::vector< std::vector<std::complex<double> > > SorF_t;

//-----------
//std::vector<std::vector<int>> intvec;
//intvec filled somehow

//std::vector<std::vector<double>> doublevec;
//doublevec.reserve(intvec.size());
//for (auto&& v : intvec) doublevec.emplace_back(std::begin(v), std::end(v));
//-----------

// various structures

// internal state - list of internal k'sb_type
// ext_vars and external_variable_list : loaded from a file and populated with external parameters to be evaluated
// ami_vars - actual input to ami routines for constructing solutions

struct internal_state{

internal_state(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
}

internal_state(){}

void initialize(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
	
}

k_vect_list_t internal_k_list_;
int order_;
int dim_;
double prefactor_;
// R_t R_array_;
// P_t P_array_;
// S_t S_array_;

};


struct ext_vars{

ext_vars(int dim){
	
KDIM_=dim;
external_k_vector_.assign(dim,0.0);	
	
}
ext_vars(int dim, double beta, std::complex<double> mu){
	
KDIM_=dim;
external_k_vector_.assign(dim,0.0);	

BETA_=beta;
MU_=mu;
	
}

ext_vars(){}


k_vector_t external_k_vector_;	
int KDIM_;
frequency_t external_freq_;
double BETA_;
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






// IO function for external variables
void read_external(std::string filename, external_variable_list &extern_list);
void read_text_S_solutions(std::string filename, S_t &s_array);
void read_text_P_solutions(std::string eps_filename,std::string alpha_filename, P_t &p_array);
void read_text_R_solutions(std::string eps_filename,std::string alpha_filename, R_t &r_array, int size);
double load_prefactor(std::string filename, int order);
void write_S_readable(S_t &s_array);
void write_P_readable(P_t &p_array);
void write_R_readable(R_t &r_array);


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

ami_vars construct_random_example_J(std::mt19937 &rng);

// helper functions - could move to helper.hpp
void print_g_prod_info(g_prod_t g);
void print_g_struct_info(g_struct g);
void print_epsilon_info(epsilon_t eps);
void print_alpha_info(alpha_t alpha);

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

pole_array_t find_poles(int index, g_prod_t &R);

// Residue Functions

g_prod_t simple_residue(g_prod_t G_in, pole_struct pole);
g_struct update_G_pole(g_struct g_in,pole_struct pole);


// Functions for signs
sign_t find_signs(int index, g_prod_t &R);

// Functions for R's P's and S's

void update_gprod_simple(int index, R_t &R_array, P_t &P_array, S_t &S_array);
void update_gprod_general(int index, R_t &R_array, P_t &P_array, S_t &S_array);

// Functions for Evaluation

std::complex<double> star(ami_parms &parms, SorF_t H, Ri_t R, ami_vars external);
std::complex<double> star(SorF_t H, Ri_t R);
SorF_t fermi(ami_parms &parms, Pi_t pi, ami_vars external);
std::complex<double>  fermi_pole(ami_parms &parms, pole_struct pole, ami_vars external);
SorF_t cross(SorF_t left, SorF_t right);
SorF_t dot(Si_t Si, SorF_t fermi);

std::complex<double> get_energy_from_pole(pole_struct pole, ami_vars external);
std::complex<double>  get_energy_from_g(g_struct g, ami_vars external);
//std::complex<double> eval_gprod(g_prod_t g_prod, ami_vars external);
std::complex<double> eval_gprod(ami_parms &parms, g_prod_t g_prod, ami_vars external);



energy_t random_energy(int N, std::mt19937 &rng);


// Multipole functions
// TODO Multipole construction not working

bool pole_equiv(pole_struct pole1, pole_struct pole2);

void evaluate_general_residue(g_prod_t G_in, pole_struct pole, Ri_t &Ri, pole_array_t &poles, sign_t &signs);
void take_derivatives(AmiCalc::Ri_t &Wi, AmiCalc::pole_struct pole, AmiCalc::pole_array_t &poles, AmiCalc::sign_t &signs);
double get_starting_sign(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole);
g_prod_t reduce_gprod(AmiCalc::g_prod_t G_in, AmiCalc::pole_struct pole);


public:

  ///setup of parameters
  AmiCalc();
  AmiCalc(ami_parms &parms);
  ///define parameter defaults
//  static void define_parameters(alps::params &p);  // this is needed for the evaluate stage
  ///the main calculation
//  void run(alps::params& parms);
  /// The construction

  void construct(ami_parms &parms);  // needs to know how many integrals to do
  void construct(ami_parms &parms,  g_prod_t R0, R_t &R_array, P_t &P_array, S_t &S_array);

  /// The evaluation 
  void evaluate(ami_parms &parms);
  //void evaluate(alps::params &parms, R_t &R_array, P_t &P_array, S_t &S_array);

  std::complex<double> evaluate(ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, ami_vars &external);
  std::complex<double> evaluate_multi_random(int NDAT, ami_parms &parms, R_t &R_array, P_t &P_array, S_t &S_array, std::mt19937 &rng);
  

private:




}; 




