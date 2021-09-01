
#pragma once
#include "ami_base.hpp"
#include "ami_calc.hpp"

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
NewAmiCalc ami;

double xi_cutoff=10;
std::vector<double> AMI_spec_se_Im_vector;
std::vector<double> AMI_spec_se_Im_err_vector;
std::vector<double> AMI_spec_se_Re_vector;
std::vector<double> AMI_spec_se_Re_err_vector;
std::vector<double> AMI_spec_freq_vector;
std::vector<double> AMI_spec_ky_vector;
std::vector<double> AMI_spec_kx_vector;
std::vector<double> AMI_spec_freq_vector_simplified;
std::vector<double> AMI_spec_kx_vector_simplified;
std::vector<double> AMI_spec_ky_vector_simplified;
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


typedef AmiBase::g_struct delta_t;
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

bool root=true;
int delta_count=0;

};

typedef std::vector<ami_sp_term> sp_terms;




std::complex<double> A_eval(std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E);

std::complex<double> get_A(A_struct &A, double this_x, NewAmiCalc::k_vector_t k);
std::complex<double> eval_Aprod(A_prod_t &Ap, xi_t &xi, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu);
std::complex<double> eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu);


std::complex<double> get_X(X_t &Xsym, xi_t &xi);
std::complex<double> get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps);
std::complex<double> get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X);
void find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec);
std::vector<double> return_simple_grid_vector(std::vector<double> &in_vector);
void read_self_energy(std::string file_name);

std::complex<double> construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu);

void generate_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0);
void R0_to_Aprod(AmiBase::g_prod_t &R0, A_prod_t &Ap);

// void reduce_deltas(ami_sp_term &term);
void resolve_deltas( sp_terms &sp_terms);
void resolve_deltas(ami_sp_term &sp_term);

void replace_xi(int i, AmiBase::pole_array_t &pv, ami_sp_term &sp_term);
void update_spec_pole(AmiBase::pole_struct &source_pole, AmiBase::alpha_t &target_alpha, AmiBase::epsilon_t &target_eps);


AmiSpec(double xc);

// io functions for debugging
void print_delta_prod_t(delta_prod_t &delta_prod);
void print_delta_t(delta_t &delta);
void print_a_prod_t(A_prod_t &Aprod);
void print_a_struct(A_struct &A);
void print_sp_term(ami_sp_term &term);
void print_sp_terms(sp_terms &sp_terms);
void print_int_vec(std::vector<int> vec);



private:


};
