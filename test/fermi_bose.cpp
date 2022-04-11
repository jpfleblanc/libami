

#include "../src/ami_base.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"
// include other headers if necessary   #include ../src/ami.hpps

TEST(fermi_bose, fermi_test){

//(int m, double sigma, double beta, std::complex<double> E)
// m: mth derivative
// sigma: +1 for fermion, -1 for boson
// beta: inverse template
// E: complex energy 

AmiBase obj;

std::complex<double> from_func, from_analytic;

std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(0, 1.0,  beta, E);

from_analytic=1.0/(std::exp(beta*E)+1.0);


ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());

}	

TEST(fermi_bose, fermi_test_deriv1){
	
AmiBase obj;

// std::complex<double> from_func, from_analytic;
std::complex<double> from_func(0,0);
std::complex<double> from_analytic(0,0);
std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(1, 1.0,  beta, E);

from_analytic=-( beta*std::exp(E*beta) / ( std::pow( (std::exp(E*beta) +1.0), 2.0)) );

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());	
	
}

TEST(fermi_bose, fermi_test_deriv2_james){
	
AmiBase obj;

// std::complex<double> from_func, from_analytic;
std::complex<double> from_func(0,0);
std::complex<double> from_analytic(0,0);
std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(2, 1.0,  beta, E);

from_analytic=2.0*(std::pow(beta,2)*std::exp(2.0*beta*E))/ ( std::pow( (std::exp(E*beta) +1.0), 3.0)) - std::pow(beta,2)*std::exp(beta*E)/ ( std::pow( (std::exp(E*beta) +1.0), 2.0));

ASSERT_NEAR(from_func.real(),from_analytic.real(),1e-13);	
	
}


//testing second derivative of fermi function
TEST(fermi_bose, fermi_test_deriv2){
 

AmiBase obj;

// std::complex<double> from_func, from_analytic;
std::complex<double> from_func(0,0);
std::complex<double> from_analytic(0,0);
std::complex<double> from_analytic2(0,0);
std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(2, 1.0,  beta, E);

from_analytic=( std::pow( beta, 2.0) * (std::exp(E*beta)-1.0)*std::exp(E*beta) / ( std::pow( (std::exp(E*beta) +1.0), 3.0)) );

// from_analytic2=2.0*(std::pow(beta,2)*std::exp(2.0*beta*E))/ ( std::pow( (std::exp(E*beta) +1.0), 3.0)) - std::pow(beta,2)*std::exp(beta*E)/ ( std::pow( (std::exp(E*beta) +1.0), 2.0));

// ASSERT_DOUBLE_EQ(from_analytic.real(),from_analytic2.real());

ASSERT_NEAR(from_func.real(),from_analytic.real(), 1e-13);

}

// testing fifth derivative of fermi function
TEST(fermi_bose, fermi_test_deriv5){
 

AmiBase obj;

std::complex<double> from_func, from_analytic;

std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(5, 1.0,  beta, E);

from_analytic= std::pow( beta, 5.0) * std::exp(E*beta) * (-std::exp(4.0*beta*E)+26.0*std::exp(3.0*E*beta) - 66.0*std::exp(2.0*E*beta) + 26.0*std::exp(E*beta) - 1.0 ) / (std::pow((std::exp(E*beta) +1.0),6.0)) ;

ASSERT_NEAR(from_func.real(),from_analytic.real(),1e-12);

}






TEST(fermi_bose, bose_test){


AmiBase obj;

std::complex<double> from_func, from_analytic;

std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(0, -1.0,  beta, E);

from_analytic=1.0/(-std::exp(beta*E)+1.0);

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());

}	

TEST(fermi_bose, bose_test_deriv1){


AmiBase obj;

std::complex<double> from_func, from_analytic;

std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(1, -1.0,  beta, E);

from_analytic=beta*std::exp(beta*E)/std::pow(std::exp(beta*E)-1.0,2);

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());

}	


//second bose derivative
TEST(fermi_bose, bose_test_deriv2){
 

AmiBase obj;

// std::complex<double> from_func, from_analytic;
std::complex<double> from_func(0,0);
std::complex<double> from_analytic(0,0);
std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(2, -1.0,  beta, E);

from_analytic=-( std::pow( beta, 2.0) * (std::exp(E*beta)+1.0)*std::exp(E*beta) / ( std::pow( (std::exp(E*beta)-1.0), 3.0)) );

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());

}