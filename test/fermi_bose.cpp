

#include "../src/ami_base.hpp"
#include "../src/ami_calc.hpp"
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


EXPECT_EQ(from_func,from_analytic);

}	
