

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

TEST(get_energy, energy_test){

AmiBase obj;

std::complex<double> from_func, from_analytic;

AmiBase::pole_struct pole;
AmiBase::epsilon_t epsilon_1={1,1,1,1,1};
pole.eps_=epsilon_1;


AmiBase::ami_vars external;
AmiBase::energy_t epspole={1,1.01,0,1.02,1.03};

external.energy_=epspole;

from_func=obj.get_energy_from_pole(pole, external);

from_analytic=4.06;


EXPECT_EQ(from_func,from_analytic);

}	
