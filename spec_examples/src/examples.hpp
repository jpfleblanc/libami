#include "ami_base.hpp"
#include "ami_calc.hpp"

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>

#include <chrono>
#include <thread>

void example2();
void example4();
void example4hm();
void example9();

void example2_spec();
void example_simple_spec();
void example_simple_nospec();
AmiBase::g_prod_t construct_simple_example();
AmiBase::ami_vars ext_example_simple(double w, double gamma);

AmiBase::ami_vars construct_spec_example();

/// Various functions

AmiBase::g_prod_t construct_multipole_example();
AmiBase::g_prod_t construct_hm_multipole_example();
AmiBase::g_prod_t construct_example();
AmiBase::g_prod_t construct_example_Y();
AmiBase::g_prod_t construct_example_J();


// These are depricated
AmiBase::ami_vars construct_ext_example_Y();
AmiBase::ami_vars construct_ext_example_J();
AmiBase::ami_vars construct_ext_example_Sigma();
AmiBase::ami_vars construct_ext_multipole_example();
AmiBase::ami_vars construct_6ord_ext_multipole_example();
AmiBase::ami_vars construct_4ord_ext_multipole_example();

AmiBase::ami_vars construct_random_example_J(std::mt19937 &rng);
AmiBase::energy_t random_energy(int N, std::mt19937 &rng);


