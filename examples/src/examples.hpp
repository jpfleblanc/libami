#include "ami_base.hpp"


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
void example1_bose();
void example4();
void example4hm();
void example9();

void example_debug();
AmiBase::g_prod_t construct_example_debug();
AmiBase::ami_vars construct_ext_example_debug();

void example2_spec();
AmiBase::ami_vars construct_spec_example();

/// Various functions

AmiBase::g_prod_t construct_multipole_example();
AmiBase::g_prod_t construct_hm_multipole_example();
AmiBase::g_prod_t construct_example2();
AmiBase::g_prod_t construct_example1_bose();

AmiBase::g_prod_t construct_example_Y();
AmiBase::g_prod_t construct_example_J();


// These are depricated
AmiBase::ami_vars construct_ext_example_Y();
AmiBase::ami_vars construct_ext_example_J();
AmiBase::ami_vars construct_ext_example2();
AmiBase::ami_vars construct_ext_example1_bose();
AmiBase::ami_vars construct_ext_multipole_example();
AmiBase::ami_vars construct_6ord_ext_multipole_example();
AmiBase::ami_vars construct_4ord_ext_multipole_example();

AmiBase::ami_vars construct_random_example_J(std::mt19937 &rng);
AmiBase::energy_t random_energy(int N, std::mt19937 &rng);


