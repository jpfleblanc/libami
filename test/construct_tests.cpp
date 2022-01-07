
#include "../src/ami_base.hpp"
#include "../src/ami_calc.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"


AmiBase::g_prod_t construct_example_6();


TEST(construct_tests, o6_ex_count){
	

AmiBase ami;

// Problem setup (see ami_example.cpp)
AmiBase::g_prod_t R0=construct_example_6(); // Sets initial integrand 

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects for S,P,R 
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;

AmiBase::terms amiterms;

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=6;  // Number of Matsubara sums to perform
AmiBase::ami_parms test_amiparms(N_INT, E_REG);

//Construction Stage
ami.construct(test_amiparms, R0, R_array, P_array, S_array); 
int first=R_array[N_INT].size();

// std::cout<< R_array.size()<< " "<<R_array[0].size()<< " "<< R_array[1].size()<< " "<<R_array[2].size()<< " "<<R_array[3].size()<< " "<<R_array[4].size()<< " "<<R_array[5].size()<< " "<<R_array[6].size()<<std::endl;

ami.construct(test_amiparms, R0, amiterms);
int second=amiterms.size();

// std::cout<<first<<" "<<second<<std::endl;


ASSERT_EQ(amiterms.size(),R_array.back().size());



}




AmiBase::g_prod_t construct_example_6(){

AmiBase::g_prod_t g;


// Setting up G array
// defining alpha's

AmiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
AmiBase::alpha_t alpha_2={0,1,0,0,0,0,0};
AmiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
AmiBase::alpha_t alpha_4={0,0,0,1,0,0,0};
AmiBase::alpha_t alpha_5={0,0,0,0,1,0,0};
AmiBase::alpha_t alpha_6={0,0,0,0,0,1,0};
AmiBase::alpha_t alpha_7={-1,1,0,0,0,0,1};
AmiBase::alpha_t alpha_8={-1,1,1,0,0,0,1};
AmiBase::alpha_t alpha_9={-1,1,1,1,0,0,1};
AmiBase::alpha_t alpha_10={-1,1,1,1,1,0,1};
AmiBase::alpha_t alpha_11={-1,1,1,1,1,1,1};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_4={0,0,0,1,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_5={0,0,0,0,1,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_6={0,0,0,0,0,1,0,0,0,0,0};
AmiBase::epsilon_t epsilon_7={0,0,0,0,0,0,1,0,0,0,0};
AmiBase::epsilon_t epsilon_8={0,0,0,0,0,0,0,1,0,0,0};
AmiBase::epsilon_t epsilon_9={0,0,0,0,0,0,0,0,1,0,0};
AmiBase::epsilon_t epsilon_10={0,0,0,0,0,0,0,0,0,1,0};
AmiBase::epsilon_t epsilon_11={0,0,0,0,0,0,0,0,0,0,1};


AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);
AmiBase::g_struct g4(epsilon_4,alpha_4);
AmiBase::g_struct g5(epsilon_5,alpha_5);
AmiBase::g_struct g6(epsilon_6,alpha_6);
AmiBase::g_struct g7(epsilon_7,alpha_7);
AmiBase::g_struct g8(epsilon_8,alpha_8);
AmiBase::g_struct g9(epsilon_9,alpha_9);
AmiBase::g_struct g10(epsilon_10,alpha_10);
AmiBase::g_struct g11(epsilon_11,alpha_11);


AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);
R0.push_back(g8);
R0.push_back(g9);
R0.push_back(g10);
R0.push_back(g11);





return R0;

}
