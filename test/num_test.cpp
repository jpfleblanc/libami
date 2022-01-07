
#include "../src/ami_base.hpp"
#include "../src/ami_calc.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"

AmiBase::ami_vars construct_4ord_ext_multipole_example();
AmiBase::g_prod_t construct_multipole_example();

bool example4(std::complex<double> &value);


TEST(num_test, o4_test){
	
std::complex<double> from_benchmark;
ASSERT_TRUE(example4(from_benchmark));	

std::complex<double> from_func(-0.0009888959048396373,0.00075325682931293159);

// std::cout<<from_benchmark<<" "<<from_func<<std::endl;
	
ASSERT_DOUBLE_EQ(from_func.real(),from_benchmark.real());
ASSERT_DOUBLE_EQ(from_func.imag(),from_benchmark.imag());
}





bool example4(std::complex<double> &value){

//START example	
// class instance 
AmiBase ami;

// Problem setup
AmiBase::g_prod_t R0=construct_multipole_example();
AmiBase::ami_vars avars=construct_4ord_ext_multipole_example();

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects 
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;

// Integration/Evaluation parameters
double E_REG=0;
int N_INT=4;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);

//Construction Stage
ami.construct(test_amiparms, R0, R_array, P_array, S_array);

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	
	
	//timing info
	auto t_end=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::duration<double> diff2=t_end-t2;

	std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

//END Example


// START Same Problem with Optimization 
	auto t3=std::chrono::high_resolution_clock::now();

// Storage Structures
AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

// Take existing solution and factorize it 
ami.factorize_Rn(R_array.back(), unique, rref, eval_list);

	//timing info 
	auto t4=std::chrono::high_resolution_clock::now();	

// Evaluate optimized solution - unique, rref and eval_list are not empty.	
std::complex<double> opt_calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars, unique, rref, eval_list);	

	//timing info
	auto t5=std::chrono::high_resolution_clock::now();	

	std::chrono::duration<double> diff3=t4-t3;
	std::chrono::duration<double> diff4=t5-t4;

	std::chrono::microseconds d3=std::chrono::duration_cast<std::chrono::microseconds>(diff3);
	std::chrono::microseconds d4=std::chrono::duration_cast<std::chrono::microseconds>(diff4);	

//END Optimization example

//START Term Example
// Same Problem, term storage type

	//timing info
	auto t6=std::chrono::high_resolution_clock::now();

//simplified storage type 
AmiBase::terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

	//timing info 
	auto t7=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff5=t7-t6;
	std::chrono::microseconds d5=std::chrono::duration_cast<std::chrono::microseconds>(diff5);

	//timing info 
	auto t8=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
//END Term Example 	

// Start Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list);
	auto t11=std::chrono::high_resolution_clock::now();


	auto t12=std::chrono::high_resolution_clock::now();
std::complex<double> opt_term_val=ami.evaluate(test_amiparms, amiterms, avars,unique, rref, eval_list);
	auto t13=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff7=t11-t10;
	std::chrono::duration<double> diff8=t13-t12;

	std::chrono::microseconds d7=std::chrono::duration_cast<std::chrono::microseconds>(diff7);
	std::chrono::microseconds d8=std::chrono::duration_cast<std::chrono::microseconds>(diff8);	

value=calc_result;

//std::cout<<term_val<<" "<<calc_result<<std::endl;
if(std::abs(term_val-calc_result)<1e-8){return true;}else{return false;}
	
}



AmiBase::g_prod_t construct_multipole_example(){

AmiBase::g_prod_t g;


// Setting up G array
// defining alpha's
AmiBase::alpha_t alpha_1={1,0,0,1,-1};
AmiBase::alpha_t alpha_2={0,0,0,1,0};
AmiBase::alpha_t alpha_3={1,0,0,0,0};
AmiBase::alpha_t alpha_4={1,0,0,0,0};
AmiBase::alpha_t alpha_5={0,1,0,0,0};
AmiBase::alpha_t alpha_6={-1,1,1,0,0};
AmiBase::alpha_t alpha_7={0,0,1,0,0};


//defining epsilon's
AmiBase::epsilon_t epsilon_1={0,0,0,0,0,1,0};
AmiBase::epsilon_t epsilon_2={0,0,0,0,1,0,0};
AmiBase::epsilon_t epsilon_3={1,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_4={1,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_5={0,1,0,0,0,0,0};
AmiBase::epsilon_t epsilon_6={0,0,0,1,0,0,0};
AmiBase::epsilon_t epsilon_7={0,0,1,0,0,0,0};


AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);
AmiBase::g_struct g4(epsilon_4,alpha_4);
AmiBase::g_struct g5(epsilon_5,alpha_5);
AmiBase::g_struct g6(epsilon_6,alpha_6);
AmiBase::g_struct g7(epsilon_7,alpha_7);


AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);



return R0;
}



AmiBase::ami_vars construct_4ord_ext_multipole_example(){



AmiBase::energy_t energy={1,1.1,1.2,1.31,1.4,0.01, 0.1}; //{1,1.1,1.2,1.3,1.4,0.01, 0.1};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};


AmiBase::ami_vars external(energy, frequency);

external.BETA_=1.0;

return external;

}