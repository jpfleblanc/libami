
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
// You are going to create two pole_struct (similar to g_struct) then,
// test both epsilons have the same size, and both have the same alpha
// size. Then, check each entry in eps and alpha if they are the same
// for both poles.


// Checks if both alpha vectors are equal.
TEST(pole_equiv, alpha_size){
	// Create two alpha objects
	AmiBase::alpha_t alpha_1={0,0,0,1,0};
	AmiBase::alpha_t alpha_2={0,0,0,1,0,0};
	
	// Create two epsilon objects
	AmiBase::epsilon_t epsilon_1={0,0,0,0,0,1,0};
	AmiBase::epsilon_t epsilon_2={0,0,0,0,0,1,0};
	
	// Create the pole_struct objects
	AmiBase::pole_struct pole_test_1;
	AmiBase::pole_struct pole_test_2;
	
	// Attach the eps and alpha to the pole_struct
	pole_test_1.alpha_ = alpha_1;
	pole_test_1.eps_ = epsilon_1;
	
	pole_test_2.alpha_ = alpha_2;
	pole_test_2.eps_ = epsilon_2;
	
	
	// Testing part
	AmiBase obj;
	ASSERT_FALSE(obj.pole_equiv(pole_test_1, pole_test_2)) << "Alpha size of poles are different" ;

// Test must tell the user where is the mistake (either the size
// not the same or which entry in the array was different)
}


// Checks if both epsilon vectors are equal
TEST(pole_equiv, epsilon_size){
	// Create two alpha objects
	AmiBase::alpha_t alpha_1={0,0,0,1,0};
	AmiBase::alpha_t alpha_2={0,0,0,1,0};
	
	// Create two epsilon objects
	AmiBase::epsilon_t epsilon_1={0,0,0,0,1,0};
	AmiBase::epsilon_t epsilon_2={0,0,0,0,0,1,0};
	
	// Create the pole_struct objects
	AmiBase::pole_struct pole_test_1;
	AmiBase::pole_struct pole_test_2;
	
	// Attach the eps and alpha to the pole_struct
	pole_test_1.alpha_ = alpha_1;
	pole_test_1.eps_ = epsilon_1;
	
	pole_test_2.alpha_ = alpha_2;
	pole_test_2.eps_ = epsilon_2;
	
	
	// Testing part
	AmiBase obj;
	ASSERT_FALSE(obj.pole_equiv(pole_test_1, pole_test_2)) << "Epsilon size of poles are different" ;
}


//Checks that elements in alpha vector are not the same################
TEST(pole_equiv, alpha_vector_diff){
	// Create two alpha objects
	AmiBase::alpha_t alpha_1={0,0,1,1,0};
	AmiBase::alpha_t alpha_2={0,0,0,1,0};
	
	// Create two epsilon objects
	AmiBase::epsilon_t epsilon_1={0,0,0,0,0,1,0};
	AmiBase::epsilon_t epsilon_2={0,0,0,0,0,1,0};
	
	// Create the pole_struct objects
	AmiBase::pole_struct pole_test_1;
	AmiBase::pole_struct pole_test_2;
	
	// Attach the eps and alpha to the pole_struct
	pole_test_1.alpha_ = alpha_1;
	pole_test_1.eps_ = epsilon_1;
	
	pole_test_2.alpha_ = alpha_2;
	pole_test_2.eps_ = epsilon_2;
	
	
	// Testing part
	AmiBase obj;
	ASSERT_FALSE(obj.pole_equiv(pole_test_1, pole_test_2)) << "Vector alpha entries are different for both poles." ;
	
}	


//Checks that elements in epsilon vector are not the same
TEST(pole_equiv, epsilon_vector_diff){
	// Create two alpha objects
	AmiBase::alpha_t alpha_1={0,0,0,1,0};
	AmiBase::alpha_t alpha_2={0,0,0,1,0};
	
	// Create two epsilon objects
	AmiBase::epsilon_t epsilon_1={0,0,2,0,0,1,0};
	AmiBase::epsilon_t epsilon_2={0,0,0,0,0,1,0};
	
	// Create the pole_struct objects
	AmiBase::pole_struct pole_test_1;
	AmiBase::pole_struct pole_test_2;
	
	// Attach the eps and alpha to the pole_struct
	pole_test_1.alpha_ = alpha_1;
	pole_test_1.eps_ = epsilon_1;
	
	pole_test_2.alpha_ = alpha_2;
	pole_test_2.eps_ = epsilon_2;
	
	
	// Testing part
	AmiBase obj;
	ASSERT_FALSE(obj.pole_equiv(pole_test_1, pole_test_2)) << "Vector epsilon entries are different." ;
	
}	


//Checks that the function will spit TRUE when poles are identical.
TEST(pole_equiv, everything_works){
	// Create two alpha objects
	AmiBase::alpha_t alpha_1={0,0,0,1,0};
	AmiBase::alpha_t alpha_2={0,0,0,1,0};
	
	// Create two epsilon objects
	AmiBase::epsilon_t epsilon_1={0,0,0,0,0,1,0};
	AmiBase::epsilon_t epsilon_2={0,0,0,0,0,1,0};
	
	// Create the pole_struct objects
	AmiBase::pole_struct pole_test_1;
	AmiBase::pole_struct pole_test_2;
	
	// Attach the eps and alpha to the pole_struct
	pole_test_1.alpha_ = alpha_1;
	pole_test_1.eps_ = epsilon_1;
	
	pole_test_2.alpha_ = alpha_2;
	pole_test_2.eps_ = epsilon_2;
	
	
	// Testing part
	AmiBase obj;
	ASSERT_TRUE(obj.pole_equiv(pole_test_1, pole_test_2)) << "Alpha and Epsilon vectors are identical for both poles." ;
	
}	


