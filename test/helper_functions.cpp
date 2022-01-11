

#include "../src/ami_base.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"
// include other headers if necessary   #include ../src/ami.hpps

TEST(helper_functions, energy_test){

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




// Checks if both alpha vectors are equal.
TEST(helper_functions, alpha_size){
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
TEST(helper_functions, epsilon_size){
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
TEST(helper_functions, alpha_vector_diff){
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
TEST(helper_functions, epsilon_vector_diff){
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
TEST(helper_functions, everything_works){
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



