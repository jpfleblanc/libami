
#include "../src/ami_base.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"


AmiBase::g_prod_t construct_example_6();
void print_g_struct_info(AmiBase::g_struct g);
void print_g_prod_info(AmiBase::g_prod_t g);
void print_epsilon_info(AmiBase::epsilon_t eps);
void print_alpha_info(AmiBase::alpha_t alpha);
void print_pole_struct_info(AmiBase::pole_struct g);


TEST(construct_tests, term_deriv){
	
AmiBase ami;

//take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

AmiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
AmiBase::alpha_t alpha_2={0,1,0,0,0,0,0};
AmiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
	
AmiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);




AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);


AmiBase::pole_array_t poles=ami.find_poles(0,R0);
AmiBase::pole_struct pole=poles[0];

AmiBase::Ri_t ri_out;
AmiBase::pole_array_t pi_out;
AmiBase::sign_t si_out;
ami.take_derivative_gprod(R0,pole,1.0,ri_out, pi_out,si_out);
	

ASSERT_EQ(ri_out.size(),2);
ASSERT_EQ(pi_out.size(),2);
ASSERT_EQ(si_out.size(),2);
	
ASSERT_EQ(ri_out[0].size(),3);
ASSERT_EQ(ri_out[1].size(),4);	
	
// print_pole_struct_info(pi_out[0]);
ASSERT_EQ(pi_out[0].der_,1);
ASSERT_EQ(pi_out[1].der_,0);
ASSERT_EQ(ami.pole_equiv(pi_out[0],pi_out[1]),true);	
	
}


TEST(construct_tests, term_deriv2){
	
AmiBase ami;

//take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

AmiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
AmiBase::alpha_t alpha_2={1,1,0,0,0,0,0};
AmiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
	
AmiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);




AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);


AmiBase::pole_array_t poles=ami.find_poles(0,R0);
AmiBase::pole_struct pole=poles[0];

AmiBase::Ri_t ri_out;
AmiBase::pole_array_t pi_out;
AmiBase::sign_t si_out;
ami.take_derivative_gprod(R0,pole,1.0,ri_out, pi_out,si_out);
	
ASSERT_EQ(ri_out.size(),3);
ASSERT_EQ(pi_out.size(),3);
ASSERT_EQ(si_out.size(),3);
	
// print_g_prod_info(ri_out[0]);
// std::cout<<"--"<<std::endl;	
// print_g_prod_info(ri_out[1]);
// std::cout<<"--"<<std::endl;
// print_g_prod_info(ri_out[2]);
// std::cout<<"--"<<std::endl;	

ASSERT_EQ(pi_out[0].der_,1);
ASSERT_EQ(pi_out[1].der_,0);
ASSERT_EQ(pi_out[2].der_,0);

ASSERT_EQ(ami.g_equiv(ri_out[1][0],ri_out[1][1]),true);
ASSERT_EQ(ami.g_equiv(ri_out[2][1],ri_out[2][2]),true);
	
}



TEST(construct_tests, update_g){
	
AmiBase ami;

AmiBase::alpha_t alpha_1={-1,1,0,0,0,0,1};	
AmiBase::epsilon_t eps_1={0,1,0,0,1,0,-1,0,0,0,0};

AmiBase::g_struct g1(eps_1,alpha_1);


AmiBase::g_prod_t R0=construct_example_6();
AmiBase::pole_array_t pole_list=ami.find_poles(0,R0);

AmiBase::g_struct updated=ami.update_G_pole(g1, pole_list[0]);

// print_g_struct_info(g1);

// print_g_struct_info(updated);

AmiBase::epsilon_t eps_result={1, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0};
AmiBase::alpha_t alpha_result={0, 1, 0, 0, 0, 0, 1};

AmiBase::pole_struct p1(eps_result,alpha_result);
AmiBase::pole_struct p2(updated.eps_,updated.alpha_);

ami.pole_equiv(p1, p2);
	
	
}

TEST(construct_tests, find_poles){
	

AmiBase ami;

// Problem setup (see ami_example.cpp)
AmiBase::g_prod_t R0=construct_example_6();

AmiBase::pole_array_t pole_list_0=ami.find_poles(0,R0);
AmiBase::pole_array_t pole_list_1=ami.find_poles(1,R0);
AmiBase::pole_array_t pole_list_2=ami.find_poles(2,R0);
AmiBase::pole_array_t pole_list_3=ami.find_poles(3,R0);
AmiBase::pole_array_t pole_list_6=ami.find_poles(6,R0);

std::cout<<"Expect warnings: { ";
AmiBase::pole_array_t pole_list_7=ami.find_poles(7,R0);
AmiBase::pole_array_t pole_list_12=ami.find_poles(12,R0);
std::cout<<"} Warnings caught"<<std::endl;
	
ASSERT_EQ(pole_list_0.size(),6);
ASSERT_EQ(pole_list_1.size(),6);
ASSERT_EQ(pole_list_2.size(),5);
ASSERT_EQ(pole_list_3.size(),4);
ASSERT_EQ(pole_list_6.size(),5);
ASSERT_EQ(pole_list_7.size(),0);
ASSERT_EQ(pole_list_12.size(),0);	
	
}


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



void print_g_struct_info(AmiBase::g_struct g){

std::cout<<"Species="<<g.species_<<" ";
std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}



void print_epsilon_info(AmiBase::epsilon_t eps){


//for (std::vector<signed char>::iterator it= eps.begin(); it != eps.end(); ++it){
for (std::vector<int>::iterator it= eps.begin(); it != eps.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}


void print_alpha_info(AmiBase::alpha_t alpha){

for (std::vector<int>::iterator it= alpha.begin(); it != alpha.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}

void print_pole_struct_info(AmiBase::pole_struct g){

std::cout<<"Der="<<g.der_<<" | ";
std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}


void print_g_prod_info(AmiBase::g_prod_t g){

std::cout<<"----------Printing G_product---------- " <<std::endl;
for (std::vector<AmiBase::g_struct>::iterator it= g.begin(); it != g.end(); ++it){

std::cout<<"----------Printing next---------- " <<std::endl;
print_g_struct_info(it[0]);

}


}