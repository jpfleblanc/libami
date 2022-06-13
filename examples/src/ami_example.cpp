//=======================================================================
// Copyright 2018 James PF LeBlanc
//=======================================================================


#include "ami_base.hpp"
#include "examples.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>



AmiBase::g_prod_t construct_example2(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


AmiBase::alpha_t alpha_1={1,0,0};
AmiBase::alpha_t alpha_2={0,1,0};
AmiBase::alpha_t alpha_3={-1,1,1};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0};
AmiBase::epsilon_t epsilon_3={0,0,1};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);

// AmiBase::g_prod_t R0={g1,g2,g3};

// OR
AmiBase::g_prod_t R0;
R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);




return R0;

}

AmiBase::ami_vars construct_ext_example2(){


AmiBase::energy_t energy={-4,0.1,-1};

AmiBase::frequency_t frequency;

for(int i=0;i<2;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}


AmiBase::g_prod_t construct_example1_bose(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


AmiBase::alpha_t alpha_1={1,0,0};
AmiBase::alpha_t alpha_2={1,1,0};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);


AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);


return R0;

}

AmiBase::ami_vars construct_ext_example1_bose(){


AmiBase::energy_t energy={-4,0.1};

AmiBase::frequency_t frequency;

for(int i=0;i<1;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));// This frequency is expected to be a bosonic matsubara or real frequency. There is no catch if this is untrue. 

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);



return external;

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



AmiBase::ami_vars construct_ext_multipole_example(){



AmiBase::energy_t energy={1,1.01,0,1.02,1.03};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}

AmiBase::ami_vars construct_6ord_ext_multipole_example(){



AmiBase::energy_t energy={1,1.1,1.2,1.3,1.4,0, 0.1, 0.2, 0.3,0.4, 0.5};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}



AmiBase::ami_vars construct_4ord_ext_multipole_example(){



AmiBase::energy_t energy={1,1.1,1.2,1.31,1.4,0.01, 0.1}; //{1,1.1,1.2,1.3,1.4,0.01, 0.1};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}




AmiBase::ami_vars construct_spec_example(){



// energy_t energy={-4,1,-1};
AmiBase::energy_t energy={-0.2, 0.35, -0.44 };

AmiBase::frequency_t frequency;

for(int i=0;i<2;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}





AmiBase::ami_vars construct_ext_example_Y(){



AmiBase::energy_t energy={4,-1,-1,-1,-2,2,1};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}


AmiBase::g_prod_t construct_example_Y(){

AmiBase::g_prod_t g;


// Setting up G array
// defining alpha's

AmiBase::alpha_t alpha_1={1,0,0,0,0};
AmiBase::alpha_t alpha_2={0,1,0,0,0};
AmiBase::alpha_t alpha_3={0,0,1,0,0};
AmiBase::alpha_t alpha_4={0,0,0,1,0};
AmiBase::alpha_t alpha_5={-1,1,0,0,1};
AmiBase::alpha_t alpha_6={0,1,1,-1,0};
AmiBase::alpha_t alpha_7={1,0,1,0,-1};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0};
AmiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0};
AmiBase::epsilon_t epsilon_4={0,0,0,1,0,0,0};
AmiBase::epsilon_t epsilon_5={0,0,0,0,1,0,0};
AmiBase::epsilon_t epsilon_6={0,0,0,0,0,1,0};
AmiBase::epsilon_t epsilon_7={0,0,0,0,0,0,1};


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


AmiBase::ami_vars construct_ext_example_J(){

AmiBase::energy_t energy={-4.64,1.02,1.04,1.05,1.06,1.07,1.08,1.09,1.11,1.23,-4.43,-4.52,1.5,1.6,1.7,1.8,1.9};
		//{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//energy_t energy={-4,1,-1,1,1,-4,1,1,1,1,1,1,1,1,1,1,1};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

//for(int i=0;i<9;i++){ frequency.push_back(std::complex<double>(0,0));}//frequency.push_back(std::complex<double>(0,0));}

//frequency.push_back(std::complex<double>(0, M_PI));//(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}


AmiBase::g_prod_t construct_example_J(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's


AmiBase::alpha_t alpha_1={1,0,0,0,0,0,0,0,0,0};
AmiBase::alpha_t alpha_2={0,1,0,0,0,0,0,0,0,0};
AmiBase::alpha_t alpha_3={1,1,0,0,0,0,0,0,0,-1};
AmiBase::alpha_t alpha_4={0,0,1,0,0,0,0,0,0,0};
AmiBase::alpha_t alpha_5={0,0,0,1,0,0,0,0,0,0};
AmiBase::alpha_t alpha_6={1,1,1,-1,0,0,0,0,0,-1};
AmiBase::alpha_t alpha_7={0,0,0,0,1,0,0,0,0,0};
AmiBase::alpha_t alpha_8={0,0,0,0,0,0,0,1,0,0};
AmiBase::alpha_t alpha_9={0,-1,0,0,1,0,0,1,0,0};
AmiBase::alpha_t alpha_10={0,0,1,-1,1,0,0,0,0,0};
AmiBase::alpha_t alpha_11={0,0,0,0,0,1,0,0,0,0};
AmiBase::alpha_t alpha_12={0,0,0,0,0,0,0,0,1,0};
AmiBase::alpha_t alpha_13={0,0,0,0,0,0,1,0,0,0};
AmiBase::alpha_t alpha_14={0,0,-1,1,-1,1,1,0,0,0};
AmiBase::alpha_t alpha_15={0,-1,-1,1,0,0,0,1,0,1};
AmiBase::alpha_t alpha_16={0,0,1,-1,1,-1,0,0,1,0};
AmiBase::alpha_t alpha_17={0,-1,0,0,1,-1,0,1,1,0};


//defining epsilon's
AmiBase::epsilon_t epsilon_1= {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_2= {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_3= {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_4= {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_5= {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_6= {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_7= {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_8= {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_9= {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_10={0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_11={0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
AmiBase::epsilon_t epsilon_12={0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0};
AmiBase::epsilon_t epsilon_13={0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
AmiBase::epsilon_t epsilon_14={0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0};
AmiBase::epsilon_t epsilon_15={0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0};
AmiBase::epsilon_t epsilon_16={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0};
AmiBase::epsilon_t epsilon_17={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};


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
AmiBase::g_struct g12(epsilon_12,alpha_12);
AmiBase::g_struct g13(epsilon_13,alpha_13);
AmiBase::g_struct g14(epsilon_14,alpha_14);
AmiBase::g_struct g15(epsilon_15,alpha_15);
AmiBase::g_struct g16(epsilon_16,alpha_16);
AmiBase::g_struct g17(epsilon_17,alpha_17);




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
R0.push_back(g12);
R0.push_back(g13);
R0.push_back(g14);
R0.push_back(g15);
R0.push_back(g16);
R0.push_back(g17);


return R0;

}




AmiBase::ami_vars construct_random_example_J(std::mt19937 &rng){

AmiBase::energy_t energy=random_energy(17, rng);

//energy_t energy={-4,1,1,1,1,1,1,1,1,1,-4,-4,1,1,1,1,1};
		//{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//energy_t energy={-4,1,-1,1,1,-4,1,1,1,1,1,1,1,1,1,1,1};

AmiBase::frequency_t frequency= {std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
			std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,0),
				std::complex<double>(0,M_PI)};

//for(int i=0;i<9;i++){ frequency.push_back(std::complex<double>(0,0));}//frequency.push_back(std::complex<double>(0,0));}

//frequency.push_back(std::complex<double>(0, M_PI));//(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}


AmiBase::energy_t random_energy(int N, std::mt19937 &rng){

std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);




AmiBase::energy_t result;

for (int i=0; i<N;i++){

double kx=dis(rng);
double ky=dis(rng);
double energy=-2.0*(cos(kx)+cos(ky));

result.push_back(energy);


std::cout<< "kx ky and energy are \n "<< kx <<" "<< ky <<" "<< energy <<std::endl;

}


return result;

}






AmiBase::g_prod_t construct_example_debug(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>
// int alpha1[]={1,0,0,0};
// int alpha2[]={0,1,0,0};
// int alpha3[]={0,0,1,0};
// int alpha4[]={1,1,0,-1};
// int alpha5[]={1,0,1,-1};

//defining epsilon's
// int epsilon1[]={1,0};
// int epsilon2[]={1,0};
// int epsilon3[]={1,0};
// int epsilon4[]={0,1};
// int epsilon5[]={0,1};

int alpha1[]={1,0,0,0};
int alpha2[]={0,1,0,0};
int alpha3[]={0,0,1,0};
int alpha4[]={-1,0,1,1};
int alpha5[]={0,-1,1,1};

//defining epsilon's
int epsilon1[]={0,1};
int epsilon2[]={0,1};
int epsilon3[]={0,1};
int epsilon4[]={1,0};
int epsilon5[]={1,0};

AmiBase::alpha_t alpha_1(alpha1, alpha1 + sizeof(alpha1) / sizeof(int) );
AmiBase::alpha_t alpha_2(alpha2, alpha2 + sizeof(alpha2) / sizeof(int) );
AmiBase::alpha_t alpha_3(alpha3, alpha3 + sizeof(alpha3) / sizeof(int) );
AmiBase::alpha_t alpha_4(alpha4, alpha4 + sizeof(alpha4) / sizeof(int) );
AmiBase::alpha_t alpha_5(alpha5, alpha5 + sizeof(alpha5) / sizeof(int) );


AmiBase::epsilon_t epsilon_1(epsilon1, epsilon1+sizeof(epsilon1)/sizeof(int));
AmiBase::epsilon_t epsilon_2(epsilon2, epsilon2+sizeof(epsilon2)/sizeof(int));
AmiBase::epsilon_t epsilon_3(epsilon3, epsilon3+sizeof(epsilon3)/sizeof(int));
AmiBase::epsilon_t epsilon_4(epsilon4, epsilon4+sizeof(epsilon4)/sizeof(int));
AmiBase::epsilon_t epsilon_5(epsilon5, epsilon5+sizeof(epsilon5)/sizeof(int));

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);
AmiBase::g_struct g4(epsilon_4,alpha_4);
AmiBase::g_struct g5(epsilon_5,alpha_5);

AmiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);




return R0;

}

AmiBase::ami_vars construct_ext_example_debug(){



// energy_t energy={-4,1,-1};
AmiBase::energy_t energy={0.6247996574321064,-0.6247996574321064};

AmiBase::frequency_t frequency;

for(int i=0;i<3;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

} 




AmiBase::g_prod_t construct_example2_safe(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


AmiBase::alpha_t alpha_1={1,0,0};
AmiBase::alpha_t alpha_2={0,1,0};
AmiBase::alpha_t alpha_3={-1,1,1};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0};
AmiBase::epsilon_t epsilon_3={0,0,1};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);

// AmiBase::g_prod_t R0={g1,g2,g3};

// OR
AmiBase::g_prod_t R0;
R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);




return R0;

}

AmiBase::ami_vars construct_ext_example2_safe(){


AmiBase::energy_t energy={1.000,-1.00,-1.00};

AmiBase::frequency_t frequency;

for(int i=0;i<2;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}

AmiBase::g_prod_t construct_example3_safe(){

AmiBase::g_prod_t g;

// Setting up G array
// defining alpha's /// std::vector<int>


AmiBase::alpha_t alpha_1={1,-1,1,0};
AmiBase::alpha_t alpha_2={1,-1,1,0};
AmiBase::alpha_t alpha_3={1,0,0,0};
AmiBase::alpha_t alpha_4={0,1,0,0};
AmiBase::alpha_t alpha_5={0,0,1,0};

//defining epsilon's
AmiBase::epsilon_t epsilon_1={1,0,0,0,0};
AmiBase::epsilon_t epsilon_2={0,1,0,0,0};
AmiBase::epsilon_t epsilon_3={0,0,1,0,0};
AmiBase::epsilon_t epsilon_4={0,0,0,1,0};
AmiBase::epsilon_t epsilon_5={0,0,0,0,1};

AmiBase::g_struct g1(epsilon_1,alpha_1);
AmiBase::g_struct g2(epsilon_2,alpha_2);
AmiBase::g_struct g3(epsilon_3,alpha_3);
AmiBase::g_struct g4(epsilon_4,alpha_4);
AmiBase::g_struct g5(epsilon_5,alpha_5);

// AmiBase::g_prod_t R0={g1,g2,g3};

// OR
AmiBase::g_prod_t R0;
R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);




return R0;

}

AmiBase::ami_vars construct_ext_example3_safe(){


AmiBase::energy_t energy={-2,-2,-3.41421,-2.82843,-1.41421};

AmiBase::frequency_t frequency;

for(int i=0;i<3;i++){ frequency.push_back(std::complex<double>(0,0));}

frequency.push_back(std::complex<double>(0,M_PI));

double BETA=1.0;
AmiBase::ami_vars external(energy, frequency,BETA);

return external;

}