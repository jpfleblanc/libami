
#include "examples.hpp"



int main(int argc, char** argv)
{
	

std::cout<<"---- Running 2nd order example -----"<<std::endl;	
example2();
std::cout<<"----"<<std::endl;
	

std::cout<<"---- Running 4th order multipole example -----"<<std::endl;	
example4();
std::cout<<"----"<<std::endl;


std::cout<<"---- Running 9th order example -----"<<std::endl;	
example9();
std::cout<<"----"<<std::endl;
		
	
	
}


void example2(){
	

AmiCalc ami;

AmiCalc::g_prod_t R0=ami.construct_example();
AmiCalc::ami_vars avars=ami.construct_ext_example_Sigma();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiCalc::R_t R_array;
AmiCalc::P_t P_array;
AmiCalc::S_t S_array;


double E_REG=0;//1e-5;
int N_INT=2;
AmiCalc::ami_parms test_amiparms(N_INT, E_REG);


ami.construct(test_amiparms, R0, R_array, P_array, S_array);
auto t2=std::chrono::high_resolution_clock::now();
// ami.print_final(N_INT, R_array, P_array, S_array);

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	
	
auto t_end=std::chrono::high_resolution_clock::now();

std::chrono::duration<double> diff1=t2-t1;
std::chrono::duration<double> diff2=t_end-t2;

std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;

std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

	
	
	
	
}


void example4(){
	

AmiCalc ami;

AmiCalc::g_prod_t R0=ami.construct_multipole_example();
AmiCalc::ami_vars avars=ami.construct_4ord_ext_multipole_example();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiCalc::R_t R_array;
AmiCalc::P_t P_array;
AmiCalc::S_t S_array;


double E_REG=0;//1e-10;
int N_INT=4;
AmiCalc::ami_parms test_amiparms(N_INT, E_REG);


ami.construct(test_amiparms, R0, R_array, P_array, S_array);
auto t2=std::chrono::high_resolution_clock::now();
// ami.print_final(N_INT, R_array, P_array, S_array);

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	
	
auto t_end=std::chrono::high_resolution_clock::now();

std::chrono::duration<double> diff1=t2-t1;
std::chrono::duration<double> diff2=t_end-t2;

std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;

std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	


	
}

void example9(){
	

AmiCalc ami;

AmiCalc::g_prod_t R0=ami.construct_example_J();
AmiCalc::ami_vars avars=ami.construct_ext_example_J();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiCalc::R_t R_array;
AmiCalc::P_t P_array;
AmiCalc::S_t S_array;


double E_REG=0;//1e-8;
int N_INT=9;
AmiCalc::ami_parms test_amiparms(N_INT, E_REG);


ami.construct(test_amiparms, R0, R_array, P_array, S_array);
auto t2=std::chrono::high_resolution_clock::now();
// ami.print_final(N_INT, R_array, P_array, S_array);

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	
	
auto t_end=std::chrono::high_resolution_clock::now();

std::chrono::duration<double> diff1=t2-t1;
std::chrono::duration<double> diff2=t_end-t2;

std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;

std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	


	
	
	
}




