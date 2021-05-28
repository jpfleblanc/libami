
#include "examples.hpp"



int main(int argc, char** argv)
{
	

std::cout<<"---- Running 2nd order example -----"<<std::endl;	
example2();
std::cout<<"----"<<std::endl;
	

std::cout<<"---- Running 4th order multipole example -----"<<std::endl;	
example4();
std::cout<<"----"<<std::endl;


std::cout<<"---- Running 4th order higher multipole example -----"<<std::endl;	
example4hm();
std::cout<<"----"<<std::endl;



std::cout<<"---- Running 9th order example -----"<<std::endl;	
example9();
std::cout<<"----"<<std::endl;
		
	
	
}


void example2(){
	

AmiBase ami;

AmiBase::g_prod_t R0=construct_example();
AmiBase::ami_vars avars=construct_ext_example_Sigma();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;


double E_REG=0;//1e-5;
int N_INT=2;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);


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



auto t3=std::chrono::high_resolution_clock::now();

AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

ami.factorize_Rn(R_array.back(), unique, rref, eval_list);

auto t4=std::chrono::high_resolution_clock::now();	
	
std::complex<double> opt_calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars, unique, rref, eval_list);	
auto t5=std::chrono::high_resolution_clock::now();	


std::chrono::duration<double> diff3=t4-t3;
std::chrono::duration<double> diff4=t5-t4;

std::chrono::microseconds d3=std::chrono::duration_cast<std::chrono::microseconds>(diff3);
std::chrono::microseconds d4=std::chrono::duration_cast<std::chrono::microseconds>(diff4);	


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
	
std::cout<<"Optimization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	

std::cout<<"Constructing AMI term by term "<<std::endl;

AmiBase::terms amiterms;

ami.construct(N_INT, R0, amiterms);

std::cout<<"Result has num_terms="<<amiterms.size()<<std::endl; //" compared to standard "<<R_array[N_INT].size()<<std::endl;

for(int i=0; i<R_array.size(); i++){
std::cout<<"R["<<i<<"]->"<<R_array[i].size()<<std::endl;	
	
	
}
	
	
	
}


void example4(){
	

AmiBase ami;

AmiBase::g_prod_t R0=construct_multipole_example();
AmiBase::ami_vars avars=construct_4ord_ext_multipole_example();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;


double E_REG=0;//1e-10;
int N_INT=4;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);


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


	
auto t3=std::chrono::high_resolution_clock::now();

AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

ami.factorize_Rn(R_array.back(), unique, rref, eval_list);

auto t4=std::chrono::high_resolution_clock::now();	
	
std::complex<double> opt_calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars, unique, rref, eval_list);	
auto t5=std::chrono::high_resolution_clock::now();	


std::chrono::duration<double> diff3=t4-t3;
std::chrono::duration<double> diff4=t5-t4;

std::chrono::microseconds d3=std::chrono::duration_cast<std::chrono::microseconds>(diff3);
std::chrono::microseconds d4=std::chrono::duration_cast<std::chrono::microseconds>(diff4);	


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
	
std::cout<<"Optimization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	
	


std::cout<<"Constructing AMI term by term "<<std::endl;

AmiBase::terms amiterms;

ami.construct(N_INT, R0, amiterms);

std::cout<<"Result has num_terms="<<amiterms.size()<<std::endl; //" compared to standard "<<R_array[N_INT].size()<<std::endl;

for(int i=0; i<R_array.size(); i++){
std::cout<<"R["<<i<<"]->"<<R_array[i].size()<<std::endl;	
	
	
}
	
	
	
	
}


void example4hm(){
	

AmiBase ami;

AmiBase::g_prod_t R0=construct_hm_multipole_example();
AmiBase::ami_vars avars=construct_4ord_ext_multipole_example();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;


double E_REG=0;//1e-10;
int N_INT=4;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);


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


	
auto t3=std::chrono::high_resolution_clock::now();

AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

ami.factorize_Rn(R_array.back(), unique, rref, eval_list);

auto t4=std::chrono::high_resolution_clock::now();	
	
std::complex<double> opt_calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars, unique, rref, eval_list);	
auto t5=std::chrono::high_resolution_clock::now();	


std::chrono::duration<double> diff3=t4-t3;
std::chrono::duration<double> diff4=t5-t4;

std::chrono::microseconds d3=std::chrono::duration_cast<std::chrono::microseconds>(diff3);
std::chrono::microseconds d4=std::chrono::duration_cast<std::chrono::microseconds>(diff4);	


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
	
std::cout<<"Optimization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	
	
}

void example9(){
	

AmiBase ami;

AmiBase::g_prod_t R0=construct_example_J();
AmiBase::ami_vars avars=construct_ext_example_J();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t t1, t2, t_end;
// clock_t tc,te;
// int top=10000;

auto t1=std::chrono::high_resolution_clock::now();
 // for(int i=0; i< top; i++)
// {
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;


double E_REG=0;//1e-8;
int N_INT=9;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);


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




