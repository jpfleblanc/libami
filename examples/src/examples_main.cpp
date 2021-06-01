
#include "examples.hpp"



int main(int argc, char** argv)
{
	

std::cout<<"---- Running 2nd order example -----"<<std::endl;	
example2();
std::cout<<"----"<<std::endl;
	
std::cout<<std::endl;
std::cout<<"---- Running 4th order multipole example -----"<<std::endl;	
example4();
std::cout<<"----"<<std::endl;

std::cout<<std::endl;
std::cout<<"---- Running 4th order higher multipole example -----"<<std::endl;	
example4hm();
std::cout<<"----"<<std::endl;


std::cout<<std::endl;
std::cout<<"---- Running 9th order example -----"<<std::endl;	
example9();
std::cout<<"----"<<std::endl;
		
	
	
}


void example2(){
	
// class instance
AmiBase ami;

// Problem setup
AmiBase::g_prod_t R0=construct_example();
AmiBase::ami_vars avars=construct_ext_example_Sigma();

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects 
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;

// Integration/Evaluation parameters
double E_REG=0;
int N_INT=2;
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

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	


// Same Problem with Optimization 
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


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	


// Same Problem, term storage type

std::cout<<"-----Constructing AMI term by term-----"<<std::endl;

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

std::cout<<"Construction took "<<d5.count()<<" microseconds"<<std::endl;
std::cout<<"Result has num_terms="<<amiterms.size()<<" compared to standard "<<R_array[N_INT].size()<<std::endl;

	//timing info 
	auto t8=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"evaluation took "<<d6.count()<<" microseconds"<<std::endl;


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list);
	auto t11=std::chrono::high_resolution_clock::now();



	auto t12=std::chrono::high_resolution_clock::now();
std::complex<double> opt_term_val=ami.evaluate(test_amiparms, amiterms, avars,unique, rref, eval_list);
	auto t13=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff7=t11-t10;
	std::chrono::duration<double> diff8=t13-t12;

	std::chrono::microseconds d7=std::chrono::duration_cast<std::chrono::microseconds>(diff7);
	std::chrono::microseconds d8=std::chrono::duration_cast<std::chrono::microseconds>(diff8);	


std::cout<<"OPT term val= "<< opt_term_val<<std::endl;
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d7.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<<d8.count()<<" microseconds"<<std::endl;
	
	
}


void example4(){
	
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
double E_REG=0;//1e-10;
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

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	


// Same Problem with Optimization 
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


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	


// Same Problem, term storage type

std::cout<<"-----Constructing AMI term by term-----"<<std::endl;

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

std::cout<<"Construction took "<<d5.count()<<" microseconds"<<std::endl;
std::cout<<"Result has num_terms="<<amiterms.size()<<" compared to standard "<<R_array[N_INT].size()<<std::endl;

	//timing info 
	auto t8=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"evaluation took "<<d6.count()<<" microseconds"<<std::endl;


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list);
	auto t11=std::chrono::high_resolution_clock::now();


	auto t12=std::chrono::high_resolution_clock::now();
std::complex<double> opt_term_val=ami.evaluate(test_amiparms, amiterms, avars,unique, rref, eval_list);
	auto t13=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff7=t11-t10;
	std::chrono::duration<double> diff8=t13-t12;

	std::chrono::microseconds d7=std::chrono::duration_cast<std::chrono::microseconds>(diff7);
	std::chrono::microseconds d8=std::chrono::duration_cast<std::chrono::microseconds>(diff8);	


std::cout<<"OPT term val= "<< opt_term_val<<std::endl;
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d7.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<<d8.count()<<" microseconds"<<std::endl;
	
	
}


void example4hm(){
	
// class instance 
AmiBase ami;

// Problem setup 
AmiBase::g_prod_t R0=construct_hm_multipole_example();
AmiBase::ami_vars avars=construct_4ord_ext_multipole_example();

	//timing info 
	auto t1=std::chrono::high_resolution_clock::now();
 
// Storage objects  
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;

//Integration/Evaluation parameters 
double E_REG=0;//1e-10;
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

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	


// Same Problem with Optimization 
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


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	


// Same Problem, term storage type

std::cout<<"-----Constructing AMI term by term-----"<<std::endl;

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

std::cout<<"Construction took "<<d5.count()<<" microseconds"<<std::endl;
std::cout<<"Result has num_terms="<<amiterms.size()<<" compared to standard "<<R_array[N_INT].size()<<std::endl;

	//timing info 
	auto t8=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"evaluation took "<<d6.count()<<" microseconds"<<std::endl;


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list);
	auto t11=std::chrono::high_resolution_clock::now();


	auto t12=std::chrono::high_resolution_clock::now();
std::complex<double> opt_term_val=ami.evaluate(test_amiparms, amiterms, avars,unique, rref, eval_list);
	auto t13=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff7=t11-t10;
	std::chrono::duration<double> diff8=t13-t12;

	std::chrono::microseconds d7=std::chrono::duration_cast<std::chrono::microseconds>(diff7);
	std::chrono::microseconds d8=std::chrono::duration_cast<std::chrono::microseconds>(diff8);	


std::cout<<"OPT term val= "<< opt_term_val<<std::endl;
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d7.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<<d8.count()<<" microseconds"<<std::endl;
	
	
}

void example9(){
	
// class instance 
AmiBase ami;

// Problem setup 
AmiBase::g_prod_t R0=construct_example_J();
AmiBase::ami_vars avars=construct_ext_example_J();

	//timing info 
	auto t1=std::chrono::high_resolution_clock::now();
 
// Storage objects  
AmiBase::R_t R_array;
AmiBase::P_t P_array;
AmiBase::S_t S_array;

// Integration/Evaluation parameters 
double E_REG=0;//1e-8;
int N_INT=9;
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

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

// Factorization takes too long for large integrands
/* 
// Same Problem with Optimization 
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


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	
 */

// Same Problem, term storage type

std::cout<<"-----Constructing AMI term by term-----"<<std::endl;

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

std::cout<<"Construction took "<<d5.count()<<" microseconds"<<std::endl;
std::cout<<"Result has num_terms="<<amiterms.size()<<" compared to standard "<<R_array[N_INT].size()<<std::endl;

	//timing info 
	auto t8=std::chrono::high_resolution_clock::now();

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"evaluation took "<<d6.count()<<" microseconds"<<std::endl;

// Factorization takes too long for large integrands
/* 
	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list);
	auto t11=std::chrono::high_resolution_clock::now();


	auto t12=std::chrono::high_resolution_clock::now();
std::complex<double> opt_term_val=ami.evaluate(test_amiparms, amiterms, avars,unique, rref, eval_list);
	auto t13=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff7=t11-t10;
	std::chrono::duration<double> diff8=t13-t12;

	std::chrono::microseconds d7=std::chrono::duration_cast<std::chrono::microseconds>(diff7);
	std::chrono::microseconds d8=std::chrono::duration_cast<std::chrono::microseconds>(diff8);	


std::cout<<"OPT term val= "<< opt_term_val<<std::endl;
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d7.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<<d8.count()<<" microseconds"<<std::endl; */
	
	
}




