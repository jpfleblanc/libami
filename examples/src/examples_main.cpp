
#include "examples.hpp"



int main(int argc, char** argv)
{
	

example2();
example1_bose();
example4();
example9();

		
//safe_example();	
	
}

void safe_example(){
  
  std::cout<<std::endl<<"-_-_-_ Example - Second Order _-_-_-"<<std::endl<<std::endl;	
	
//START example	
std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
// class instance
AmiBase ami;

// Problem setup (see ami_example.cpp)
AmiBase::g_prod_t R0=construct_example2_safe(); // Sets initial integrand 
AmiBase::ami_vars avars=construct_ext_example2_safe(); // Sets 'external' parameter values 

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects for S,P,R 
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=2;  // Number of Matsubara sums to perform
AmiBase::ami_parms test_amiparms(N_INT, E_REG);

//Construction Stage
ami.construct(test_amiparms, R0, R_array, P_array, S_array);  // Populates S,P,R with solution 

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	// Evaluate integrand for parameters in 'avars'
	
	//timing info
	auto t_end=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::duration<double> diff2=t_end-t2;

	std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	
 

auto t_otf1=std::chrono::high_resolution_clock::now(); 
 
for( double n=0.9; n<1.1; n+=0.001){ 
avars.energy_[0]=n; 
calc_result=ami.evaluate_otf(test_amiparms,R_array, P_array, S_array,  avars); 

std::cout<<n<<" "<<calc_result.real()<<" "<< calc_result.imag()<<std::endl;

}
auto t_otf2=std::chrono::high_resolution_clock::now();
std::chrono::duration<double> diffotf=t_otf2-t_otf1;
std::chrono::microseconds dotf=std::chrono::duration_cast<std::chrono::microseconds>(diffotf);
std::cout<<"Result for OTF evaluation was "<< calc_result<<std::endl;
std::cout<<"Evaluation took "<< dotf.count() <<" microseconds"<<std::endl;	
  
}


void example2(){
std::cout<<std::endl<<"-_-_-_ Example - Second Order _-_-_-"<<std::endl<<std::endl;	
	
//START example	
std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
// class instance
AmiBase ami;

// Problem setup (see ami_example.cpp)
AmiBase::g_prod_t R0=construct_example2(); // Sets initial integrand 
AmiBase::ami_vars avars=construct_ext_example2(); // Sets 'external' parameter values 

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects for S,P,R 
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=2;  // Number of Matsubara sums to perform
AmiBase::ami_parms test_amiparms(N_INT, E_REG);

//Construction Stage
ami.construct(test_amiparms, R0, R_array, P_array, S_array);  // Populates S,P,R with solution 

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	// Evaluate integrand for parameters in 'avars'
	
	//timing info
	auto t_end=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::duration<double> diff2=t_end-t2;

	std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

//END Example


// START Same Problem with Optimization 
	auto t3=std::chrono::high_resolution_clock::now();

// Storage Structures
AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

// Take existing solution from first part and factorize it 
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
//END Optimization example 

//START Term Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing AMI term by term-----"<<std::endl;

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
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case 
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d6.count()<<" microseconds"<<std::endl;


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
	
//END Term Example 	

//START Optimization for Term construction 	
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list); 
	auto t11=std::chrono::high_resolution_clock::now();



	auto t12=std::chrono::high_resolution_clock::now();
// Evaluate optimized term-by-term construction 	
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


// This example is virtually identical to the fermionic case except for a single line
// AmiBase::ami_parms test_amiparms(N_INT, E_REG); -> AmiBase::ami_parms test_amiparms(N_INT, E_REG,1);  Final argument is 0 for fermionic(default) and 1 for bosonic. 
void example1_bose(){
std::cout<<std::endl<<"-_-_-_ Example - First Order For Bosonic _-_-_-"<<std::endl<<std::endl;	
	
//START example	
std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
// class instance
AmiBase ami;

// Problem setup (see ami_example.cpp)
AmiBase::g_prod_t R0=construct_example1_bose(); // Sets initial integrand 
AmiBase::ami_vars avars=construct_ext_example1_bose(); // Sets 'external' parameter values 

	//timing info
	auto t1=std::chrono::high_resolution_clock::now();

// Storage objects for S,P,R 
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;

// Integration/Evaluation parameters
double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
int N_INT=1;  // Number of Matsubara sums to perform
AmiBase::graph_type bose=AmiBase::Pi_phuu;
AmiBase::ami_parms test_amiparms(N_INT, E_REG,bose);

//Construction Stage
ami.construct(test_amiparms, R0, R_array, P_array, S_array);  // Populates S,P,R with solution 

	//timing info 
	auto t2=std::chrono::high_resolution_clock::now();

// Evaluate the result
std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);	// Evaluate integrand for parameters in 'avars'
	
	//timing info
	auto t_end=std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff1=t2-t1;
	std::chrono::duration<double> diff2=t_end-t2;

	std::chrono::microseconds d=std::chrono::duration_cast<std::chrono::microseconds>(diff1);
	std::chrono::microseconds d2=std::chrono::duration_cast<std::chrono::microseconds>(diff2);

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

//END Example


// START Same Problem with Optimization 
	auto t3=std::chrono::high_resolution_clock::now();

// Storage Structures
AmiBase::g_prod_t unique;
AmiBase::R_ref_t rref;
AmiBase::ref_eval_t eval_list;

// Take existing solution from first part and factorize it 
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
//END Optimization example 

//START Term Example 
// Same Problem, using ami_term storage type
std::cout<<std::endl<<"-----Constructing AMI term by term-----"<<std::endl;

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
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars); // Evaluate the term-by-term result for external values in 'avars'
	
	//timing info
	auto t9=std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff6=t9-t8;
	std::chrono::microseconds d6=std::chrono::duration_cast<std::chrono::microseconds>(diff6);

std::cout<<"Term result was "<< term_val<<std::endl;
std::cout<<"Evaluation took "<<d6.count()<<" microseconds"<<std::endl;


	//timing info
	auto t10=std::chrono::high_resolution_clock::now();
	
//END Term Example 	

//START Optimization for Term construction 	
//Factorized form of term-by-term construction 
ami.factorize_terms(amiterms, unique, rref, eval_list); 
	auto t11=std::chrono::high_resolution_clock::now();



	auto t12=std::chrono::high_resolution_clock::now();
// Evaluate optimized term-by-term construction 	
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
std::cout<<std::endl<<"-_-_-_ Example - Fourth Order _-_-_-"<<std::endl<<std::endl;

//START example	
std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
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

std::cout<<"Result was "<< calc_result<<std::endl;
std::cout<<"Construction took "<<d.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d2.count()<<" microseconds"<<std::endl;	

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


std::cout<<"Optimized result was "<< opt_calc_result<<std::endl;	
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d3.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<< d4.count()<<" microseconds"<<std::endl;	
//END Optimization example

//START Term Example
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


std::cout<<"OPT term val= "<< opt_term_val<<std::endl;
std::cout<<"Factorization returned "<<unique.size()<<" unique Green's functions and took "<<d7.count()<<" microseconds"<<std::endl;
std::cout<<"Evaluation took "<<d8.count()<<" microseconds"<<std::endl;

	
}



void example9(){
std::cout<<std::endl<<"-_-_-_ Example - Ninth Order _-_-_-"<<std::endl<<std::endl;

//START example	
std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
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




/* 
void example2_spec(){
	
// class instance
NewAmiCalc amicalc;
AmiBase ami;

// Problem setup
AmiBase::g_prod_t R0=construct_example2();
AmiBase::ami_vars avars=construct_spec_example();

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
amicalc.print_final(3, R_array, P_array, S_array);

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




std::vector<double> xi_list;
int max=1;
double xcut=4.0;
int gsize=3;

avars.gamma_=0.4;
std::complex<double> sum(0,0);

// for(int n=0; n<max; n++){
// amicalc.randomize_xi(xi_list,gsize,xcut);
std::complex<double> spcalc_result(0,0);

xi_list.resize(gsize);
xi_list[0]=1.0;
xi_list[1]=0.5;
xi_list[2]=0.2;

spcalc_result=amicalc.evaluate_simple_spectral(test_amiparms,R_array, P_array, S_array,  avars, unique, rref, eval_list,  xi_list);

sum+=spcalc_result;

// }

sum=sum;///double(max)*std::pow(2.0*xcut,gsize);
std::cout<<"Got "<<sum<<std::endl;



return;





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
std::cout<<"Result has num_terms="<<amiterms.size()<<std::endl;

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
	
std::cout<<"Simple Spectral evaluation"<<std::endl;




std::complex<double> spec_sum(0,0);

for(int n=0; n<max; n++){


// AmiBase::Pi_t unique_poles;
				
		// graph.ami.collect_spectral_poles(GG_AMI_MATRIX[ord][num][gnum].Unique, GG_AMI_MATRIX[ord][num][gnum].Unique_poles);

std::vector<double> XI_LIST;

amicalc.randomize_xi(XI_LIST,gsize,xcut);

std::complex<double> this_val=amicalc.evaluate_terms_simple_spectral(test_amiparms, amiterms, avars,unique, rref, eval_list, XI_LIST);

spec_sum+= this_val;
	
	
	
}

spec_sum =spec_sum/double(max)*std::pow(2.0*xcut,gsize);

std::cout<<"Average of the simple spectral eval is "<< spec_sum<<std::endl;	
	
	
	
} */