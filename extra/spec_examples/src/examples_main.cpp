
#include "examples.hpp"



int main(int argc, char** argv)
{
	


// std::cout<<"---- Running 2nd order example -----"<<std::endl;	
// example2();
// std::cout<<"----"<<std::endl;


std::cout<<"---- Running simplest example with no spectral representation -----"<<std::endl;	
example_simple_nospec();
std::cout<<"----"<<std::endl;

std::cout<<"---- Running simplest example with spectral representation -----"<<std::endl;	
example_simple_spec();
std::cout<<"----"<<std::endl;

// std::cout<<"---- Running 2nd order spectral example -----"<<std::endl;	
// example2_spec();
// std::cout<<"----"<<std::endl;

	
// std::cout<<std::endl;
// std::cout<<"---- Running 4th order multipole example -----"<<std::endl;	
// example4();
// std::cout<<"----"<<std::endl;

// std::cout<<std::endl;
// std::cout<<"---- Running 4th order higher multipole example -----"<<std::endl;	
// example4hm();
// std::cout<<"----"<<std::endl;


// std::cout<<std::endl;
// std::cout<<"---- Running 9th order example -----"<<std::endl;	
// example9();
// std::cout<<"----"<<std::endl;
		
	
	
}


void example_simple_spec(){
	
std::ofstream file;
file.open("simple_spec_outfile.dat",  std::ofstream::out );	
	

AmiBase ami;
NewAmiCalc amicalc;

AmiBase::g_prod_t R0=construct_simple_example();	

std::vector< std::complex<double>> results;

//simplified storage type 
AmiBase::terms amiterms;
double E_REG=0;
int N_INT=0;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);
// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

// std::cout<<"Terms is of length "<<amiterms.size()<<std::endl;

double xi_cut=10;
double steps=100;

for(double w=-8; w<8; w=w+0.1){

AmiBase::ami_vars avars=ext_example_simple(w, 0.0);

avars.gamma_=0.1;

std::vector<double> xi_list;

xi_list.resize(amiterms[0].g_list.size());

double xstep=xi_cut/33331;

std::complex<double> sums(0,0);
int count=0;

for(double this_x=-xi_cut; this_x<xi_cut; this_x+=xstep){ 

xi_list[0]=this_x;
//Evaluate term-by-term solution at one xi value 
// std::complex<double> evaluate_sp_term(AmiBase::ami_parms &parms, AmiBase::term &ami_term, AmiBase::ami_vars &external, std::vector<double> &xi_list, double &xi_cut);

std::complex<double> term_val=amicalc.evaluate_sp_term(test_amiparms, amiterms[0], avars, xi_list, xi_cut);	

sums+= term_val*xstep;
std::cout<<"w x and term "<< w<<" "<<this_x<<" "<< term_val.real()<<" "<< term_val.imag()<<std::endl;

}
results.push_back(sums);	
file<< w <<" "<< sums.real()<<" "<< sums.imag()<<std::endl;	
	
}	

file.close();

	
	
}


void example_simple_nospec(){
	
std::ofstream file;
file.open("simple_nospec_outfile.dat",  std::ofstream::out );	
	

AmiBase ami;

AmiBase::g_prod_t R0=construct_simple_example();	

std::vector< std::complex<double>> results;

//simplified storage type 
AmiBase::terms amiterms;
double E_REG=0;
int N_INT=0;
AmiBase::ami_parms test_amiparms(N_INT, E_REG);
// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

// std::cout<<"Terms is of length "<<amiterms.size()<<std::endl;

for(double w=-8; w<8; w=w+0.1){

AmiBase::ami_vars avars=ext_example_simple(w, 0.2);

//Evaluate term-by-term solution 
std::complex<double> term_val=ami.evaluate(test_amiparms, amiterms, avars);	

// std::cout<<"Returned value of "<<term_val<<std::endl;
	results.push_back(term_val);
	
file<< w <<" "<< term_val.real()<<" "<< term_val.imag()<<std::endl;	
	
}	

file.close();

	
	
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


void example2_spec(){
	
// class instance
NewAmiCalc amicalc;
AmiBase ami;

// Problem setup
AmiBase::g_prod_t R0=construct_example();
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




