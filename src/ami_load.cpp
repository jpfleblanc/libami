#include "ami.hpp"



void AmiCalc::load_solutions(std::string top_directory, solution_set_matrix_t &AMI_MATRIX, int MAX_ORDER, double EREG){
	
for(int ord=0; ord< MAX_ORDER; ord++){	
	
std::string path=top_directory+"/"+std::to_string(ord+1)+"_order/";
std::vector<std::string> files;

std::string S_folder, P_epsfolder, P_alphafolder, R_epsfolder, R_alphafolder, R0_folder, f_folder, eps_folder, mul_folder;
std::string S_file, P_epsfile, P_alphafile, R_epsfile, R_alphafile, R0_file, f_file, eps_file, mul_file;
S_folder=path+"S_m_"+std::to_string(ord+1)+"_txt_files/";
P_epsfolder=path+"P_mnta_m_"+std::to_string(ord+1)+"_txt_files/";
P_alphafolder=path+"P_freq_m_"+std::to_string(ord+1)+"_txt_files/";
R_epsfolder=path+"R_mnta_m_"+std::to_string(ord+1)+"_txt_files/";
R_alphafolder=path+"R_freq_m_"+std::to_string(ord+1)+"_txt_files/";
R0_folder=path+"alpha_m_"+std::to_string(ord+1)+"_txt_files/";
f_folder=path+"f_m_"+std::to_string(ord+1)+"_txt_files/";
eps_folder=path+"epsilon_m_"+std::to_string(ord+1)+"_txt_files/";
mul_folder=path+"mul_m_"+std::to_string(ord+1)+"_txt_files/";

files.clear();
for (const auto & entry : std::experimental::filesystem::directory_iterator(S_folder)){
std::cout << entry.path() << std::endl;	
files.push_back(entry.path());
}
//std::cout<<files.size()<<std::endl;




for(int num=0; num<files.size(); num++){

// file names
S_file=S_folder+"S_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
P_epsfile=P_epsfolder+"P_mnta_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
P_alphafile=P_alphafolder+"P_freq_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
R_epsfile=R_epsfolder+"R_mnta_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
R_alphafile=R_alphafolder+"R_freq_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
R0_file=R0_folder+"alpha_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
f_file=f_folder+"f_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
eps_file=eps_folder+"epsilon_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";
mul_file=mul_folder+"mul_m_"+std::to_string(ord+1)+"_num_"+std::to_string(num+1)+".txt";

// actual variables 
AmiCalc::S_t test_S;
AmiCalc::P_t test_P;
AmiCalc::R_t test_R;
AmiCalc::g_prod_t test_R0;
double test_prefactor;
double test_ereg;
int test_graph_order;
AmiCalc::external_variable_list test_ext;

test_ereg=EREG; // TODO: This needs to already be somewhere rather than hard-coded
test_graph_order=ord+1;

read_text_S_solutions(S_file, test_S);
//std::cout<<"Test_S has size "<< test_S.size()<<std::endl;
//g.ami.write_S_readable(test_S);


read_text_P_solutions(P_epsfile, P_alphafile, test_P);
//write_P_readable(test_P);
//g.ami.print_P(4,test_P);


read_text_R_solutions(R_epsfile, R_alphafile, test_R, test_graph_order);
//write_R_readable(test_R);


read_text_R0(R0_file, eps_file, test_R0);

test_prefactor=load_prefactor(f_file,mul_file, test_graph_order);

// at this point everything is constructed?

AmiCalc::ami_parms test_amiparms(test_graph_order, test_ereg);
//
// me being lazy
AmiCalc::solution_set solution(test_R0, test_S, test_P, test_R, test_amiparms,  test_prefactor);

AMI_MATRIX[ord].push_back(solution);


}	
	
	
}
}



void AmiCalc::read_external(std::string filename, external_variable_list &extern_list){
	
std::ifstream infile_stream;

infile_stream.open(filename);

ext_vars line_variables;

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
std::getline(infile_stream,line);	

while (std::getline(infile_stream, line))
{

//std::cout<<line<<std::endl;

	std::stringstream ss(line);
    double beta, realmu,imagmu;
	double H;
	double kdim;
	double realW, imagW;
	double kx, ky;
	std::string laststring;
	bool read = bool(ss >> laststring);
	if(read){
	beta=std::stod(laststring);
	ss >>  realmu>>imagmu >> H >> kdim;// >> kx >> ky >> realW>> imagW;
	
	// std::cout<<beta<<" "<<mu<<" "<<kdim<<std::endl;//<<" "<<kx<<" "<<ky<<" "<<realW<<" "<<imagW<<" "<<std::endl;
	
    // if (ss >> beta >> mu >> kdim)
    // {
		
		line_variables.BETA_=beta;
		std::complex<double> mu(realmu,imagmu);
		line_variables.MU_=mu;
		line_variables.H_=H;
		line_variables.KDIM_=int(kdim);
		line_variables.external_k_vector_.assign(kdim,0.0);
    // }
	 for (int i=0; i<kdim;i++){
		 
		 ss >> line_variables.external_k_vector_[i];
		 // std::cout<<line_variables.external_k_vector_[i]<<std::endl;;
	 }
	 
	 ss >> realW>>imagW;
	 std::complex<double> freq(realW,imagW);
	 line_variables.external_freq_[0]=freq;
	 
	 
	//std::cout<<"Push back"<<std::endl;
	extern_list.push_back(line_variables);
	}
}
	
}

void AmiCalc::read_text_S_solutions(std::string filename, S_t &s_array){
	
	std::ifstream infile_stream;

infile_stream.open(filename);

Si_t collect_si;
sign_t collect_signs;

ext_vars line_variables;

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
//std::getline(infile_stream,line);

int si_index=0;	
int sign_index=0;

while (std::getline(infile_stream, line))
{
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);
std::string left="[";
std::string right="]";

std::string item;	

while( ! ss.eof()){


ss >> item;
//std::cout<<"Item is "<< item <<std::endl;

if(item!=left && item!=right){
	
collect_signs.push_back(std::stod(item));
	
	
}



if(item==right){
	if(collect_signs.size()!=0){
	collect_si.push_back(collect_signs);
	collect_signs.clear();
	}
}
}


if(collect_si.size() != 0){	
si_index++;	
s_array.push_back(collect_si);
collect_si.clear();
}
	
	
	
}
	
	
infile_stream.close();	
	
}


void AmiCalc::read_text_P_solutions(std::string eps_filename,std::string alpha_filename, P_t &p_array){
	
std::ifstream infile_stream;

// read epsilons first
infile_stream.open(eps_filename);

Pi_t collect_pi;
pole_array_t collect_pole_array;
pole_struct collect_eps;

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<eps_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
//std::getline(infile_stream,line);

int pi_index=0;	
int pa_index=0;

while (std::getline(infile_stream, line))
{
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);
std::string left="[";
std::string right="]";

std::string current, next;	
ss >> current;

while( ! ss.eof()){

ss >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;

if(current!="[" && current !="]"){ 

collect_eps.eps_.push_back( std::stoi(current));

}
if(current=="]" && next=="["){

if(collect_eps.eps_.size()!=0){	
collect_pole_array.push_back(collect_eps);
collect_eps.eps_.clear();	
}
}

if(current=="]" && next=="]"){
	if(collect_eps.eps_.size()!=0){
collect_pole_array.push_back(collect_eps);
collect_eps.eps_.clear();
	}
	if(collect_pole_array.size()!=0){
collect_pi.push_back(collect_pole_array);
collect_pole_array.clear();
	}
}


current=next;
}


pi_index++;

if(collect_pi.size()!=0){
p_array.push_back(collect_pi);
collect_pi.clear();	
}
}

	
	
infile_stream.close();	

/////// epsilons now loaded. now iterate through a collect alphas

int ln, bn, an;
ln=0;
bn=0;
an=0;

infile_stream.open(alpha_filename);

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<alpha_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
while (std::getline(infile_stream, line))
{
	
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);
std::string left="[";
std::string right="]";

std::string current, next;	
ss >> current;


while( ! ss.eof()){

ss >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;

if(current!="[" && current !="]"){ 

collect_eps.alpha_.push_back( std::stoi(current));

}
if(current=="]" && next=="["){

if(collect_eps.alpha_.size()!=0){
//std::cout<<ln<<" "<<bn<<" "<< an<<std::endl;	
p_array[ln][bn][an].alpha_=collect_eps.alpha_;
collect_eps.alpha_.clear();	

an++;
}
}

if(current=="]" && next=="]"){
	if(collect_eps.alpha_.size()!=0){
//		std::cout<<ln<<" "<<bn<<" "<< an<<std::endl;
p_array[ln][bn][an].alpha_=collect_eps.alpha_;
collect_eps.alpha_.clear();	

an=0;
bn++;
}
}


current=next;

}
ln++;
an=0;
bn=0;
}





	
}

double AmiCalc::load_prefactor(std::string filename, std::string mul_file, int order){
	
std::ifstream infile_stream;

infile_stream.open(filename);	
	
if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
//std::getline(infile_stream,line);

double fermi_loops;

while (std::getline(infile_stream, line))
{
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);	
	
ss >> 	fermi_loops;
	
}

double mult=load_mul(mul_file);


double output=pow(-1, fermi_loops + double(order))*mult;

return output;	
}

double AmiCalc::load_mul(std::string filename){
	
double output;

std::ifstream infile_stream;

//std::cout<<"Opening mul file "<<filename<<std::endl;

infile_stream.open(filename);	
	
if(infile_stream.fail()) // checks to see if file opended 
    { 
	output=1.0;
//	std::cout<<"Didn't open mul_file"<<std::endl;
    } else{	

//std::cout<<"Opened mul_file"<<std::endl;	
std::string line;
//std::getline(infile_stream,line);


while (std::getline(infile_stream, line))
{
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);	
	
ss >> 	output;
	
}
	}	
	
	
std::cout<<"Returning mul="<<output<<std::endl;	
	return output;
}

void AmiCalc::read_text_R0(std::string alpha_filename,std::string eps_filename, g_prod_t &R0){

g_struct collect_g;

std::ifstream infile_stream;
infile_stream.open(alpha_filename);

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<alpha_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;

std::getline(infile_stream, line);
std::stringstream ss(line);

std::string current, next;	
ss >> current;


while( ! ss.eof()){

ss >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;
if(current!="[" && current !="]"){ 

collect_g.alpha_.push_back( std::stoi(current));

}

if(current=="]" && next=="["){
	if(collect_g.alpha_.size()!=0){
R0.push_back(collect_g);
	collect_g.alpha_.clear();}	
}

if(current=="]" && next=="]"){
	if(collect_g.alpha_.size()!=0){
R0.push_back(collect_g);
	collect_g.alpha_.clear();	}
}

current=next;

}

infile_stream.close();

std::cout<<"Actual"<<std::endl;
std::cout<<eps_filename<<std::endl;

//std::cout<<"R0 has "<< R0.size()<<" Green's functions"<<std::endl;

// for(int i=0; i< R0.size(); i++){

// R0[i].eps_.resize(R0.size(),0);
// R0[i].eps_[i]=1;	
	
// }

// for(int i=0; i<R0.size(); i++){
// for(int j=i; j< R0.size(); j++){

// if(i!=j){

// if(R0[j].alpha_==R0[i].alpha_){
// std::cout<<"This triggered"<<std::endl;
// R0[j].eps_=R0[i].eps_;
// }	
	
// }	

// }	
// }

// epsilon loading


infile_stream.open(eps_filename);

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<eps_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	


std::getline(infile_stream, line);
std::stringstream ss2(line);


ss2 >> current;

int ind=0;
while( ! ss2.eof()){

ss2 >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;
if(current!="[" && current !="]"){ 

R0[ind].eps_.push_back( std::stoi(current));

}

if(current=="]" && next=="["){
	ind++;
}

current=next;

}



}

void AmiCalc::read_text_R_solutions(std::string eps_filename,std::string alpha_filename, R_t &r_array, int size){

g_prod_array_t collect_ri;
g_prod_t collect_gprod;
g_struct collect_g;

// fill dummy R elements that don't  get evaluated 


for(int i=0; i<size; i++){
	
r_array.push_back(collect_ri); // this is just pushing back empty structures that don't get looked at 
	
}
	
std::ifstream infile_stream;

// read epsilons first
infile_stream.open(eps_filename);



if(infile_stream.fail()) // checks to see if file opended 
    { 
	  std::cout<<eps_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
//std::getline(infile_stream,line);


while (std::getline(infile_stream, line))
{
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);

std::string current, next;	
ss >> current;

while( ! ss.eof()){

ss >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;

if(current!="[" && current !="]"){ 

collect_g.eps_.push_back( std::stoi(current));

}
if(current=="]" && next=="["){

if(collect_g.eps_.size()!=0){	
collect_gprod.push_back(collect_g);
collect_g.eps_.clear();	
}
}

if(current=="]" && next=="]"){

if(collect_g.eps_.size()!=0){	
collect_gprod.push_back(collect_g);
collect_g.eps_.clear();	
}
}


current=next;
}


if(collect_gprod.size()!=0){
collect_ri.push_back(collect_gprod);
collect_gprod.clear();	
}





}

r_array.push_back(collect_ri);
collect_ri.clear();



//std::cout<<"Loading R-alphas"<<std::endl;



	
infile_stream.close();	

/////// epsilons now loaded. now iterate through a collect alphas

int ln, bn, an;
ln=0;
bn=0;
an=0;

int last=size;

infile_stream.open(alpha_filename);

if(infile_stream.fail()) // checks to see if file opended 
    { 
	  std::cout<<alpha_filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
while (std::getline(infile_stream, line))
{
	
	
//std::cout<<line<<std::endl;

std::stringstream ss(line);

std::string current, next;	
ss >> current;

while( ! ss.eof()){

ss >> next;

//std::cout<<"Current and next: "<< current<<" "<<next<<std::endl;

if(current!="[" && current !="]"){ 
//r_array[r_array.size()-1][i][j].eps_[k]
collect_g.alpha_.push_back( std::stoi(current));

}
if(current=="]" && next=="["){

if(collect_g.alpha_.size()!=0){
//std::cout<<ln<<" "<<bn<<" "<< an<<std::endl;	
r_array[last][ln][bn].alpha_=collect_g.alpha_;
collect_g.alpha_.clear();	

bn++;
}
}

if(current=="]" && next=="]"){
	if(collect_g.alpha_.size()!=0){
//		std::cout<<ln<<" "<<bn<<" "<< an<<std::endl;	
r_array[last][ln][bn].alpha_=collect_g.alpha_;
collect_g.alpha_.clear();	

bn=0;
}
}


current=next;

}
ln++;
bn=0;
}


	infile_stream.close();
}








