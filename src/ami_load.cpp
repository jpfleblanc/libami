#include "ami.hpp"



void AmiCalc::read_external(std::string filename, external_variable_list &extern_list){
	
std::ifstream infile_stream;

infile_stream.open(filename);

ext_vars line_variables;

if(infile_stream.fail()) // checks to see if file opended 
    { 
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
std::getline(infile_stream,line);	

while (std::getline(infile_stream, line))
{

//std::cout<<line<<std::endl;

	std::stringstream ss(line);
    double beta, mu;
	int kdim;
	double realW, imagW;
	double kx, ky;
	
	
	ss >> beta >> mu >> kdim;// >> kx >> ky >> realW>> imagW;
	
//	std::cout<<beta<<" "<<mu<<" "<<kdim<<std::endl;//<<" "<<kx<<" "<<ky<<" "<<realW<<" "<<imagW<<" "<<std::endl;
	
    // if (ss >> beta >> mu >> kdim)
    // {
		
		line_variables.BETA_=beta;
		line_variables.MU_=mu;
		line_variables.KDIM_=kdim;
		line_variables.external_k_vector_.assign(kdim,0.0);
    // }
	 for (int i=0; i<kdim;i++){
		 
		 ss >> line_variables.external_k_vector_[i];
	 }
	 
	 ss >> realW>>imagW;
	 line_variables.external_freq_.push_back(std::complex<double>(realW,imagW));
	 
	 
//	 std::cout<<line_variables.BETA_<<" "<<line_variables.MU_<<" "<<line_variables.KDIM_<<" "<<line_variables.external_k_vector_[0]<<" "<<line_variables.external_k_vector_[1]<<" "<<line_variables.external_freq_[0].real()<<" "<<line_variables.external_freq_[0].imag()<<" "<<std::endl;
	
	 // if (ss >>	line_variables.external_k_vector_[i]){
		
	 // }
	 // }
	 
	 // if ( ss>> realW >> imagW){
		 
		 // line_variables.external_freq_.push_back(std::complex<double>(realW,imagW));
		 
		 
	 // }
	
	extern_list.push_back(line_variables);
	
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
      throw std::runtime_error("Could not open input file");
    } 	
	
while (std::getline(infile_stream, line))
{
	
	
std::cout<<line<<std::endl;

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

double AmiCalc::load_prefactor(std::string filename, int order){
	
std::ifstream infile_stream;

infile_stream.open(filename);	
	
if(infile_stream.fail()) // checks to see if file opended 
    { 
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


double output=pow(-1, fermi_loops + double(order));

return output;	
}


void AmiCalc::read_text_R_solutions(std::string eps_filename,std::string alpha_filename, R_t &r_array, int size){

g_prod_array_t collect_ri;
g_prod_t collect_gprod;
g_struct collect_g;

// fill dummy R elements that don't  get evaluated 


for(int i=0; i<size-1; i++){
	
r_array.push_back(collect_ri); // this is just pushing back empty structures that don't get looked at 
	
}
	
std::ifstream infile_stream;

// read epsilons first
infile_stream.open(eps_filename);



if(infile_stream.fail()) // checks to see if file opended 
    { 
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

int last=size-1;

infile_stream.open(alpha_filename);

if(infile_stream.fail()) // checks to see if file opended 
    { 
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








