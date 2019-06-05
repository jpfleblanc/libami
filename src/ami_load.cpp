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

std::cout<<line<<std::endl;

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

