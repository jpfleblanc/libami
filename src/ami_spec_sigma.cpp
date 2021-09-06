#include "ami_spec.hpp"



void AmiSpec::find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec){
  std::vector<double>  comparison_vector;
  comparison_vector=vec;
  for (int n=0;n<vec.size();n++){
    comparison_vector[n]=std::abs(comparison_vector[n]-point);
  }
  //for(int n=0;n<comparison_vector.size();n++){std::cout << "comparison: " << comparison_vector[n]<<std::endl;}
  auto minimum_dif_iter = std::min_element(comparison_vector.begin(),comparison_vector.end());
  double minimum_dif=*minimum_dif_iter;//gives minimum difference between point and contents of vec
  //std::cout<< "\nminimum_dif:  " << minimum_dif << std::endl;
  int minimum_dif_index=0;


  for (int n=0;n<vec.size();n++){
    if (comparison_vector[n]==minimum_dif) minimum_dif_index=n;}
	double minimum2=comparison_vector[0];//what if index=0 is closest point?
	int  minimum2_dif_index=0;
	if(minimum_dif_index==0) {
		minimum2=comparison_vector[1];
		minimum2_dif_index=1;
	}
	for (int n=0;n<vec.size();n++){
    if(comparison_vector[n]<minimum2 && n!=minimum_dif_index) {minimum2=comparison_vector[n]; minimum2_dif_index=n;}

  }

  if (vec[minimum_dif_index]>vec[minimum2_dif_index]) {closest_lt=vec[minimum2_dif_index];closest_gt=vec[minimum_dif_index];}
  else {closest_gt=vec[minimum2_dif_index];closest_lt=vec[minimum_dif_index];}

}
std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X)
{
//if(abs(X)>xi_cutoff) return std::complex<double> (0.0,0.0);
NewAmiCalc::k_vector_t k_copy=k;
for(int for_counter=0;for_counter<k_copy.size();for_counter++){
	if(abs(k_copy[for_counter])>M_PI){
		double modulo_shifted=k_copy[for_counter]+M_PI;
		k_copy[for_counter]=-M_PI+modulo_shifted-floor(modulo_shifted/(2*M_PI))*2*M_PI;//fmod(k_copy[for_counter])
	}

}

  //std::vector<double> point_wanted_modulo;
  //if(X>M_PI) point_wanted_modulo.push_back(fmod())
  //identify 8 closest points
  double freq_lt=AMI_spec_freq_vector_simplified[0];
  double freq_gt=AMI_spec_freq_vector_simplified[0];
  double kx_lt=AMI_spec_kx_vector_simplified[0];
  double kx_gt=AMI_spec_kx_vector_simplified[0];
  double ky_lt=AMI_spec_ky_vector_simplified[0];
  double ky_gt=AMI_spec_ky_vector_simplified[0];

    //if(freq_vector_simplified[n] > freq_lt && abs(X - freq_vector_simplified[n]) < abs(X - freq_lt)) freq_lt = freq_vector_simplified[n];
    //if(freq_vector_simplified[n] < freq_gt && abs(X - freq_vector_simplified[n]) < abs(X - freq_gt)) freq_gt = freq_vector_simplified[n];


  find_closest_points_in_vector(freq_lt,freq_gt, X.real(), AMI_spec_freq_vector_simplified);

  find_closest_points_in_vector(kx_lt,kx_gt, k[0], AMI_spec_kx_vector_simplified);

  find_closest_points_in_vector(ky_lt,ky_gt, k[1], AMI_spec_ky_vector_simplified);
  //find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec)
  int lll_corner;
  int llg_corner;
  int lgl_corner;
  int lgg_corner;
  int gll_corner;
  int glg_corner;
  int ggl_corner;
  int ggg_corner;
  for (int n=0;n<AMI_spec_kx_vector.size();n++){
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_lt) lll_corner=n;
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_lt) llg_corner=n;
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_gt) lgl_corner=n;
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_gt) lgg_corner=n;
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_lt) gll_corner=n;
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_lt) glg_corner=n;
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_gt) ggl_corner=n;
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_gt) ggg_corner=n;

  }
  /*std::cout << "lll: "<<AMI_spec_freq_vector[lll_corner]<<"  "<<AMI_spec_kx_vector[lll_corner]<<"  "<<AMI_spec_ky_vector[lll_corner] << "  Real: " << AMI_spec_se_Re_vector[lll_corner] << "  Imag: " << AMI_spec_se_Im_vector[lll_corner] <<std::endl;
  std::cout << "ggg: "<<AMI_spec_freq_vector[ggg_corner]<<"  "<<AMI_spec_kx_vector[ggg_corner]<<"  "<<AMI_spec_ky_vector[ggg_corner] << "  Real: " << AMI_spec_se_Re_vector[ggg_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggg_corner] <<std::endl;
	std::cout << "llg: "<<AMI_spec_freq_vector[llg_corner]<<"  "<<AMI_spec_kx_vector[llg_corner]<<"  "<<AMI_spec_ky_vector[llg_corner] << "  Real: " << AMI_spec_se_Re_vector[llg_corner] << "  Imag: " << AMI_spec_se_Im_vector[llg_corner] <<std::endl;
	std::cout << "ggl: "<<AMI_spec_freq_vector[ggl_corner]<<"  "<<AMI_spec_kx_vector[ggl_corner]<<"  "<<AMI_spec_ky_vector[ggl_corner] << "  Real: " << AMI_spec_se_Re_vector[ggl_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggl_corner] <<std::endl;
	std::cout << "lgl: "<<AMI_spec_freq_vector[lgl_corner]<<"  "<<AMI_spec_kx_vector[lgl_corner]<<"  "<<AMI_spec_ky_vector[lgl_corner] << "  Real: " << AMI_spec_se_Re_vector[lgl_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgl_corner] <<std::endl;
  std::cout << "gll: "<<AMI_spec_freq_vector[gll_corner]<<"  "<<AMI_spec_kx_vector[gll_corner]<<"  "<<AMI_spec_ky_vector[gll_corner] << "  Real: " << AMI_spec_se_Re_vector[gll_corner] << "  Imag: " << AMI_spec_se_Im_vector[gll_corner] <<std::endl;
	std::cout << "glg: "<<AMI_spec_freq_vector[glg_corner]<<"  "<<AMI_spec_kx_vector[glg_corner]<<"  "<<AMI_spec_ky_vector[glg_corner] << "  Real: " << AMI_spec_se_Re_vector[glg_corner] << "  Imag: " << AMI_spec_se_Im_vector[glg_corner] <<std::endl;
	std::cout << "lgg: "<<AMI_spec_freq_vector[lgg_corner]<<"  "<<AMI_spec_kx_vector[lgg_corner]<<"  "<<AMI_spec_ky_vector[lgg_corner] << "  Real: " << AMI_spec_se_Re_vector[lgg_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgg_corner] <<std::endl;*/

  //linearinterpolation along freq axis
  double se_Im_ll_corner = AMI_spec_se_Im_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Im_vector[llg_corner] - AMI_spec_se_Im_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
  double se_Im_lg_corner = AMI_spec_se_Im_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Im_vector[lgg_corner] - AMI_spec_se_Im_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
  double se_Im_gl_corner = AMI_spec_se_Im_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Im_vector[glg_corner] - AMI_spec_se_Im_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
  double se_Im_gg_corner = AMI_spec_se_Im_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Im_vector[ggg_corner] - AMI_spec_se_Im_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

	double se_Re_ll_corner = AMI_spec_se_Re_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Re_vector[llg_corner] - AMI_spec_se_Re_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
	double se_Re_lg_corner = AMI_spec_se_Re_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Re_vector[lgg_corner] - AMI_spec_se_Re_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
	double se_Re_gl_corner = AMI_spec_se_Re_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Re_vector[glg_corner] - AMI_spec_se_Re_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
	double se_Re_gg_corner = AMI_spec_se_Re_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Re_vector[ggg_corner] - AMI_spec_se_Re_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

  //std::cout<<"\nse_Re_ll_corner:  "<<se_Re_ll_corner<<std::endl;
  //std::cout<<"se_Re_gg_corner:  "<<se_Re_gg_corner<<std::endl;

  //bilinear interpolation in momentum plane
  double bilinear_prefactor=1/((kx_gt-kx_lt)*(ky_gt-ky_lt));
  //std::cout<<"prefactor:  "<<bilinear_prefactor<<std::endl;
  //std::cout<<"se_Re_ll_corner*(kx_gt-k[0])*(ky_gt-k[1]):  "<<se_Re_ll_corner*(kx_gt-k[0])*(ky_gt-k[1])<<std::endl;
  //std::cout << "se_Re_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt):  " << se_Re_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt)<<std::endl;
  double bilinear_main_Re=se_Re_ll_corner*(kx_gt-k[0])*(ky_gt-k[1]) + se_Re_lg_corner*(kx_gt-k[0])*(k[1]-ky_lt) + se_Re_gl_corner*(k[0]-kx_lt)*(ky_gt-k[1]) + se_Re_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt);
  double  se_Re_interpolated=bilinear_prefactor*bilinear_main_Re;

	double bilinear_main_Im=se_Im_ll_corner*(kx_gt-k[0])*(ky_gt-k[1]) + se_Im_lg_corner*(kx_gt-k[0])*(k[1]-ky_lt) + se_Im_gl_corner*(k[0]-kx_lt)*(ky_gt-k[1]) + se_Im_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt);
  double  se_Im_interpolated=bilinear_prefactor*bilinear_main_Im;

  return std::complex<double> (se_Re_interpolated,se_Im_interpolated);

}

std::vector<double> AmiSpec::return_simple_grid_vector(std::vector<double> &in_vector){
  std::vector<double> out_vector;
  out_vector.push_back(in_vector[0]);
  for (int n=0;n<in_vector.size();n++){
    if (in_vector[n] != out_vector.back()){
      if (std::count(out_vector.begin(),out_vector.end(),in_vector[n] )) return out_vector;
      else out_vector.push_back(in_vector[n]);
    }
  }
}

//std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X){}

void AmiSpec::read_self_energy(std::string file_name){
  std::ifstream MyReadFile(file_name);
  std::string myText;
  std::stringstream ss;
  //std::ifstream MyReadFile(argv[1]);
  //std::cout << argv[1] << std::endl;
  //std::vector<float> frequency;
  //std::vector<float> gf;
  //int count = 0;
  //std::vector<double> se_Im_vector;
  //std::vector<double> se_Im_err_vector;
  //std::vector<double> se_Re_vector;
  //std::vector<double> se_Re_err_vector;
  //std::vector<double> freq_vector;
  //std::vector<double> ky_vector;
  //std::vector<double> kx_vector;


  std::vector<std::string> output_vector;
  std::string a_string;
  getline (MyReadFile, myText);
  while (getline (MyReadFile, myText)) {
    //std::cout <<"myText  " <<myText << std::endl;
    ss.str( myText);

   output_vector.clear();
   //   std::cout << "contents of ss: " << ss.str() << std::endl;
   //std::cout << "a string OUTSIDE: " << a_string << std::endl;
     a_string="";
   while(ss >> a_string){
     //std::cout << "a string: " << a_string << std::endl;
     output_vector.push_back(a_string);
   }
   //std::cout<<"got to here "<< output_vector.back()<<std::endl;
   AMI_spec_kx_vector.push_back(std::stod(output_vector[0]));
   AMI_spec_ky_vector.push_back(std::stod(output_vector[1]));
   AMI_spec_freq_vector.push_back(std::stod(output_vector[2]));
   AMI_spec_se_Re_vector.push_back(std::stod(output_vector[3]));
  AMI_spec_se_Im_vector.push_back(std::stod(output_vector[4]));
   AMI_spec_se_Re_err_vector.push_back(std::stod(output_vector[5]));
   AMI_spec_se_Im_err_vector.push_back(std::stod(output_vector[6]));
   ss.clear();
  }
  AMI_spec_ky_vector_simplified = return_simple_grid_vector(AMI_spec_ky_vector);
  AMI_spec_kx_vector_simplified = return_simple_grid_vector(AMI_spec_kx_vector);
 AMI_spec_freq_vector_simplified = return_simple_grid_vector(AMI_spec_freq_vector);
//for (std::vector<double>::const_iterator i = ky_vector_simplified.begin(); i != ky_vector_simplified.end(); ++i)
//    std::cout << *i << ' ';

// std::cout<<"\n";
// for (std::vector<double>::const_iterator i = freq_vector_simplified.begin(); i != freq_vector_simplified.end(); ++i)    std::cout << *i << ' ';
// std::cout<<"\n";
   //tk::spline spline_I(freq_vector,se_Re_vector);
   //std::cout<< spline_I(std::stod(argv[2]))<<std::endl;
 //std::vector<double> point_wanted={0.0,0.314159,0.0};//0.314159
 //std::cout <<   get_sigma(point_wanted, kx_vector, ky_vector, freq_vector, ky_vector_simplified,  kx_vector_simplified, freq_vector_simplified,se_Im_vector, se_Re_vector) << std::endl;
  MyReadFile.close();

}
