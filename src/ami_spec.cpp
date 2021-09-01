#include "ami_spec.hpp"


AmiSpec::AmiSpec(double xc){

xi_cutoff=xc;

}

std::complex<double> AmiSpec::get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps){
	if(ei.size()!= eps.size()){
	throw std::runtime_error("In A_t the epsilon_t and energy_t do not match in size - exiting");
	}
	std::complex<double> output(0,0);
	for(int i=0; i< ei.size(); i++){

		output+=ei[i]*(double)eps[i];

	}

	return output;
}

std::complex<double> AmiSpec::get_X(X_t &Xsym, xi_t &xi){
	if(Xsym.size()!= xi.size()){
	throw std::runtime_error("In A_t the X_t and xi_t do not match in size - exiting");
	}
	std::complex<double> output(0,0);
	for(int i=0; i< xi.size(); i++){

		output+=(double)Xsym[i]*xi[i];

	}

	return output;
}

std::complex<double> AmiSpec::eval_Aprod(A_prod_t &Ap, xi_t &xi, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

	std::complex<double> output(1,0);

	// A(Sigma, X, E)  : Sigma: self-energy, X: is frequency, E: energy from k_vector
	for(int i=0; i< Ap.size(); i++){

		std::complex<double> this_X=get_X( Ap[i].x_, xi);

		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
		std::complex<double> this_E=eval_tb(1.,0., this_k, mu);

		std::complex<double> this_sigma=get_sigma(this_k, this_X);


		output=output*A_eval(this_sigma, this_X, this_E);



	}


	return output;

}

//TODO: connor please do this part :)
/*std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X){


}*/

void AmiSpec::find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec){
  std::vector<double>  comparison_vector;
  comparison_vector=vec;
  for (int n=0;n<vec.size();n++){
    comparison_vector[n]=std::abs(comparison_vector[n]-point);
  }
  //for(int n=0;n<comparison_vector.size();n++){std::cout << "comparison: " << comparison_vector[n]<<std::endl;}
  auto minimum_dif_iter = std::min_element(comparison_vector.begin(),comparison_vector.end());
  double minimum_dif=*minimum_dif_iter;//gives minimum difference between point and contents of vec
  std::cout<< "\nminimum_dif:  " << minimum_dif << std::endl;
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
if(abs(X)>xi_cutoff) return std::complex<double> (0.0,0.0);
NewAmiCalc::k_vector_t k_copy=k;
for(int for_counter=0;for_counter<k_copy.size();for_counter++){
	if(abs(k_copy[for_counter])>M_PI){
		double modulo_shifted=k_copy[for_counter]+M_PI;
		k_copy[for_counter]=M_PI+modulo_shifted-floor(modulo_shifted/(2*M_PI))*2*M_PI;//fmod(k_copy[for_counter])
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
  std::cout << "lll: "<<AMI_spec_freq_vector[lll_corner]<<"  "<<AMI_spec_kx_vector[lll_corner]<<"  "<<AMI_spec_ky_vector[lll_corner] << "  Real: " << AMI_spec_se_Re_vector[lll_corner] << "  Imag: " << AMI_spec_se_Im_vector[lll_corner] <<std::endl;
  std::cout << "ggg: "<<AMI_spec_freq_vector[ggg_corner]<<"  "<<AMI_spec_kx_vector[ggg_corner]<<"  "<<AMI_spec_ky_vector[ggg_corner] << "  Real: " << AMI_spec_se_Re_vector[ggg_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggg_corner] <<std::endl;
	std::cout << "llg: "<<AMI_spec_freq_vector[llg_corner]<<"  "<<AMI_spec_kx_vector[llg_corner]<<"  "<<AMI_spec_ky_vector[llg_corner] << "  Real: " << AMI_spec_se_Re_vector[llg_corner] << "  Imag: " << AMI_spec_se_Im_vector[llg_corner] <<std::endl;
	std::cout << "ggl: "<<AMI_spec_freq_vector[ggl_corner]<<"  "<<AMI_spec_kx_vector[ggl_corner]<<"  "<<AMI_spec_ky_vector[ggl_corner] << "  Real: " << AMI_spec_se_Re_vector[ggl_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggl_corner] <<std::endl;
	std::cout << "lgl: "<<AMI_spec_freq_vector[lgl_corner]<<"  "<<AMI_spec_kx_vector[lgl_corner]<<"  "<<AMI_spec_ky_vector[lgl_corner] << "  Real: " << AMI_spec_se_Re_vector[lgl_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgl_corner] <<std::endl;
  std::cout << "gll: "<<AMI_spec_freq_vector[gll_corner]<<"  "<<AMI_spec_kx_vector[gll_corner]<<"  "<<AMI_spec_ky_vector[gll_corner] << "  Real: " << AMI_spec_se_Re_vector[gll_corner] << "  Imag: " << AMI_spec_se_Im_vector[gll_corner] <<std::endl;
	std::cout << "glg: "<<AMI_spec_freq_vector[glg_corner]<<"  "<<AMI_spec_kx_vector[glg_corner]<<"  "<<AMI_spec_ky_vector[glg_corner] << "  Real: " << AMI_spec_se_Re_vector[glg_corner] << "  Imag: " << AMI_spec_se_Im_vector[glg_corner] <<std::endl;
	std::cout << "lgg: "<<AMI_spec_freq_vector[lgg_corner]<<"  "<<AMI_spec_kx_vector[lgg_corner]<<"  "<<AMI_spec_ky_vector[lgg_corner] << "  Real: " << AMI_spec_se_Re_vector[lgg_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgg_corner] <<std::endl;

  //linearinterpolation along freq axis
  double se_Im_ll_corner = AMI_spec_se_Im_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Im_vector[llg_corner] - AMI_spec_se_Im_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
  double se_Im_lg_corner = AMI_spec_se_Im_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Im_vector[lgg_corner] - AMI_spec_se_Im_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
  double se_Im_gl_corner = AMI_spec_se_Im_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Im_vector[glg_corner] - AMI_spec_se_Im_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
  double se_Im_gg_corner = AMI_spec_se_Im_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Im_vector[ggg_corner] - AMI_spec_se_Im_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

	double se_Re_ll_corner = AMI_spec_se_Re_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Re_vector[llg_corner] - AMI_spec_se_Re_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
	double se_Re_lg_corner = AMI_spec_se_Re_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Re_vector[lgg_corner] - AMI_spec_se_Re_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
	double se_Re_gl_corner = AMI_spec_se_Re_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Re_vector[glg_corner] - AMI_spec_se_Re_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
	double se_Re_gg_corner = AMI_spec_se_Re_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Re_vector[ggg_corner] - AMI_spec_se_Re_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);



  //bilinear interpolation in momentum plane
  double bilinear_prefactor=1/((kx_gt-kx_lt)*(ky_gt-ky_lt));
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



// todo: probably don't need to pass A if this function takes in X and E already
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;

}

void AmiSpec::generate_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0){

// first lets identify which G's have external frequencies -
AmiBase::g_prod_t innert;
AmiBase::g_prod_t nert;

// std::vector<int> pole_indices;
// ami_term.g_list
for(int i=0; i<start_term.g_list.size(); i++){

	if(start_term.g_list[i].alpha_.back()!=0){
		nert.push_back(start_term.g_list[i]);
	}else{
		innert.push_back(start_term.g_list[i]);
	}

}

// Next we want to create the product of ( a +d)*(b+d)*(c+d)... So we want to specify how many deltas are in each one
int max_deltas=nert.size();

std::vector<std::vector<int>> pp_v;
ami.get_pp_comb(max_deltas,pp_v);

sp_terms new_terms;
new_terms.resize(pp_v.size()); // one new term for each combination in pp_v

for(int i=0; i< pp_v.size(); i++){
	std::cout<<"Term:";
	for(int m=0; m< pp_v[i].size(); m++){
	std::cout<<pp_v[i][m]<<" ";

	}
	std::cout<<std::endl;
}

for(int i=0; i< pp_v.size(); i++){
	for(int m=0; m< pp_v[i].size(); m++){

	if(pp_v[i][m]==0){

		new_terms[i].dprod_.push_back(nert[m]);
		new_terms[i].delta_count++;
	}

	if(pp_v[i][m]==1){

	new_terms[i].ami_term_.g_list.push_back(nert[m]);
	}

	}

}

// Generate the A product from R0

A_prod_t this_Ap;
R0_to_Aprod(R0, this_Ap);

// std::cout<<"This ap has size "<<this_Ap<<std::endl;

// now attach the innert G's to each new term

for(int i=0; i< new_terms.size(); i++){

new_terms[i].ami_term_.g_list.insert(new_terms[i].ami_term_.g_list.end(), innert.begin(), innert.end());
new_terms[i].ami_term_.p_list= start_term.p_list;
new_terms[i].aprod_=this_Ap;

}

// no update to the sign needed - the delta_count will tell about prefactors when evaluating
// each ami_term.p_list is unchanged from the original
// What is missing is how to handle the A_prod_t - which at least at the start is defined by R0. and the same for every term.



new_sp_terms=new_terms;



}


void AmiSpec::resolve_deltas( sp_terms &sp_terms){

for(int i=0; i< sp_terms.size(); i++){
resolve_deltas(sp_terms[i]);
}

}


// TODO: if the prefactor is not +1 or -1, then need to modify the sign for the term to account for the delta rule delta(2*x)=delta(x)/|2|

// TODO: if there is a delta function that is squared, then it is the same as the single delta function but still gives the volume element due to normalization in the larger space.
void AmiSpec::resolve_deltas(ami_sp_term &sp_term){


std::cout<<"Number of deltas is "<<sp_term.dprod_.size()<<std::endl;

if(	sp_term.dprod_.size()==0){ return;}

// look at each delta. create a pole for each, and assign an index_ to it that represents which xi will be replaced.  do this for each and make sure each delta gets a unique xi to replace
// for each pole with respect to its xi values generate the actual pole
// for each pole replace the xi in every Aterm, g_prod, fermi_pole AND other remaining deltas. once done mark the delta for removal
// for removal define an empty integer vec of size delta.size(). initialize to zero. set to 1 for each resolved delta.  check at the end that all the entries are 1.  and then remove all of the deltas. otherwise throw an error.  OR remove only the entries that are resolved.

AmiBase::pole_array_t pv;
pv.resize(sp_term.dprod_.size());

int xi_size=sp_term.dprod_[0].eps_.size();
int delta_size=sp_term.dprod_.size();

std::vector<int> used, assigned;
used.resize(xi_size,0);
assigned.resize(delta_size,0);

for(int i=0; i< xi_size; i++){

// search through deltas to find a non-zero and then break

	for(int m=0; m< sp_term.dprod_.size(); m++){
		if( assigned[m]==0){
		if(sp_term.dprod_[m].eps_[i]!=0 ){
			pv[m].index_=i;
			used[i]=1;
			assigned[m]=1;
			break;

		}
		}

	}

if( count(used.begin(), used.end(),1)==delta_size){
	break;
}


}

if(count(used.begin(), used.end(),1)!= delta_size){
	throw std::runtime_error("Didn't resolve deltas correctly - exiting");
	exit(0);
}

// at this point we should have a list of xi we are going to replace. so generate the poles

 for(int i=0; i< sp_term.dprod_.size(); i++){

	pv[i].eps_.resize( sp_term.dprod_[i].eps_.size());

	pv[i].alpha_.resize( sp_term.dprod_[i].alpha_.size());

	int this_index=pv[i].index_;
	int this_prefactor=sp_term.dprod_[i].eps_[this_index];


	for( int m=0; m< pv[i].alpha_.size(); m++){
		pv[i].alpha_[m]=sp_term.dprod_[i].alpha_[m]*(-this_prefactor);
	}

	for(int m=0; m< pv[i].eps_.size(); m++){
		if(m==this_index){pv[i].eps_[m]=0;}
		else{

			pv[i].eps_[m]= sp_term.dprod_[i].eps_[m]*(-this_prefactor);
		}

	}


}

// now at this point we should have all of the poles accounted for
// std::cout<<std::count(used.begin(), used.end(),1)<<std::endl;

std::cout<<"BEFORE: Printing poles "<<std::endl;

for(int i=0; i< pv.size(); i++){

std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
print_int_vec(pv[i].alpha_);
std::cout<<" | eps-";
print_int_vec(pv[i].eps_);
std::cout<<std::endl;


}


// so now replace poles
for(int i=0; i< pv.size(); i++){

replace_xi(i,pv, sp_term);


}

std::cout<<"AFTER: Printing poles "<<std::endl;

for(int i=0; i< pv.size(); i++){

std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
print_int_vec(pv[i].alpha_);
std::cout<<" | eps-";
print_int_vec(pv[i].eps_);
std::cout<<std::endl;


}


// we should have a catch in case something goes wrong.
sp_term.dprod_.clear();



// for(int i=0; i< sp_term.dprod_.size(); i++){
	// pv[i].eps_=sp_term.dprod_[i].eps_;
	// pv[i].alpha_=sp_term.dprod_[i].alpha_;
// }








}

void AmiSpec::replace_xi(int i, AmiBase::pole_array_t &pv, ami_sp_term &sp_term){


AmiBase::pole_struct this_pole=pv[i];

// update other poles j>i
for(int j=0; j< pv.size(); j++){
	if(j==i){ continue;}
	update_spec_pole(pv[i], pv[j].alpha_, pv[j].eps_);

}

// update the fermi poles in the sp_term

for(int j=0; j< sp_term.ami_term_.p_list.size(); j++){

	update_spec_pole(pv[i],sp_term.ami_term_.p_list[j].alpha_,sp_term.ami_term_.p_list[j].eps_);

}

// update the g_prod
for(int j=0; j< sp_term.ami_term_.g_list.size(); j++){

	update_spec_pole(pv[i], sp_term.ami_term_.g_list[j].alpha_, sp_term.ami_term_.g_list[j].eps_);

}

// update the Aprod
// note that Aprod is defined by the X_t x_ values not eps_
for(int j=0; j< sp_term.aprod_.size(); j++){

	update_spec_pole(pv[i], sp_term.aprod_[j].alpha_, sp_term.aprod_[j].x_);

}





}

void AmiSpec::update_spec_pole(AmiBase::pole_struct &source_pole, AmiBase::alpha_t &target_alpha, AmiBase::epsilon_t &target_eps){

// std::cout<<"Updating : starting with: alpha - ";
// print_int_vec(target_alpha);
// std::cout<<" | and eps=";
// print_int_vec(target_eps);
// std::cout<<std::endl;

// std::cout<<"USING pole: alpha - ";
// print_int_vec(source_pole.alpha_);
// std::cout<<" | and eps=";
// print_int_vec(source_pole.eps_);
// std::cout<<std::endl;


int index=source_pole.index_;
if(target_eps[index]==0){ return;} // nothing to do

for(int m=0; m< target_alpha.size(); m++){
	target_alpha[m]+=source_pole.alpha_[m];
}

for(int m=0; m< target_eps.size(); m++){
	if(m==index){ target_eps[m]=0;}
	else{
	target_eps[m]+=source_pole.eps_[m];
	}
}




// std::cout<<"Updating : resulting in: alpha - ";
// print_int_vec(target_alpha);
// std::cout<<" | and eps=";
// print_int_vec(target_eps);
// std::cout<<std::endl;


	return;


}



// epsilon_i either by index or by value
// int eps_index=-1;
// std::complex<double> eps_val;

// AmiBase::epsilon_t eps_;
// X_t x_;
// AmiBase::alpha_t alpha_;
// AmiBase::species_t species_;


// };
void AmiSpec::R0_to_Aprod(AmiBase::g_prod_t &R0, A_prod_t &Ap){

// R0 is g_prod_t
// each g_prod_t has an epsilon and an alpha

Ap.clear();
Ap.resize(R0.size());

for(int i=0; i< R0.size(); i++){

X_t this_X;
this_X.resize(R0.size(),0);
this_X[i]=1;

Ap[i].alpha_=R0[i].alpha_;
Ap[i].species_=R0[i].species_;
Ap[i].x_=this_X;
Ap[i].eps_index=i;
// Ap[i].eps_=R0[i].eps_; // not sure this is needed

}


}

std::complex<double> AmiSpec::eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu){

	std::complex<double> output(0,0);

	for(int i=0; i<k.size();i++){

	output+=-2.0*t*cos(k[i]);

	}

	// uncomment this block if want t' in dispersion
	// double term=-4.0*tp;
	// for(int i=0; i<k.size(); i++){

		// term=term*cos(k[i]);
	// }
	// output+=term;


	output -= mu;

	return output;

}


std::complex<double> AmiSpec::construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

std::complex<double> result=0;

NewAmiCalc::k_vector_t this_k=ami.construct_k(alpha, klist);
result=eval_tb(1.,0., this_k, mu);


return result;

}



// void generate_sp_terms(ami_term &start_term, sp_terms &new_sp_terms);
// void reduce_deltas(ami_sp_term &term);
