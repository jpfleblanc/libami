#include "ami_base.hpp"

// TODO: since there is some duplication of code in these two evaluate commands
// they should be separated
std::complex<double> AmiBase::evaluate(ami_parms &parms, R_t &R_array,
                                       P_t &P_array, S_t &S_array,
                                       ami_vars &external, g_prod_t &unique_g,
                                       R_ref_t &Rref, ref_eval_t &Eval_list) {
  overflow_detected = false;

  if (Rref.size() == 0 || unique_g.size() == 0 || Eval_list.size() == 0) {
    return evaluate(parms, R_array, P_array, S_array, external);
  }

  int dim = parms.N_INT_;

  SorF_t SorF_result;

  if (dim == 0) {
    std::complex<double> gprod;

    gprod = eval_gprod(parms, R_array[0][0], external);

    return gprod;
  }

  if (dim == 1) {
    SorF_t SF_left, SF_right;
    SorF_t S_double_left, S_double_right;

    SF_left = dot(S_array[0], fermi(parms, P_array[0], external));

    SorF_result = SF_left;

  }

  for (int i = 0; i < dim - 1; i++) {
    SorF_t SF_left, SF_right;
    SorF_t S_double_left, S_double_right;


    if (i == 0) {
      // do dot operation
      SF_left = dot(S_array[i], fermi(parms, P_array[i], external));
    } else {
      SF_left = SorF_result;
    }

    // do dot
    SF_right = dot(S_array[i + 1], fermi(parms, P_array[i + 1], external));

    

    SorF_result = cross(SF_left, SF_right);
    
  }

  std::complex<double> final_result;

  final_result =
      optimized_star(parms, SorF_result, unique_g, Rref, Eval_list, external);


  return final_result;
}

std::complex<double> AmiBase::optimized_star(ami_parms &parms, SorF_t K,
                                             g_prod_t &unique_g, R_ref_t &Rref,
                                             ref_eval_t &Eval_list,
                                             ami_vars external) {
  

#ifdef BOOST_MP
  std::complex<boost::multiprecision::float128> term;
  std::complex<boost::multiprecision::float128> gprod;
  std::complex<boost::multiprecision::float128> output(0, 0);
  boost::multiprecision::float128 prefactor = external.prefactor;
  
  std::vector<std::complex<boost::multiprecision::float128>> unique_vals;
#else
  std::complex<double> term;
  std::complex<double> output = 0;
  double prefactor = external.prefactor;
  std::vector<std::complex<double>> unique_vals;
  std::complex<double> gprod;
#endif
  

  

  for (int i = 0; i < unique_g.size(); i++) {
    // Context: we pretend each G is a g_prod and use the same function as the normal
    // star operation for evaluating products of G's
    g_prod_t this_term;
    this_term.push_back(unique_g[i]);

#ifdef BOOST_MP
	gprod =
        eval_gprod_mp(parms, this_term, external) * prefactor;


	if ((boost::multiprecision::abs(std::real(gprod)) > precision_cutoff) ||
        (boost::multiprecision::abs(std::imag(gprod)) > precision_cutoff)) {
      overflow_detected = true;
    }
#else
	gprod =
			eval_gprod(parms, this_term, external) * prefactor;
	

    if ((std::abs(std::real(gprod)) > precision_cutoff) ||
        (std::abs(std::imag(gprod)) > precision_cutoff)) {
      overflow_detected = true;
    }
#endif

    unique_vals.push_back(gprod); // This removes the overall prefactor for
                                  // each g. we then add that back later
  }


  // now we have the numerical values of our unique G's

  for (int i = 0; i < Eval_list.size(); i++) {
#ifdef BOOST_MP
    std::complex<boost::multiprecision::float128> ksum(0, 0);
#else
    std::complex<double> ksum(0, 0);
#endif
    for (int j = 0; j < Eval_list[i].size(); j++) {
 

#ifdef BOOST_MP
      std::complex<boost::multiprecision::float128> this_k(
          K[0][Eval_list[i][j].first].real(),
          K[0][Eval_list[i][j].first].imag());
      boost::multiprecision::float128 factor = double(Eval_list[i][j].second);

#else
      std::complex<double> this_k = K[0][Eval_list[i][j].first];
      double factor = double(Eval_list[i][j].second);
#endif

      ksum += this_k * factor;
    }
    // in principle, every entry in eval_list has the same Rref terms
    ref_v_t pair_vec = Rref[i]; // just grab the first one

#ifdef BOOST_MP
    std::complex<boost::multiprecision::float128> this_gprod(1, 0);
    for (int j = 0; j < pair_vec.size(); j++) {
      
      std::complex<boost::multiprecision::float128> this_val =
          unique_vals[pair_vec[j].first];
      this_gprod = this_gprod * this_val;
	  
	  if(verbose){
		std::cout<<j<<" "<< this_val<<" "<< this_gprod<<std::endl;  
		  }
	  
    }
    std::complex<boost::multiprecision::float128> ksum_mp(ksum.real(),
                                                          ksum.imag());
    term = ksum_mp * this_gprod * prefactor;
	
	
    if ((boost::multiprecision::abs(std::real(term)) > precision_cutoff) ||
        (boost::multiprecision::abs(std::imag(term)) > precision_cutoff)) {
      overflow_detected = true;
    }
#else

    std::complex<double> this_gprod(1, 0);

    
    for (int j = 0; j < pair_vec.size(); j++) {
      
      this_gprod = this_gprod * unique_vals[pair_vec[j].first];
    }

    
    term = ksum * this_gprod *
           prefactor; // add back the overall prefactor for this term

    if ((std::abs(std::real(term)) > precision_cutoff) ||
        (std::abs(std::imag(term)) > precision_cutoff)) {
      overflow_detected = true;
    }

#endif

  
    output += term;
  }

 
#ifdef BOOST_MP
  std::complex<double> final_output(output.real().convert_to<double>(),
                                    output.imag().convert_to<double>());


#else
  std::complex<double> final_output = output;
#endif

 
  return final_output;
}

/**
 * This is a primary AMI symbolic evaluation function.  It takes a solution
 * defined by S, P, and R arrays and the external parameters and returns a
 * `complex<double>` restult.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] R_array: Input `R_t`.
 * @param[in] P_array : Input `P_t`.
 * @param[in] S_array : Input `S_t`.
 * @param[in] external : Input `ami_vars` containing point to evaluate.
 * @return Result is returned as single value of type `std::complex<double>`.
 */
std::complex<double> AmiBase::evaluate(ami_parms &parms, R_t &R_array,
                                       P_t &P_array, S_t &S_array,
                                       ami_vars &external) {
  overflow_detected = false;

  int dim = parms.N_INT_;

  SorF_t SorF_result;

  if (dim == 0) {
    std::complex<double> gprod;

    gprod = eval_gprod(parms, R_array[0][0], external);

    return gprod;
  }

  if (dim == 1) {
    SorF_t SF_left, SF_right;
    SorF_t S_double_left, S_double_right;

    SF_left = dot(S_array[0], fermi(parms, P_array[0], external));

    SorF_result = SF_left;
  }

  for (int i = 0; i < dim - 1; i++) {
    // std::cout<<"i is "<<i<<" of dim:"<<dim<<std::endl;
    SorF_t SF_left, SF_right;
    SorF_t S_double_left, S_double_right;

    if (i == 0) {
      // do dot operation
      SF_left = dot(S_array[i], fermi(parms, P_array[i], external));
 
    } else {
      SF_left = SorF_result;
    }

    // do dot

    SF_right = dot(S_array[i + 1], fermi(parms, P_array[i + 1], external));



    SorF_result = cross(SF_left, SF_right);
  
  }

  std::complex<double> final_result;

  
  final_result = star(parms, SorF_result, R_array[dim], external);

  return final_result;
}

std::complex<double> AmiBase::star(ami_parms &parms, SorF_t K, Ri_t R,
                                   ami_vars external) {
  std::complex<double> output = 0;
  std::complex<double> term;
  std::complex<double> gprod;

  std::complex<double> t0, t10;

  // These are necessary debugging lines 
  // std::ofstream file;
  // file.open("outfile.dat",  std::ofstream::out | std::ofstream::app);
  bool print_output = false;

  for (int i = 0; i < K[0].size(); i++) {
    print_output = false;

    gprod = eval_gprod(parms, R[i], external);
    term = K[0][i] * gprod;

    // if overflow flag not set then check for overflow
    if (!overflow_detected) {
      if (std::abs(std::real(K[0][i])) > precision_cutoff) {
        overflow_detected = true;
      }

      if ((std::abs(std::real(term)) > precision_cutoff) ||
          (std::abs(std::imag(term)) > precision_cutoff)) {
        overflow_detected = true;
      } else {
        // edit -> G can be exactly equal to 1 without indicating an
        // overflow
        if ((std::floor(std::abs(std::real(gprod))) ==
             std::abs(std::real(gprod))) &&
            std::abs(std::real(gprod)) != 1) {
         
          overflow_detected = true;
        }
      }
    }

	// These are necessary debugging lines that should stay 
    // std::cout<<"In star K[]*R"<<std::endl;
    // if(std::abs( std::real(term))>1000){
    // print_output=true;
    // std::cout<< std::setprecision(20)<< i<<" "<< K[0][i] <<" "<<
    // std::real(gprod)<<" "<<std::imag(gprod)<< " "<<std::real(term)<<" "<<
    // std::imag(term) <<" CO="<<output <<std::endl;


    output += term;

 
  }

  // file.close();


  return output;
}

/// Evaluation dot-product from AMI paper.
AmiBase::SorF_t AmiBase::dot(Si_t Si, SorF_t fermi) {
  SorF_t output;
  output.reserve(Si.size());

  for (int i = 0; i < Si.size(); i++) {
    std::vector<std::complex<double>> line;
    line.reserve(Si[i].size());
    for (int j = 0; j < Si[i].size(); j++) {
      line.push_back(Si[i][j] * fermi[i][j]);
    }

    output.push_back(line);
  }
  
  // std::cout<<"Returning dot output of size "<<output.size()<<std::endl;

  return output;
}

/// Evaluation cross operator from AMI paper.
AmiBase::SorF_t AmiBase::cross(SorF_t left, SorF_t right) {
  SorF_t output;

  std::vector<std::complex<double>> line;
  
  line.reserve(left[0].size() *
               right.back().size()); // just guess that the last one is
                                     // the biggest for reservation

  for (int i = 0; i < left[0].size(); i++) {
    for (int rj = 0; rj < right[i].size(); rj++) {
      line.push_back(left[0][i] * right[i][rj]);
    }
  }

  output.push_back(line);
  return output;
}

/// The fermi/bose operator. Translates a vector of poles, into a vector of
/// Fermi/Bose functions (and their derivatives) with specific numerical values
/// dependent upon input `ami_vars`.
AmiBase::SorF_t AmiBase::fermi(ami_parms &parms, Pi_t Pi, ami_vars external) {
  SorF_t output;
  output.reserve(Pi.size());

  for (int i = 0; i < Pi.size(); i++) {
    std::vector<std::complex<double>> group;
    group.reserve(Pi[i].size());

    for (int j = 0; j < Pi[i].size(); j++) {

      group.push_back(fermi_pole(parms, Pi[i][j], external));
    }

    output.push_back(group);
  }

  return output;
}

/// Evaluation of a single pole in Fermi/Bose functions for a given `ami_vars`.
std::complex<double> AmiBase::fermi_pole(ami_parms &parms, pole_struct pole,
                                         ami_vars external) {
  if (verbose) {
    std::cout << "Working on pole" << std::endl;
    print_pole_struct_info(pole);
  }

  std::complex<double> output, term;
  int eta = 0;

  double beta = external.BETA_;
  double E_REG = parms.E_REG_;

  // Spectral evaluation only
  std::complex<double> freq_shift(0, 0);
  if (pole.x_alpha_.size() != 0) {
    for (int i = 0; i < pole.x_alpha_.size(); i++) {
      freq_shift = external.frequency_[i] * (double)pole.x_alpha_[i];
    }
  }

  // Future Development: 
  // In order to generalize to have fermi and bose lines, from here to 'sigma' needs
  // to be considered.

  // example fixed
  //
  // 1) create a stat map for the frequencies std::vector<int> stat_map:
  // populate 1 for fermi and 0 for bose. length is equal to alpha.size() 2)
  // simply replace eta=eta + 1*stat_map[i]
  //

  // alternate fix.  parms.TYPE_ is 0 for sigma, 1 for Pi etc.  So if
  // parms.TYPE_==1 and pole.alpha_.back()==1 (or -1), don't add one. else add
  // one to eta.


  for (int i = 0; i < pole.alpha_.size() - 1; i++) {
   
    if (pole.alpha_[i] != 0) {
      eta++;
    }
  }

  // handle external based on graph type
  if (pole.alpha_.back() != 0 && parms.TYPE_ != AmiBase::Pi_phuu &&
      parms.TYPE_ != AmiBase::Pi_phud && parms.TYPE_ != AmiBase::doubleocc &&
      parms.TYPE_ != AmiBase::Pi_ppuu && parms.TYPE_ != AmiBase::Pi_ppud &&
      parms.TYPE_ != AmiBase::FORCE) {
    eta++;
   
  }

  // if this is a double occupancy graph then the external line is bosonic. so
  // it is a bosonic matsubara integral. so eta needs to be incremented IF the
  // pole is for the last integration index
  if (parms.TYPE_ == AmiBase::doubleocc) {
    if (pole.index_ == pole.alpha_.size() - 1) {
      eta++;
    
    }
  }

  // END TODO

  // could put infor into ami_vars external as to what the state type of the
  // external variables is.
  std::complex<double> E = get_energy_from_pole(pole, external);

  // In the case of spectral poles the freq_shift might not be zero
  E = E + freq_shift;


  double sigma = pow(-1.0, double(eta));

  std::complex<double> zero(0, 0);
  std::complex<double> im(0, 1);

  // If energy denominator would be zero attempts to regulate if bosonic
  // function
  if (E == zero && sigma == -1) { // && pole.der_==0    Not sure if pole
    // derivative plays any role

    if (drop_bosonic_diverge) {

      return zero; // This is a dangerous approximation. 
    }
    E += E_REG;
  } else {
    if (sgn(E.real()) != 0) {
      E += E_REG * sgn(E.real());
    } else {
      E += E_REG;
    }
  }

  if (drop_der && pole.der_ != 0) {
    return zero;
  }

  int m = pole.der_;

  // compute m'th derivative
  output = fermi_bose(m, sigma, beta, E);

  if (verbose) {
    std::cout << "Fermi Pole Term evaluated to " << output << " at energy " << E
              << " with sigma " << sigma << " betaE is " << beta * E
              << " in exponent " << std::exp(beta * (E)) << std::endl;
			  
	std::cout<<"Energy list is :(";
		for(int ii=0; ii< external.energy_.size(); ii++){
			std::cout<<std::setprecision(20)<<" "<<external.energy_[ii]<<" ,";
		}
		std::cout<<")"<<std::endl;
  }

  if (parms.TYPE_ == AmiBase::doubleocc) {
    output = -1.0 * output;
  }

  return output;
}

/// This computes the mth order derivative of the Fermi function or the
/// negative of the Bose distribution functions given by \f$\frac{1}{\sigma
/// \exp^{\beta E}+1} \f$ at \f$ \beta\f$, for energy \f$ E\f$. \f$
/// \sigma=1.0\f$ for Fermi and -1.0 for Bose.
std::complex<double> AmiBase::fermi_bose(int m, double sigma, double beta,
                                         std::complex<double> E) {
  std::complex<double> output, term;
  output = 0.0;

  if (m == 0) {
    double arg = std::real(beta * E);
    double arg_amp = std::abs(arg);

    if (arg > exp_max_arg) {
      double arg_sign = (double)sgn(arg);
      if (arg_sign > 0) {
        output = 0;
      } else {
        output = 1;
      }
    } else {
      output = 1.0 / (sigma * std::exp(beta * (E)) + 1.0);
    }
    // return output;
  } else { // compute m'th derivative

    for (int k = 0; k < m + 1; k++) {
      // depricated: Original format
      // term = frk(m, k) * std::exp(k * beta * (E)) * std::pow(sigma, k) *
             // std::pow(-1.0, k + 1) /
             // std::pow(sigma * std::exp(beta * (E)) + 1.0, k + 1);
	// Reformatted with fewer exponentials
      term= frk(m,k)*std::pow(sigma,k) *std::pow(-1.0,
      k+1)*(1.0/(sigma*std::exp(beta*(E))+1.0)/std::pow(sigma
      +std::exp(-beta*(E)), k)) ;
      output += term;

      if (verbose) {
        std::cout << "On k'th derivative " << k << std::endl;

        std::cout << "Fermi-bose Term evaluated to " << term << " at energy "
                  << E << " with sigma " << sigma << " betaE is " << beta * E
                  << " in exponent " << std::exp(beta * (E)) << std::endl;
      }
    }

    output = output * std::pow(beta, m) * (-1.0);

    if ((std::abs(std::real(output)) > precision_cutoff)|| (std::abs(std::imag(output)) > precision_cutoff)) {
    
      overflow_detected = true;
    }
  }

  return output;
}

std::complex<double> AmiBase::get_energy_from_pole(pole_struct pole,
                                                   ami_vars external) {
  std::complex<double> output(0, 0);

  // Evaluating energies for pole
  for (int i = 0; i < pole.eps_.size(); i++) {
    output += double(pole.eps_[i]) * external.energy_[i];
  }

  return output;
}

std::complex<double> AmiBase::get_energy_from_g(g_struct g, ami_vars external) {
  std::complex<double> output = 0;

  for (int i = 0; i < g.eps_.size(); i++) {
    output += double(g.eps_[i]) * external.energy_[i];
  }

  return output;
}

/**
 *
 * Numerical evaluation of a product of Green's functions. Used both in `terms`
 * and `R_t` constructions.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] g_prod : `g_prod_t` a list of `g_struct` that is interpretted as
 * \f$ \prod{G_i} \f$.
 * @param[in] external : Input external variables in a `ami_vars` struct.
 * @return Single value for product of Green's functions.
 */
std::complex<double> AmiBase::eval_gprod(ami_parms &parms, g_prod_t g_prod,
                                         ami_vars external) {
  std::complex<double> output(0, 0);

  std::complex<double> denom_prod(1, 0);
  double prefactor = external.prefactor;


  double E_REG = parms.E_REG_;

  bool verbose = false;


  for (int i = 0; i < g_prod.size(); i++) {
    std::complex<double> alphadenom(0, 0);
    std::complex<double> epsdenom(0, 0);

    for (int a = 0; a < g_prod[i].alpha_.size(); a++) {
      alphadenom += double(g_prod[i].alpha_[a]) * external.frequency_[a];
    }

    std::complex<double> zero(0, 0);
    std::complex<double> im(0, 1);

    for (int a = 0; a < g_prod[i].eps_.size(); a++) {
      epsdenom += double(g_prod[i].eps_[a]) * external.energy_[a];

    }

  
    denom_prod = denom_prod * (alphadenom + epsdenom);
   
  }

  output = 1.0 / denom_prod * prefactor;

  return output;
}

#ifdef BOOST_MP

std::complex<boost::multiprecision::float128>
AmiBase::eval_gprod_mp(ami_parms &parms, g_prod_t g_prod, ami_vars external) {
  std::complex<boost::multiprecision::float128> output(0, 0);

  std::complex<boost::multiprecision::float128> denom_prod(1, 0);
  boost::multiprecision::float128 prefactor = external.prefactor;
  

  double E_REG = parms.E_REG_;

  bool verbose = false;


  for (int i = 0; i < g_prod.size(); i++) {
    std::complex<boost::multiprecision::float128> alphadenom(0, 0);
    std::complex<boost::multiprecision::float128> epsdenom(0, 0);

    for (int a = 0; a < g_prod[i].alpha_.size(); a++) {
      alphadenom += double(g_prod[i].alpha_[a]) * external.frequency_[a];
    }

    std::complex<double> zero(0, 0);
    std::complex<double> im(0, 1);
    for (int a = 0; a < g_prod[i].eps_.size(); a++) {
      epsdenom += double(g_prod[i].eps_[a]) * external.energy_[a];
    }


    std::complex<boost::multiprecision::float128> alphadenom_mp(
        alphadenom.real(), alphadenom.imag());
    std::complex<boost::multiprecision::float128> epsdenom_mp(epsdenom.real(),
                                                              epsdenom.imag());

    denom_prod = denom_prod * (alphadenom_mp + epsdenom_mp);
    
  }

  std::complex<boost::multiprecision::float128> one(1, 0);

  output = one / denom_prod * prefactor;

  return output;
}

#endif

/// Using notation to match https://doi.org/10.1103/PhysRevB.99.035120.
/// They produced coefficients to the fermi functions and put them in a table.
/// We derive a general expression for those coefficients - we believe this to
/// be general but have only checked up to 6th order.
double AmiBase::frk(int r, int k) {
  double output = 0.0;

  for (int m = 0; m < k + 1; m++) {
    output += binomialCoeff(k, m) * std::pow(m, r) * (std::pow(-1, k - m));
  }

  
  return output;
}

/// Returns value of Binomial Coefficient C(n, k).
int AmiBase::binomialCoeff(int n, int k) {
  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if (k > n - k)
    k = n - k;

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

void AmiBase::get_triggers(ami_parms &parms, R_t &R_array,P_t &P_array, S_t &S_array,ami_vars &external, std::vector< std::vector<int>> &triggers){
triggers.clear(); 


for(int i=0; i<R_array.size(); i++){
  for(int j=0; j<R_array[i].size(); j++){
    for(int m=0; m<R_array[i][j].size(); m++){
  
  std::complex<double> epsdenom(0, 0);
   std::complex<double> zero(0, 0);
   bool not_symbolic_zero=false;
  for (int a = 0; a < R_array[i][j][m].eps_.size(); a++) {
      epsdenom += double(R_array[i][j][m].eps_[a]) * external.energy_[a];

    if(R_array[i][j][m].eps_[a]!=0){
      not_symbolic_zero=true;
    }

    }
    // std::cout<<"On i,j: "<<i<<" "<<j<<" "<<m<<" eps denom gave "<<epsdenom<<std::endl;
  
    if(epsdenom==zero && not_symbolic_zero){
      std::vector<int> issue={i,j,m};
      triggers.push_back(issue);
      
      
      
    }
  
    // std::cout<<"On i,j: "<<i<<" "<<j<<"Evaluation gave "<<eval_gprod(parms, R_array[i][j], external)<<std::endl;
    }
    
    
    
  }
  
  if(triggers.size()!=0){
    break;
  }
  
}
 


  
  
  
}



void AmiBase::find_equal_values(ami_parms &parms, R_t &R_array, ami_vars &external, std::vector<std::vector<int>> &equal_pairs){
equal_pairs.clear(); 

double tol=parms.tol_;

for (int ord=0; ord<R_array.size(); ord++){
  for(int num=0; num<R_array[ord].size(); num++){
    for(int m=0; m<R_array[ord][num].size(); m++){
      for(int n=m+1; n<R_array[ord][num].size(); n++){
        
      if(!isEqual(R_array[ord][num][m].eps_,R_array[ord][num][n].eps_)){
        std::complex<double> eps_m=0.0;
        std::complex<double> eps_n=0.0;
        // if G's are not identical already
        for (int a = 0; a < R_array[ord][num][m].eps_.size(); a++) {
      eps_m += double(R_array[ord][num][m].eps_[a]) * external.energy_[a];
        }
        
        for (int a = 0; a < R_array[ord][num][n].eps_.size(); a++) {
      eps_n += double(R_array[ord][num][n].eps_[a]) * external.energy_[a];
        }
        
        if(std::abs(eps_m-eps_n)<tol){
        
        std::vector<int> this_pair{ord,num,m,n};
        equal_pairs.push_back(this_pair);
        }
        
      }
      }
    }
    
    
}

if(equal_pairs.size()!=0){
  break;
}

} 
 

  
  
  
}

void AmiBase::print_final( int dim, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array){

print_S(dim, S_array);
print_P( dim, P_array);
print_R( dim, R_array);


}


void AmiBase::print_S( int dim, AmiBase::S_t &S_array){

std::cout<<"Printing S"<<std::endl;
 for (int i=0; i< S_array.size();i++){
std::cout<<"S["<<i<<"]=(";
  for (int j=0; j< S_array[i].size(); j++){
    std::cout<<"{";
    for (int k=0; k< S_array[i][j].size(); k++){

    std::cout<<S_array[i][j][k]<<" ";
    }
    std::cout<<"}";
  }
std::cout<<")";
std::cout<<std::endl;
 }
std::cout<<std::endl;



}

void AmiBase::print_P( int dim, AmiBase::P_t &P_array){

for (int i=0; i< P_array.size(); i++){
std::cout<< "P["<< i<<"]=";

print_Pi(dim, P_array[i]);

}


}


void AmiBase::print_Pi( int dim, AmiBase::Pi_t &Pi_array){

for (int i=0; i< Pi_array.size(); i++){
std::cout<<"[";
print_pole_array(Pi_array[i]);
std::cout<<"]";
std::cout<<std::endl;
}
std::cout<<std::endl;

}

void AmiBase::print_g_prod_array(AmiBase::g_prod_array_t g_array){

for (int i=0; i< g_array.size(); i++){
print_g_prod_info(g_array[i]);
}


}



void AmiBase::print_pole_array(AmiBase::pole_array_t g){

for (int i=0; i< g.size(); i++){
print_pole_struct_info(g[i]);

}


}




void AmiBase::print_g_prod_info(AmiBase::g_prod_t g){

std::cout<<"----------Printing G_product---------- " <<std::endl;
for (std::vector<AmiBase::g_struct>::iterator it= g.begin(); it != g.end(); ++it){

std::cout<<"----------Printing next---------- " <<std::endl;
print_g_struct_info(it[0]);

}


}


void AmiBase::print_R( int dim, AmiBase::R_t &R_array){

for (int i=0; i< R_array.size(); i++){
std::cout<< "R["<< i<<"]=";

print_g_prod_array(R_array[i]);

}

}

std::complex<double> AmiBase::detect_otf_trigger(ami_parms &parms, R_t &R_array,P_t &P_array, S_t &S_array,ami_vars &external){

std::vector< std::vector<int>> triggers;
std::vector<std::vector<int>> pairs;

// std::cout<<"Starting otf loop"<<std::endl;


get_triggers(parms,R_array, P_array, S_array,external, triggers);
// std::cout<<"finished triggers with size: "<<triggers.size()<<" on to equals "<<std::endl;
find_equal_values(parms,R_array,external,pairs);
// std::cout<<"Trigger size is "<<triggers.size()<<std::endl;
// std::cout<<"Pairs size is "<<pairs.size()<<std::endl;

if(triggers.size()==0 && pairs.size()==0){
  // std::cout<<"Evaluating"<<std::endl;
  // std::cout<<"Sizes "<<R_array.size()<<" "<<P_array.size()<<" "<<S_array.size()<<std::endl;
 
print_final(parms.N_INT_, R_array, P_array, S_array); 
  
std::complex<double> calc_result=evaluate(parms,R_array, P_array, S_array,  external);

// std::cout<<"Returning result for OTF evaluation was "<< calc_result<<std::endl;
  return calc_result;
}

int triggers_start;//=triggers[0][0];
int pairs_start;//=pairs[0][0];

int problem_start;

if(triggers.size()!=0 && pairs.size()!=0){
triggers_start=triggers[0][0];
pairs_start=pairs[0][0];
problem_start=std::min(triggers_start, pairs_start);
}

if(triggers.size()==0 && pairs.size()!=0){
pairs_start=pairs[0][0];
problem_start=pairs_start;
}
if(triggers.size()!=0 && pairs.size()==0){
triggers_start=triggers[0][0];
problem_start=triggers_start;
}



// std::cout<<"At end triggers contains "<<std::endl;
// for(int i=0; i< triggers.size(); i++){
  // std::cout<<triggers[i][0]<<" "<<triggers[i][1]<<" "<<triggers[i][2]<<std::endl;  
// }

// std::cout<<"At end pairs contains "<<std::endl;
// for(int i=0; i< pairs.size(); i++){
  // std::cout<<pairs[i][0]<<" "<<pairs[i][1]<<" "<<pairs[i][2]<<" "<<pairs[i][3]<< std::endl;  
// }


R_t R_otf;
P_t P_otf;
S_t S_otf;
//

// std::cout<<"Problem start is "<<problem_start<<std::endl;

// std::cout<<"R_otf size is "<<R_otf.size()<<std::endl;
for(int i=0; i<problem_start+1; i++){
// std::cout<<"Pushing back Ri for i="<<i<<std::endl;
R_otf.push_back(R_array[i]);
if(i<problem_start){
P_otf.push_back(P_array[i]);
S_otf.push_back(S_array[i]);
}

}


if(problem_start==triggers_start){

for(int i=0; i<triggers.size(); i++){


int a=triggers[i][0];
int b=triggers[i][1];
int c=triggers[i][2];  
  // here I just zero the epsilon 
      std::fill(R_otf[a][b][c].eps_.begin(), R_otf[a][b][c].eps_.end(),0.0);
 
}
}

if(problem_start==pairs_start){

for(int i=0; i< pairs.size(); i++){
  
int a=pairs[i][0];
int b=pairs[i][1];
int m=pairs[i][2];
int n=pairs[i][3];

R_otf[a][b][m].eps_=R_otf[a][b][n].eps_;  
  
}
}


// std::cout<<"R_otf size is "<<R_otf.size()<<std::endl;
// Now have the original starting point. 
int dim = parms.N_INT_;

 for (int index = R_otf.size()-1; index < dim; index++) {
   // std::cout<<"On index:"<< index<<std::endl;
   // std::cout<<"Size of R before "<<R_otf.size()<<std::endl;
    update_gprod_general(index, index, R_otf, P_otf, S_otf);
    // std::cout<<"Size of R after "<<R_otf.size()<<std::endl;
    // if the last entry of R has size==0 then there were no suitable poles and the result is just zero 
    if(R_otf.back().size()==0){ return 0.0;}
  }

return detect_otf_trigger(parms, R_otf, P_otf, S_otf, external);


}

std::complex<double> AmiBase::evaluate_otf(ami_parms &parms, R_t &R_array,
                                       P_t &P_array, S_t &S_array,
                                       ami_vars &external) {
  overflow_detected = false;
  
  bool otf=false;
  return detect_otf_trigger(parms, R_array, P_array, S_array, external);
  
}

