#include "ami_base.hpp"



std::complex<double> AmiBase::evaluate(ami_parms &parms, terms &ami_terms,
                                       ami_vars &external, g_prod_t &unique_g,
                                       R_ref_t &Rref, ref_eval_t &Eval_list) {
  overflow_detected = false;

  if (Rref.size() == 0 || unique_g.size() == 0 || Eval_list.size() == 0) {
	  return evaluate(parms, ami_terms, external);
  }

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



  // std::complex<double> output = 0;
  // std::complex<double> term;

  // std::vector<std::complex<double>> unique_vals;

  // std::cout<<"External prefactor is "<<external.prefactor<<std::endl;

  for (int i = 0; i < unique_g.size(); i++) {
    // we pretend each G is a g_prod and use the same function as the normal
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

    unique_vals.push_back(gprod); // This removes the overall prefactor for each g. we then add that back later
  }

  // now evaluate each term in list

  for (int i = 0; i < Eval_list.size(); i++) {
#ifdef BOOST_MP
    std::complex<boost::multiprecision::float128> ksum(0, 0);
#else
    std::complex<double> ksum(0, 0);
#endif
    for (int j = 0; j < Eval_list[i].size(); j++) {
      if (verbose) {
        std::cout << j << " "
                  << eval_fprod(parms, ami_terms[Eval_list[i][j].first].p_list,
                                external) *
                         ami_terms[Eval_list[i][j].first].sign
                  << " " << double(Eval_list[i][j].second) << std::endl;
      }
      std::complex<double> kterm =
          eval_fprod(parms, ami_terms[Eval_list[i][j].first].p_list, external);

	if (verbose) {
		std::cout<<"Kterm: "<<j<<": "<< kterm<<std::endl;
	}


      ksum += kterm * ami_terms[Eval_list[i][j].first].sign *
              double(Eval_list[i][j].second);
    }

    ref_v_t pair_vec = Rref[i];
    
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
 * This is the primary AMI evaluation function for the term-by-term
 * construction.
 * @param[in] parms : `ami_parms` object, basic parameters for AMI. Typically
 * same as construction.
 * @param[in] ami_terms : `terms` result from construction function.
 * @param[in] external : Input external variables in a `ami_vars` struct.
 * @return Result is returned as single value of type `std::complex<double>`.
 */
std::complex<double> AmiBase::evaluate(ami_parms &parms, terms &ami_terms,
                                       ami_vars &external) {
  overflow_detected = false;
  
#ifdef BOOST_MP
  std::complex<boost::multiprecision::float128> term;
 
  std::complex<boost::multiprecision::float128> output(0, 0);
 #else
  std::complex<double> term;
  std::complex<double> output(0,0);
  
#endif
  

  for (int i = 0; i < ami_terms.size(); i++) {
   
    #ifdef BOOST_MP
	
	term =
        evaluate_term_mp(parms, ami_terms[i], external);
	
	#else
    term =
        evaluate_term(parms, ami_terms[i], external);
	#endif
    output += term;

    if (verbose) {
		print_term(ami_terms[i]);
      std::cout <<std::setprecision(32)<< "Term gave " << i << " " << term
                << " and current total is : " << output << std::endl;
    }
  }

#ifdef BOOST_MP
  std::complex<double> final_output(output.real().convert_to<double>(),
                                    output.imag().convert_to<double>());
#else
  std::complex<double> final_output = output;
#endif

 
  return final_output;
}

std::complex<double> AmiBase::evaluate_term(ami_parms &parms, term &ami_term,
                                            ami_vars &external) {
  std::complex<double> gprod;

  gprod = eval_gprod(parms, ami_term.g_list, external);

  // check for overflow
  if ((std::floor(std::abs(std::real(gprod))) == std::abs(std::real(gprod))) &&
      std::abs(std::real(gprod)) != 0) {
    overflow_detected = true;
  }

  std::complex<double> fprod;

  fprod = eval_fprod(parms, ami_term.p_list, external);

  std::complex<double> output(0, 0);

  output = ami_term.sign * gprod * fprod;

  if (verbose) {
    std::cout << "Term gave " << ami_term.sign * fprod << " " << gprod << " "
              << output << std::endl;
  }
  return output;
}

#ifdef BOOST_MP
std::complex<boost::multiprecision::float128> AmiBase::evaluate_term_mp(ami_parms &parms, term &ami_term,
                                            ami_vars &external) {
  std::complex<boost::multiprecision::float128> gprod, fprod;
  std::complex<boost::multiprecision::float128> output(0,0);

  gprod = eval_gprod_mp(parms, ami_term.g_list, external);

  // check for overflow
  if ((boost::multiprecision::floor(boost::multiprecision::abs(std::real(gprod))) == boost::multiprecision::abs(std::real(gprod))) &&
      boost::multiprecision::abs(std::real(gprod)) != 0) {
    overflow_detected = true;
  }

  // std::complex<double> fprod;

  fprod = eval_fprod(parms, ami_term.p_list, external);

  // std::complex<double> output(0, 0);
  boost::multiprecision::float128 dsign=ami_term.sign;
  std::complex<boost::multiprecision::float128> sign=dsign;
  output = sign * gprod * fprod;

  if (verbose) {
    std::cout << "MP Term gave " << sign * fprod << " " << gprod << " "
              << output << std::endl;
  }
  return output;
}
#endif



/**
 *
 * @param[in] parms : `ami_parms` object, basic parameters for AMI.
 * @param[in] p_list : `pole_array_t` a list of poles that is interpretted as
 * \f$ \prod{f(p_i)} \f$.  If poles contain derivatives then these are resolved
 * at this stage.
 * @param[in] external : Input external variables in a `ami_vars` struct.
 * @return Single value for product of fermi/bose functions.
 */
std::complex<double> AmiBase::eval_fprod(ami_parms &parms, pole_array_t &p_list,
                                         ami_vars &external) {
  std::complex<double> output(1, 0);

  for (int i = 0; i < p_list.size(); i++) {
    std::complex<double> pole_val = fermi_pole(parms, p_list[i], external);
    if (verbose) {
      std::cout << "On pole " << i << " returned value " << pole_val
                << std::endl;
    }

    output = output * pole_val; // fermi_pole(parms, p_list[i], external);
  }

  return output;
}
