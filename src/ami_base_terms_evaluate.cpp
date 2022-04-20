#include "ami_base.hpp"

std::complex<double> AmiBase::evaluate(ami_parms &parms, terms &ami_terms,
                                       ami_vars &external, g_prod_t &unique_g,
                                       R_ref_t &Rref, ref_eval_t &Eval_list) {
  overflow_detected = false;

  if (Rref.size() == 0 || unique_g.size() == 0 || Eval_list.size() == 0) {
    return evaluate(parms, ami_terms, external);
  }

  std::complex<double> output = 0;
  std::complex<double> term;

  std::vector<std::complex<double>> unique_vals;

  // std::cout<<"External prefactor is "<<external.prefactor<<std::endl;

  for (int i = 0; i < unique_g.size(); i++) {
    // we pretend each G is a g_prod and use the same function as the normal
    // star operation for evaluating products of G's
    g_prod_t this_term;
    this_term.push_back(unique_g[i]);

    unique_vals.push_back(
        eval_gprod(parms, this_term, external) *
        external.prefactor); // This removes the overall prefactor for each
                             // g. we then add that back later
  }

  // now evaluate each term in list

  for (int i = 0; i < Eval_list.size(); i++) {
    std::complex<double> ksum(0, 0);
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

      ksum += kterm * ami_terms[Eval_list[i][j].first].sign *
              double(Eval_list[i][j].second);
    }

    ref_v_t pair_vec = Rref[i];
    std::complex<double> this_gprod(1, 0);
    for (int j = 0; j < pair_vec.size(); j++) {
      this_gprod = this_gprod * unique_vals[pair_vec[j].first];
    }

    if ((std::floor(std::abs(std::real(this_gprod))) ==
         std::abs(std::real(this_gprod))) &&
        std::abs(std::real(this_gprod)) != 0) {
      overflow_detected = true;
    }

    term = ksum * this_gprod * external.prefactor;

    if (verbose) {
      std::cout << "In term evaluate K[]*R" << std::endl;
      std::cout << std::setprecision(20) << i << " " << ksum << " "
                << std::real(this_gprod) << " " << std::imag(this_gprod) << " "
                << std::real(term) << " " << std::imag(term) << std::endl;
    }

    output += term;
  }

  return output;
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

  std::complex<double> output(0, 0);

  for (int i = 0; i < ami_terms.size(); i++) {
   
    std::complex<double> this_term =
        evaluate_term(parms, ami_terms[i], external);
    output += this_term;

    if (verbose) {
      std::cout << "Term gave " << i << " " << this_term
                << " and current total is : " << output << std::endl;
    }
  }

  return output;
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
