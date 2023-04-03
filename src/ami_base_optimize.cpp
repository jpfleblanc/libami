#include "ami_base.hpp"


void AmiBase::derivative_opt(g_prod_t &unique_g, R_ref_t &Rref,
                             ref_eval_t &Eval_list) {
  R_ref_t nRref;
  ref_eval_t nEval_list;


  for (int n = 0; n < Rref.size(); n++) {

    // going to take derivative of this term
    // for each first - check unique_g to see if the external freq is 0 or
    // +-1. If it is not zero then

    for (int pair_ind = 0; pair_ind < Rref[n].size(); pair_ind++) {
 
      ref_v_t this_prod;
      ref_v_t this_eval_list = Eval_list[n];
      double sign_change = 1.0;

      g_struct this_g = unique_g[Rref[n][pair_ind].first];


      if (this_g.alpha_.back() != 0) {
        for (int k = 0; k < Rref[n].size(); k++) {
          
          if (k == pair_ind) {
            
            this_prod.push_back(Rref[n][pair_ind]);
            ref_t new_pair = Rref[n][pair_ind];
            new_pair.second = new_pair.second * (-1) * this_g.alpha_.back();
            sign_change = (-1) * this_g.alpha_.back();
            this_prod.push_back(new_pair);
          } else {
            
            this_prod.push_back(Rref[n][k]);
          }

          // now put all the other terms in
        }
      }

      if (this_prod.size() != 0) {
        for (int m = 0; m < this_eval_list.size(); m++) {
          this_eval_list[m].second = this_eval_list[m].second * sign_change;
        }

        nRref.push_back(this_prod);
        nEval_list.push_back(this_eval_list);
      }
    }
  }

  Rref = nRref;
  Eval_list = nEval_list;
}

// this version kicks out the extra Rref entries
void AmiBase::reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list) {
  std::vector<int> used;
  

  for (int i = 0; i < Rref.size(); i++) {
    
    ref_v_t this_vec;
    if (std::find(used.begin(), used.end(), i) != used.end()) {
      continue;
    } else {
      used.push_back(i);
      int this_sign = 1;
      for (int pair = 0; pair < Rref[i].size(); pair++) {
        
        this_sign = this_sign * Rref[i][pair].second;
      }
      ref_t this_ref = std::make_pair(i, this_sign);
      this_vec.push_back(this_ref);
    }

    for (int j = i + 1; j < Rref.size(); j++) {
      
      if (std::find(used.begin(), used.end(), j) != used.end()) {
        continue;
      } else {
        int sign1, sign2;
        bool ditto = pair_v_equiv(Rref[i], Rref[j], sign1, sign2);
        
        if (ditto) {
          used.push_back(j);
          ref_t this_ref = std::make_pair(j, sign2);
          
          this_vec.push_back(this_ref);
          
        }
      }
    }
    // don't want empty though don't think it is possible
    if (this_vec.size() > 0) {
      Eval_list.push_back(this_vec);
    }
  }

  // now make a minimal Rref_v

  R_ref_t newref;

  for (int i = 0; i < Eval_list.size(); i++) {
    ref_v_t newentry;

    for (int j = 0; j < Rref[Eval_list[i][0].first].size(); j++) {
      ref_t this_pair;
      this_pair.first = Rref[Eval_list[i][0].first][j].first;
      this_pair.second = 1;
      newentry.push_back(this_pair);
    }

    newref.push_back(newentry);
  }

  Rref = newref;
}
// Deprecated function 
/*
void AmiBase::reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list){

std::vector<int> used;
std::cout<<"Rref size is "<<Rref.size()<<std::endl;

for (int i=0; i< Rref.size(); i++){
        // std::cout<<"On i="<<i<<std::endl;
        ref_v_t this_vec;
        if (std::find(used.begin(), used.end(),i)!=used.end()){
                continue;
        }else{
                used.push_back(i);
                int this_sign=1;
                for (int pair=0; pair< Rref[i].size(); pair++){
                        this_sign=this_sign*Rref[i][pair].second;
                }
                ref_t this_ref=std::make_pair(i,this_sign);
                this_vec.push_back(this_ref);
        }

        for(int j=i+1; j<Rref.size(); j++){
                // std::cout<<"On j="<<j<<std::endl;
                if(std::find(used.begin(), used.end(), j)!=used.end()){
                        continue;
                }
                else{
                        int sign1,sign2;
                        bool ditto=pair_v_equiv(Rref[i],Rref[j],sign1,sign2);
                        // std::cout<<"compare gave"<<ditto<<" "<<sign1<<"
"<<sign2<<std::endl; if(ditto){ used.push_back(j); ref_t
this_ref=std::make_pair(j,sign2);
                                // if(this_ref.size()>0){
                                this_vec.push_back(this_ref);
                                // }
                        }



                }





        }
        // don't want empty though don't think it is possible
        if(this_vec.size()>0){
        Eval_list.push_back(this_vec);
        }
}



} */



// This function needs update.  It does not catch when the order of the R's is different.
// APR 2023 Update

bool AmiBase::pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign) {
  r1sign = 0;
  r2sign = 0;

  if (r1.size() != r2.size()) {
    return false;
  }
  
  std::vector<int> v2;
  // v1.resize(r1.size(),1);
  v2.resize(r2.size(),1);
  

  int this_sign = 1;
  int this_sign2 = 1;
  for(int i=0; i<r1.size(); i++){
    bool found=false;
    for(int j=0; j< r2.size(); j++){
    
    if ( (r1[i].first == r2[j].first) && v2[j]==1){
      found=true;
      v2[j]=0;
      continue; //
      
    }
    
    
    }
    
    if(!found){return false;}
    
    this_sign = this_sign * r1[i].second;
    this_sign2 = this_sign2 * r2[i].second;
  }
  
  // if it gets here then the two list are equivalent within the relative signs 
  r1sign = this_sign;
  r2sign = this_sign2;
  return true;
  
  

  // int this_sign = 1;
  // int this_sign2 = 1;
  // for (int i = 0; i < r1.size(); i++) {
    // if (r1[i].first != r2[i].first) {
      // return false;
    // }
    // this_sign = this_sign * r1[i].second;
    // this_sign2 = this_sign2 * r2[i].second;
  // }

  // r1sign = this_sign;
  // r2sign = this_sign2;
  // return true;
}


// bool AmiBase::pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign) {
  // r1sign = 0;
  // r2sign = 0;

  // if (r1.size() != r2.size()) {
    // return false;
  // }

  // int this_sign = 1;
  // int this_sign2 = 1;
  // for (int i = 0; i < r1.size(); i++) {
    // if (r1[i].first != r2[i].first) {
      // return false;
    // }
    // this_sign = this_sign * r1[i].second;
    // this_sign2 = this_sign2 * r2[i].second;
  // }

  // r1sign = this_sign;
  // r2sign = this_sign2;
  // return true;
// }

void AmiBase::factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,
                           ref_eval_t &Eval_list) {
  unique_g.clear();
  Rref.clear();
  Eval_list.clear();

  // i is the entry of Rn
  // first step is to generate
  for (int i = 0; i < Rn.size(); i++) {
    ref_v_t ref_v;
    // Rn[i] is a g_prod_t, so Rn[i][j] is a G -> g_struct
    for (int j = 0; j < Rn[i].size(); j++) {
      ref_t ref;

      // if unique list is empty then add entry to it.
      if (unique_g.size() == 0) {
        unique_g.push_back(Rn[i][j]);
        ref = std::make_pair(0, 1);
        ref_v.push_back(ref);
        continue;
      }
      // for each G search through unique_g
      bool found = false;
      for (int m = 0; m < unique_g.size(); m++) {
        int this_sign = 0;
        bool eq = g_struct_equiv(Rn[i][j], unique_g[m], this_sign);

        if (eq) {
          ref = std::make_pair(m, this_sign);
          ref_v.push_back(ref);
          found = true;
          break; // break out of the unique_g loop since we found an
                 // equivalent Rn[i][j] g_struct
        }
      }

      // if made it here it means this is a unique G and so needs to be
      // stored and a ref created
      if (!found) {
        unique_g.push_back(Rn[i][j]);
        ref = std::make_pair(unique_g.size() - 1, 1);
        ref_v.push_back(ref);
      }
    }

    if (ref_v.size() != 0) {
      Rref.push_back(ref_v);
    }
  }


  reduce_rref(Rref, Eval_list);

}

// Deprecated
/*
void AmiBase::factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t
&Rref,ref_eval_t &Eval_list){
//i is the entry of Rn
for(int i=0; i< Rn.size(); i++){

        ref_v_t ref_v;
        // Rn[i] is a g_prod_t, so Rn[i][j] is a G -> g_struct
        for(int j=0; j<Rn[i].size(); j++){
        ref_t ref;

        // if unique list is empty then add entry to it.
        if(unique_g.size()==0){
                unique_g.push_back(Rn[i][j]);
                ref=std::make_pair(0,1);
                ref_v.push_back(ref);
                continue;
        }
        // for each G search through unique_g
        bool found=false;
        for(int m=0; m< unique_g.size();m++){

                int this_sign=0;
                bool eq=g_struct_equiv(Rn[i][j], unique_g[m], this_sign);

                if(eq){

                ref=std::make_pair(m,this_sign);
                ref_v.push_back(ref);
                found=true;
                break;// break out of the unique_g loop since we found an
equivalent Rn[i][j] g_struct
                }


        }



        // if made it here it means this is a unique G and so needs to be
stored and a ref created if(!found){ unique_g.push_back(Rn[i][j]);
                ref=std::make_pair(unique_g.size()-1,1);
                ref_v.push_back(ref);
        }






        }


        Rref.push_back(ref_v);

}


reduce_rref(Rref,Eval_list);

}
 */


bool AmiBase::g_struct_equiv(g_struct &g1, g_struct &g2, int &sign) {


  sign = 0; // by default zero sign means they are not equiv

  bool result = true;

  if (g1.eps_.size() != g2.eps_.size()) {
    return false;
  }
  if (g1.alpha_.size() != g2.alpha_.size()) {
    return false;
  }

  if (g1.species_ != g2.species_) {
    return false;
  }

  std::vector<int> signs;

  for (int i = 0; i < g1.eps_.size(); i++) {
    if (std::abs(g1.eps_[i]) != std::abs(g2.eps_[i])) {
      return false;
      break;
    } else {
      if (g1.eps_[i] != 0) {
        signs.push_back(sgn(g1.eps_[i] / g2.eps_[i]));
      }
    }
  }

  for (int i = 0; i < g1.alpha_.size(); i++) {
    if (std::abs(g1.alpha_[i]) != std::abs(g2.alpha_[i])) {
      return false;
      break;
    } else {
      if (g1.alpha_[i] != 0) {
        signs.push_back(sgn(g1.alpha_[i] / g2.alpha_[i]));
      }
    }
  }

  // in principle if the code gets here the size of signs is >0 so don't need
  // to check that


  if (std::adjacent_find(signs.begin(), signs.end(),
                         std::not_equal_to<int>()) == signs.end()) {
    sign = signs.back();
  } else {
    return false;
  }


  return result;
}

void AmiBase::print_g_struct_info(g_struct g) {
  std::cout << "Species=" << g.species_ << " ";
  std::cout << "Eps=(";
  print_epsilon_info(g.eps_);
  std::cout << ")";
  std::cout << std::endl;
  std::cout << "Alpha=(";
  print_alpha_info(g.alpha_);
  std::cout << ")";
  std::cout << std::endl;
}

void AmiBase::print_pole_struct_info(pole_struct g) {
  std::cout << "Eps=(";
  print_epsilon_info(g.eps_);
  std::cout << ")";
  std::cout << std::endl;
  std::cout << "Alpha=(";
  print_alpha_info(g.alpha_);
  std::cout << ")";
  std::cout << "X_Alpha=(";
  print_alpha_info(g.x_alpha_);
  std::cout << ")";
  std::cout << std::endl;
  std::cout << "Der: " << g.der_ << std::endl;
}

void AmiBase::print_epsilon_info(AmiBase::epsilon_t eps) {

  for (std::vector<int>::iterator it = eps.begin(); it != eps.end(); ++it) {
    std::cout << *it << ' ';
  }
 
}

void AmiBase::print_alpha_info(AmiBase::alpha_t alpha) {
  for (std::vector<int>::iterator it = alpha.begin(); it != alpha.end(); ++it) {
    std::cout << *it << ' ';
  }
  
}
