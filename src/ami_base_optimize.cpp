#include "ami_base.hpp"


void AmiBase::reduce_rref(R_ref_t &Rref, ref_eval_t &Eval_list){
	
std::vector<int> used;	
	
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
			// std::cout<<"compare gave"<<ditto<<" "<<sign1<<" "<<sign2<<std::endl;
			if(ditto){
				used.push_back(j);
				ref_t this_ref=std::make_pair(j,sign2);
				this_vec.push_back(this_ref);
			}
			
			
			
		}
		
		
		
		
		
	}
	// don't want empty though don't think it is possible 
	if(this_vec.size()>0){
	Eval_list.push_back(this_vec);
	}
}
	
	
	
}

bool AmiBase::pair_v_equiv(ref_v_t &r1, ref_v_t &r2, int &r1sign, int &r2sign){
	r1sign=0;
	r2sign=0;
	
	if(r1.size()!=r2.size()){return false;}
	int this_sign=1;
	int this_sign2=1;
	for(int i=0; i< r1.size(); i++){
		
	if(r1[i].first!=r2[i].first){return false;}
		this_sign=this_sign*r1[i].second;
		this_sign2=this_sign2*r2[i].second;
		
		
	}
	
	r1sign=this_sign;
	r2sign=this_sign2;
	return true;
	
	
}

void AmiBase::factorize_Rn(Ri_t &Rn, g_prod_t &unique_g, R_ref_t &Rref,ref_eval_t &Eval_list){
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
		break;// break out of the unique_g loop since we found an equivalent Rn[i][j] g_struct 
		}
		
		
	}
	// if made it here it means this is a unique G and so needs to be stored and a ref created 
	if(!found){
		unique_g.push_back(Rn[i][j]);
		ref=std::make_pair(unique_g.size()-1,1);
		ref_v.push_back(ref);		
	}
	
	
	
	
	
	
	}	
	
	Rref.push_back(ref_v);
	
}	


reduce_rref(Rref,Eval_list);
	
}

/* 
g_struct(){} 

epsilon_t eps_;
// std::vector<int> eps_indices_;
alpha_t alpha_;
stat_type stat_;
species_t species_;

};
 */

bool AmiBase::g_struct_equiv(g_struct &g1, g_struct &g2, int &sign){
	
// std::cout<<"Comparing graphs "<<std::endl;
// print_g_struct_info(g1);
// std::cout<<"g1 stat "<<g1.stat_<<std::endl;
// std::cout<<"Vs"<<std::endl;
// print_g_struct_info(g2);
// std::cout<<"g2 stat "<<g2.stat_<<std::endl;
// std::cout<<"------"<<std::endl;
	
sign=0;// by default zero sign means they are not equiv

bool result=true;

if(g1.eps_.size()!=g2.eps_.size()){ return false;}
if(g1.alpha_.size()!=g2.alpha_.size()){ return false;}
// if(g1.stat_ != g2.stat_){ return false;} // the stat type of the g has no meaning at this stage since the frequencies determine the statistics 
if(g1.species_ != g2.species_){ return false;}
	
std::vector<int> signs;	
	
for( int i=0; i< g1.eps_.size(); i++){
	
	if( std::abs(g1.eps_[i]) != std::abs(g2.eps_[i])){ return false; break;}else{
		if(g1.eps_[i]!=0){
		signs.push_back(sgn(g1.eps_[i]/g2.eps_[i]));}
	}
	
}


for( int i=0; i< g1.alpha_.size(); i++){
	
	if( std::abs(g1.alpha_[i]) != std::abs(g2.alpha_[i])){ return false; break;}else{
		if(g1.alpha_[i]!=0){
		signs.push_back(sgn(g1.alpha_[i]/g2.alpha_[i]));}
	}
	
}

// in principle if the code gets here the size of signs is >0 so don't need to check that

if ( std::adjacent_find( signs.begin(), signs.end(), std::not_equal_to<int>() ) == signs.end() )
{
    sign=signs.back();
}

// std::cout<<"Found equiv"<<std::endl;
// std::cout<<"-----"<<std::endl;
	
return result;	
}


void AmiBase::print_g_struct_info(g_struct g){

std::cout<<"Species="<<g.species_<<" ";
std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}


void AmiBase::print_epsilon_info(AmiBase::epsilon_t eps){


//for (std::vector<signed char>::iterator it= eps.begin(); it != eps.end(); ++it){
for (std::vector<int>::iterator it= eps.begin(); it != eps.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}


void AmiBase::print_alpha_info(AmiBase::alpha_t alpha){

for (std::vector<int>::iterator it= alpha.begin(); it != alpha.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}



// bool AmiBase::pole_equiv (pole_struct pole1, pole_struct pole2){

// bool result=true;

// if(pole1.eps_.size()!= pole2.eps_.size()){ return false;}
// if(pole1.alpha_.size()!= pole2.alpha_.size()){ return false;}

// for (int i=0; i< pole1.eps_.size(); i++)
// {
// if (pole1.eps_[i] != pole2.eps_[i]){ return false; break;}

// }

// for (int i=0; i< pole1.alpha_.size(); i++)
// {
// if (pole1.alpha_[i] != pole2.alpha_[i]){ return false; break;}

// }



// return result;
// }