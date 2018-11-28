#include "ami.hpp"
#include <iomanip>



void AmiCalc::write_S(S_t &S_array){

for (int i=0; i< S_array.size(); i++){

std::ofstream file;
std::stringstream filename;
filename<<"S_"<< i <<".dat";

file.open(filename.str());

for (int j=0; j< S_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<S_array[i][j].size();k++){

file<<S_array[i][j][k]<<" ";
}
file <<"]"<<std::endl;

}



file.close();

}

}


void AmiCalc::write_P(P_t &P_array){

for (int i=0; i< P_array.size(); i++){

std::ofstream file;
std::stringstream filename;
filename<<"P_"<<i<<".dat";

file.open(filename.str());

for (int j=0; j< P_array[i].size(); j++){

file<<"[ ";
for (int k=0; k<P_array[i][j].size();k++){
std::stringstream  pole_string;
file<<" ( ";
for (int r=0; r<P_array[i][j][k].alpha_.size(); r++){
pole_string<<P_array[i][j][k].alpha_[r]<<" ";

}

file<< pole_string.str();

file<<" ) ";
}
file <<"]"<<std::endl;

}



file.close();

}

}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



void AmiCalc::print_S( int dim, S_t &S_array){

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

void AmiCalc::print_P( int dim, P_t &P_array){

for (int i=0; i< P_array.size(); i++){
std::cout<< "P["<< i<<"]=";

print_Pi(dim, P_array[i]);

}


}

void AmiCalc::print_Pi( int dim, Pi_t &Pi_array){

for (int i=0; i< Pi_array.size(); i++){
std::cout<<"[";
print_pole_array(Pi_array[i]);
std::cout<<"]";
std::cout<<std::endl;
}
std::cout<<std::endl;

}



void AmiCalc::print_R( int dim, R_t &R_array){

for (int i=0; i< R_array.size(); i++){
std::cout<< "R["<< i<<"]=";

print_g_prod_array(R_array[i]);

}

}

void AmiCalc::print_g_prod_array(g_prod_array_t g_array){

for (int i=0; i< g_array.size(); i++){
print_g_prod_info(g_array[i]);
}


}

void AmiCalc::print_final( int dim, R_t &R_array, P_t &P_array, S_t &S_array){

print_S(dim, S_array);
print_P( dim, P_array);
print_R( dim, R_array);


}




void AmiCalc::print_g_prod_info(AmiCalc::g_prod_t g){

std::cout<<"----------Printing G_product---------- " <<std::endl;
for (std::vector<AmiCalc::g_struct>::iterator it= g.begin(); it != g.end(); ++it){

std::cout<<"----------Printing next---------- " <<std::endl;
print_g_struct_info(it[0]);

}


}


void AmiCalc::print_pole_array(AmiCalc::pole_array_t g){

for (int i=0; i< g.size(); i++){
print_pole_struct_info(g[i]);

}


}


void AmiCalc::print_g_struct_info(AmiCalc::g_struct g){

std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}

void AmiCalc::print_pole_struct_info(AmiCalc::pole_struct g){

std::cout<<"Eps=(";
print_epsilon_info(g.eps_);
std::cout<<")";
std::cout<<std::endl;
std::cout<<"Alpha=(";
print_alpha_info(g.alpha_);
std::cout<<")";
std::cout<<std::endl;

}

void AmiCalc::print_sign_array(AmiCalc::sign_array_t signs){

std::cout<<"Printing signs"<<std::endl;
std::cout<<"[";
 for (int i=0; i< signs.size();i++){
std::cout<<"(";
  for (int j=0; j< signs[i].size(); j++){
    std::cout<<signs[i][j]<<" ";

  }
std::cout<<")";
 }

std::cout<<"]";


}



void AmiCalc::print_signs(AmiCalc::sign_t signs){

std::cout<<"(";
 for (int i=0; i< signs.size();i++){

     std::cout<<signs[i]<<" ";

  }
std::cout<<")";


}



void AmiCalc::print_epsilon_info(AmiCalc::epsilon_t eps){


for (std::vector<int>::iterator it= eps.begin(); it != eps.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}


void AmiCalc::print_alpha_info(AmiCalc::alpha_t alpha){

for (std::vector<int>::iterator it= alpha.begin(); it != alpha.end(); ++it){

std::cout<< *it << ' ';

}
//std::cout << '\n';


}
