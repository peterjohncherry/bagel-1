#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<bitset>
#include<array>
#include<vector>
#include<algorithm>
#include<utility>
#include<tuple>
#include<string>
#include<memory>
#include<numeric>
#include<map>
#include<list>
#include<sstream>
#include <iostream>
#include <fstream>
#include "latex_eqn.h"

using namespace std; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  latex_eqn::init_outfile( string name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

   ofstream outfile;
   outfile.open (name);
   outfile << "\\documentclass[12pt]{article}" << endl;
   outfile << "\\usepackage[utf8x]{inputenc}"<<endl;
   outfile << "\\usepackage[english]{babel}"<<endl;
   outfile << "\\usepackage[T1]{fontenc}"<<endl;
   outfile << "\\usepackage{color}"<<endl;
   outfile << "\\usepackage{amsmath}"<<endl;
   outfile << "\\usepackage{amssymb}"<<endl;
   outfile << "\\usepackage[font=small,format=plain,labelfont=bf,up,textfont=it,up]{caption}"<<endl;
   outfile << "\\usepackage{calc}"<<endl;
   outfile << "\\usepackage{relsize}"<<endl;
   outfile << "\\usepackage{bbold}"<<endl;
   outfile << "\\usepackage{mathtext}"<<endl;
   outfile << "\\linespread{1.3}"<<endl;
   outfile << "\\begin{document}"<<endl;
 return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void latex_eqn::latex_eqn_out(string filename) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int kk = 0; 
  while ( kk != allops->size()){

   ofstream eqfile;
   eqfile.open (filename, std::fstream::in | std::fstream::out | std::fstream::app);
   static int linelen;
   int eqlen =latex_rdm_string(allids->at(kk), allops->at(kk) ).length() ;

   if( kk == 0){
     eqfile << "\\begin{equation*}"<<endl;
     linelen=0;
   }

   linelen += eqlen; 

   if ( (linelen > 150) && (eqlen < 150) ) {
     eqfile << endl<< "\\end{equation*}"<<endl;
     eqfile << "\\begin{equation*}"<<endl ;
     linelen = latex_rdm_string(allids->at(kk), allops->at(kk) ).length() ;
   }  
   

   if ( allsigns->at(kk)== 1 ){
      if(kk !=0) eqfile << " + " ;
   } else if (allsigns->at(kk)== -1){
      eqfile << " - " ;
   } else {
     eqfile << to_string(allsigns->at(kk)) << " " ; 
   }

   auto supers = make_shared<vector<string>>();
   auto opnames = make_shared<vector<char>>();
   opnames->push_back('f');
   opnames->push_back('T');
   opnames->push_back('L');

   eqfile << latex_rdm_string(allids->at(kk), allops->at(kk) ) ;
   eqfile << latex_delta_string(alldeltas->at(kk));
   eqfile << latex_genop_string(orig_ids, supers, opnames ) ;

   if( (kk== allids->size()-1) ) {
      eqfile << endl << "\\end{equation*}" << endl;
      eqfile << "\\end{document}";
   }
   eqfile.close();

   kk++;
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string latex_eqn::latex_delta_string(pair<string,string> deltas) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto delvec = make_shared<vector<pair<string,string>>>();
  delvec->push_back(deltas);
  return latex_delta_string(delvec) ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string latex_eqn::latex_delta_string(shared_ptr<vector<pair<string,string>>> deltas) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto spin = [](string idx) {
    string sp;
      if ((idx.back() == 'A') || ( idx.back() == 'a')){
       sp =  "_{\\alpha}";
      } else if ( (idx.back() == 'B') || ( idx.back() == 'b')) {
       sp = "_{\\beta}";
      }
      return sp;
    };

  string dels= "";
  for ( auto ii = 0 ; ii != deltas->size(); ii++ ){
    string idx1 = deltas->at(ii).first;
    string idx2 = deltas->at(ii).second;
    dels+= "\\delta_{";
    if(spinfree){
      dels+= idx1.substr(1,idx1.length()-1);
      dels+= idx2.substr(1,idx1.length()-1);
    }else if(!spinfree){
      dels+= idx1.substr(1, idx1.length()-2) ;
      dels+=spin(idx1);
      dels+= idx2.substr(1,idx2.length()-2); 
      dels += spin(idx2);
    }
    dels+= "} ";
  }
  cout <<"dels = " << dels << endl; 
  return dels;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string latex_eqn::latex_genop_string(shared_ptr<vector<string>> subs, shared_ptr<vector<string>> supers, shared_ptr<vector<char>> opnames ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  string ops= "";
  cout << "hello"<< endl;
  int bob;
  for (auto name : *opnames){
    cout <<" name =  " << name << endl;
    cout << "subs->size() = "<<  subs->size() << endl;
    int newop=0;
    for (auto  ii = 0 ; ii != subs->size() ; ii++){
      cout << "subs->at("<<ii<<") = "<<  subs->at(ii) << endl;
      if(subs->at(ii)[0] == name) { 
        if(newop==0){
          ops+= name;
          ops+= "_{";
          newop = 1;
        }
        auto sub = subs->at(ii);
        if (spinfree) ops+= sub.substr(1,sub.length()-1) ;
        
        if (!spinfree){
          ops+= sub.substr(1,sub.length()-2) ;
        
          if ((sub.back() == 'A') || ( sub.back() == 'a')){
           ops = ops + "_{\\alpha}";
          } else if ( (sub.back() == 'B') || ( sub.back() == 'b')) {
           ops = ops + "_{\\beta}";
          }
        } 
      }
    }
    if(newop==1)
      ops += "}";

    for (auto  ii = 0 ; ii != supers->size() ; ii++){
      if(supers->at(ii)[0] == name) { 
      cout << "supers->at("<<ii<<") = "<< subs->at(ii) << endl;
        if(newop == 0 )
          ops+= name;
          ops+= "^{";
          newop =2;
        }else if(newop == 1){
          ops+= "^{";
          newop = 2;
        }
        ops += supers->at(ii).substr(1,supers->at(ii).length()-1);
      } 
      if (newop ==2) 
        ops += "}";
      cout << "ops string = " <<  ops << endl; 
    }
    return ops;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string latex_eqn::latex_aop_string(string idx, bool dag ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  string op = "a_{";

  if (spinfree) op+= idx.substr(1,idx.length()) ;

  if (!spinfree){
    op+= idx.substr(1,idx.length()-2);

    if ((idx.back() == 'A') || ( idx.back() == 'a')){
     op = op + "_{\\alpha}";
    } else if ( (idx.back() == 'B') || ( idx.back() == 'b')) {
     op = op + "_{\\beta}";
    }
  } 

  op += "}";
  if(dag)
    op += "^{\\dagger}";
   
  return op;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string latex_eqn::latex_rdm_string(shared_ptr<vector<string>> ids, shared_ptr<vector<bool>> acs ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string rdm =  "\\langle I | "; 
  for (auto ii =0 ; ii !=acs->size(); ii++  )
    rdm +=  latex_aop_string( ids->at(ii), acs->at(ii)); 
  rdm += " | K  \\rangle ";
  return rdm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void latex_eqn::latex_print_all(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int kk = 0;
  while ( kk != allops->size()){
    auto ac  = allops->at(kk);
    auto ids = allids->at(kk);
    auto sign  = allsigns->at(kk);
    auto delta  = alldeltas->at(kk);
  }   
  return;  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////             MAIN                  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(){
  cout << "bob can dance" << endl;
  return 0;
}                                                                                                                                                                                                        
