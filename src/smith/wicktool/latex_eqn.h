#ifndef __SRC_SMITH_latexeqn_H
#define __SRC_SMITH_latexeqn_H


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
#include<sstream>
#include<map>

// class for contraction terms
class latex_eqn{ 
  private:
  // variables
  std::shared_ptr<std::vector<std::string>> orig_ids ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>>   allops ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> allids ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string,std::string>>>>> alldeltas; 
  std::shared_ptr<std::vector<int>> allsigns  ;
  bool spinfree;

  public: 
    latex_eqn( std::shared_ptr<std::vector<std::string>> orig_ids_in ,
                   std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>>   allops_in ,
                   std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> allids_in, 
                   std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string,std::string>>>>> alldeltas_in,
                   std::shared_ptr<std::vector<int>> allsigns_in,
                   bool spinfree_in) : 
                   orig_ids(orig_ids_in), allops(allops_in), allids(allids_in), alldeltas(alldeltas_in), allsigns(allsigns_in), 
                   spinfree(spinfree_in) {};

    ~latex_eqn(){};

// functions
 void latex_eqn_out(std::string filename) ;

 void latex_print_all();

 std::string latex_aop_string(std::string idx, bool dag ) ;

 std::string latex_rdm_string(std::shared_ptr<std::vector<std::string>> ids, std::shared_ptr<std::vector<bool>> acs ) ;

 void init_outfile(std::string  outfile_name ) ;


 std::string latex_delta_string(std::shared_ptr<std::vector<std::pair<std::string,std::string>>> deltas) ;

 std::string latex_delta_string(std::pair<std::string,std::string> deltas) ;

 std::string latex_genop_string(std::shared_ptr<std::vector<std::string>> subs,
                                std::shared_ptr<std::vector<std::string>> supers, 
				std::shared_ptr<std::vector<char>> opnames ) ;

};
#endif
