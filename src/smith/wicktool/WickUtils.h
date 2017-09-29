#ifndef __SRC_SMITH_WickUtils_H
#define __SRC_SMITH_WickUtils_H

#include<stdio.h>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include<utility>
#include<tuple>
#include<string>
#include<memory>
#include<map>
#include <iostream>
#include <numeric>
 
namespace WickUtils {  
  using delta_ints = std::vector<std::vector<std::pair<int,int>>>;
  using delta_strs = std::vector<std::vector<std::pair<std::string,std::string>>>;
  using delta_bools = std::vector<std::vector<std::pair<bool,bool>>>;

  using vv_ints  = std::vector< std::vector<int> >;
  using vv_strs  = std::vector< std::vector<std::string> >;
  using vv_bools = std::vector< std::vector<bool> >;

  using pint_vec = std::vector<std::pair<int,int>>;
  using pstr_vec = std::vector<std::pair<std::string,std::string>>;
  using pbool_vec = std::vector<std::pair<bool,bool>>;
 
  //routines for mimicking for loop of arbitrary depth
  void fvec_cycle(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<int>> max ) ;
  bool fvec_cycle(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<int>> max , std::shared_ptr<std::vector<int>> min) ;
  bool constrained_fvec_cycle(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<int>> max) ;

  template<class T1>
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<T1>>>> combgen( std::shared_ptr<std::vector<T1>> invec);

  std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> get_N_in_M_combsX( std::shared_ptr<std::vector<int>> vec1, int NN );

  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<int,int>>>>>
  get_unique_pairs(std::shared_ptr<std::vector<int>> ids1 , std::shared_ptr<std::vector<int>> ids2 , int num_pairs);

  template<class DT>
  void print_vec(std::vector<DT> invec , std::string vecname);
 
  void print_pvec (pint_vec pvec) ;

  std::string get_name(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_id_ranges,  std::shared_ptr<pint_vec> all_ctrs_pos) ;

  std::string get_name_rdm(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges,
                           std::shared_ptr<std::vector<std::pair<int,int>>> all_deltas_pos ); 


  std::shared_ptr<std::vector<int>> get_unc_ids_from_deltas_ids_comparison(std::shared_ptr<std::vector<int>> ids , std::shared_ptr<std::vector<std::pair<int,int>>> deltas );

  std::shared_ptr<std::vector<int>> get_unc_ids_from_deltas_ids(std::shared_ptr<std::vector<int>> ids , std::shared_ptr<std::vector<std::pair<int,int>>> deltas );

  std::string get_Aname(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges, 
                        std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos );

  std::string get_gamma_name(std::shared_ptr<std::vector<std::string>> full_idx_ranges,  std::shared_ptr<std::vector<bool>> aops_vec,
                             std::shared_ptr<std::vector<int>> idxs_pos );

  std::string get_gamma_name(std::shared_ptr<std::vector<bool>> aops_vec, std::shared_ptr<std::vector<std::string>> full_idx_ranges,  
                             std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos );

  
  bool fvec_cycle_skipper(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<int>> max,
                          std::shared_ptr<std::vector<int>> min,  int fixed_index_position );
 
  bool fvec_cycle_skipper(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<int>> max, int fixed_index_position );

  template<class DataType>
  void print_vecX(std::vector<DataType> invec, std::string name ="lazy"){
    std::cout << name <<" = [ ";
    for (auto  elem : invec)
      std::cout << elem << " " ;
    std::cout << "]  " ;
    return;
  }


  template<class DataType>
  void print_pairvec(std::vector<std::pair<DataType,DataType>> invec, std::string name ="lazy"){
    std::cout << name <<" = ( ";
    for (auto  elem : invec)
      std::cout << "(" << elem.first << "," << elem.second << ") ";
    std::cout << ")  " ;
    return;
  }

}

#endif
