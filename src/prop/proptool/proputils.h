#ifndef __SRC_PROP_PROPTOOL_PROPUTILS_H
#define __SRC_PROP_PROPTOOL_PROPUTILS_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>
#include <string>
#include <list>
#include <memory>
#include <map>
#include <iostream>
#include <numeric>
#include <complex> 
#include <cassert>
#include <functional>
#include <cctype>
namespace WickUtils {  
 
  void pair_fac_mult( const std::pair<double,double>& factor_fixed , std::pair<double,double>& factor_changing );

  std::string get_civec_name( const int state_num,  const int norb,  const int nalpha, const int nbeta);

  std::string get_det_name( char name_range1, int num_range1, char name_range2, int num_range2, int norb );

  std::string get_gamma_name( const std::vector<std::string>& full_idx_ranges, const std::vector<bool>& aops_vec,
                              const std::vector<int>& idxs_pos, std::string Bra_name, std::string Ket_name );

  std::string get_ctp_name( const std::string op_state_name, const std::vector<std::string>& idxs, const std::vector<std::string>& id_ranges,
                            const std::vector<std::pair<int,int>>& ctrs_pos );

  bool RangeCheck(const std::vector<std::string>& id_ranges, const std::vector<bool>& aops ) ;
  
  std::shared_ptr<std::vector<std::pair<int,int>>>
  standardize_delta_ordering_generic( const std::vector<std::pair<int,int>>& deltas_pos, const std::vector<std::string>& idxs  );

  unsigned int range_to_prime(char range );

  unsigned int range_to_prime_spinfree(char range );

  std::vector<char> strvec_to_chrvec( const std::vector<std::string>& strvec );
 
  std::vector<std::string> chrvec_to_strvec( const std::vector<char>& chrvec );

  std::vector<int> get_ascending_order( const std::vector<int>& scrambled_vec );

  template<typename DataType> 
  bool fvec_cycle_skipper( std::vector<DataType>& forvec, const std::vector<DataType>& max, const std::vector<DataType>& min ) {
    for(int ii = forvec.size()-1; ii!=-1 ; ii--) {
      if ( max[ii] == min[ii] ) {
        if ( ii == 0 )
          return false;
      } else if (forvec[ii] == max[ii]) {
        if (ii == 0) 
          return false;    
        forvec[ii] = min[ii];
      } else {
        forvec[ii] = forvec[ii]+ 1;
        break;
      }
    }
    return true;
  }

  template<typename DataType> 
  bool fvec_cycle_skipper_f2b( std::vector<DataType>& forvec, const std::vector<DataType>& max , const std::vector<DataType>& min ) {
    for( int ii = 0; ii != forvec.size(); ii++) {
      if ( max[ii] == min[ii] ) {
        if ( ii == forvec.size()-1 )
          return false;
      } else if (forvec[ii] == max[ii]) {
        if (ii == forvec.size()-1) 
          return false;    
        forvec[ii] = min[ii];
      } else {
        forvec[ii] = forvec[ii]+ 1;
        break;
      }
    }
    return true;
  }
  template<class DataType>
  void print_vector(std::vector<DataType> invec, std::string name =""){
    if (name != "" ) 
      std::cout << name << " ="; std::cout.flush();
    std::cout << " [ ";
    for (auto  elem : invec)
      std::cout << elem << " " ; std::cout.flush();
    std::cout << "]  " ;
    return;
  }

  template<class T1, class T2  >
  void print_pair_vector(std::vector<std::pair<T1,T2>> invec, std::string name =""){
    if (name != "" ) 
      std::cout << name << " =";  std::cout.flush();
    std::cout << " [ ";
    for (auto  elem : invec)
      std::cout << "(" << elem.first << "," << elem.second << ") ";
    std::cout << "]  " ; std::cout.flush();
    return;
  }

  template<class T1, class T2, class T3, class T4 >
  void print_pair_pair_vector( std::vector<std::pair<std::pair<T1,T2>, std::pair<T3,T4>>> ccp_vec,
                               std::string name = "" ) {
    std::cout << name << " = [ "; std::cout.flush();
    for ( auto& elem : ccp_vec ){
      std::cout << "{(" << elem.first.first << "," << elem.first.second << "),("  ;std::cout.flush();
      std::cout << elem.second.first << "," << elem.second.second << ")} "  ; std::cout.flush();
    }
    std::cout << "]" ; std::cout.flush();
    
    return;
  }

  template<typename DataType> // returns the relative order of the elements in destination and origin
  std::shared_ptr<std::vector<int>> get_pattern_match_order( std::shared_ptr<std::vector<DataType>> destination , std::shared_ptr<std::vector<DataType>> origin ) { 
    assert( destination->size() == origin->size() ); 
    std::shared_ptr<std::vector<int>> new_order = std::make_shared<std::vector<int>>(destination->size()); 
    std::vector<int>::iterator new_order_it = new_order->begin();  
    for ( int ii = 0 ; ii != destination->size() ; ii++  ) 
      for ( int jj = 0 ; jj != origin->size() ; jj++  ) 
        if ( destination->at(ii) == origin->at(jj) ) {
          *new_order_it++ = jj;
          break;
        }
    assert( new_order_it == new_order->end()); // if this trips the vectors don't have the same elements 
    return new_order;
  }

  template<typename vtype>
  std::vector<vtype> inverse_reorder_vector(const std::vector<int>& neworder, const std::vector<vtype>& origvec ) {
    std::vector<int> new_order = get_ascending_order( new_order );
    std::vector<vtype> newvec(origvec.size());
    typename std::vector<vtype>::iterator nv_it = newvec.begin();
    for( std::vector<int>::const_iterator no_it = neworder.begin(); no_it != neworder.end(); nv_it++, no_it++ )
      *nv_it = origvec[*no_it];
    return newvec;
  };

  template<typename vtype>
  std::vector<vtype> reorder_vector(const std::vector<int>& neworder , const std::vector<vtype>& origvec ) {
    std::vector<vtype> newvec(origvec.size());
    typename std::vector<vtype>::iterator nv_it = newvec.begin();
    for( std::vector<int>::const_iterator no_it = neworder.begin(); no_it != neworder.end(); nv_it++, no_it++ )
      *nv_it = origvec[*no_it];
    return newvec;
  };

  template< typename DataType>  
  void reorder_vector_inplace(const std::vector<int>& new_order, std::vector<DataType>& orig_vec ){
    std::vector<DataType> reordered_vec(orig_vec.size());
    auto rv_it = reordered_vec.begin();
    for( std::vector<int>::const_iterator no_it = new_order.begin(); no_it != new_order.end(); no_it++, rv_it++ )
       *rv_it = orig_vec[*no_it];
    orig_vec =  std::move(reordered_vec); 
    return;
  };

  template < typename DataType> 
  std::vector< DataType> get_1d_from_2d ( const std::vector<std::vector<DataType>>&  source_2d, const std::vector<int>& pos_vec ) { 
    std::vector<DataType> out_vec( pos_vec.size() );
    typename std::vector<DataType>::iterator ov_it = out_vec.begin(); 
    typename std::vector<std::vector<DataType>>::const_iterator s2_it = source_2d.begin(); 
    for ( std::vector<int>::const_iterator pv_it = pos_vec.begin(); pv_it != pos_vec.end() ; ++pv_it, ++ov_it, ++s2_it ) 
      *ov_it = (*s2_it)[*pv_it];
    return out_vec;      
  };
};
#endif
