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

  // NOTE!!: The last_iters vector contains points to the last element you want to actually use for comparison.
  //         It is almost certainly  _NEVER_ container.end(), far more likely container.end()-1 .
  template<typename IterType> 
  bool fvec_cycle_skipper_iter( std::vector<IterType>& forvec,
                                const std::vector<IterType>& last_iters,
                                const std::vector<IterType>& first_iters ) {
   
    typename std::vector<IterType>::const_reverse_iterator li_it = last_iters.crbegin();
    typename std::vector<IterType>::const_reverse_iterator fi_it = first_iters.crbegin();
    typename std::vector<IterType>::reverse_iterator id_1st_it = forvec.rend()-1;
    for( typename std::vector<IterType>::reverse_iterator  fv_it = forvec.rbegin(); fv_it !=forvec.rend(); fv_it++, li_it++ , fi_it++ ) {
      if ( *li_it == *fi_it ) {
        if ( fv_it == id_1st_it )
          return false;
      } else if (*fv_it == *li_it) {
        if ( fv_it == id_1st_it )
          return false;    
        else
        *fv_it = *fi_it;
      } else {
        ++(*fv_it);
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
  std::vector<vtype> inverse_reorder_vector(const std::vector<int>& reordering, const std::vector<vtype>& origvec ) {
    std::vector<int> inverse_reordering = get_ascending_order( reordering );
    std::vector<vtype> newvec(origvec.size());
    typename std::vector<vtype>::iterator nv_it = newvec.begin();
    for( std::vector<int>::const_iterator ir_it = inverse_reordering.begin(); ir_it != inverse_reordering.end(); nv_it++, ir_it++ )
      *nv_it = origvec[*ir_it];
    return newvec;
  };

  template<typename vtype>
  std::vector<vtype> reorder_vector(const std::vector<int>& reordering , const std::vector<vtype>& origvec ) {
    std::vector<vtype> newvec(origvec.size());
    typename std::vector<vtype>::iterator nv_it = newvec.begin();
    for( std::vector<int>::const_iterator r_it = reordering.begin(); r_it != reordering.end(); nv_it++, r_it++ )
      *nv_it = origvec[*r_it];
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

  template < typename DataType> 
  std::vector< DataType> get_subvector ( const std::vector<DataType>& main_vec, const std::vector<int>& sub_vec_pos ) { 
    std::vector<int>::const_iterator svp_it = sub_vec_pos.begin();
    std::vector<DataType> sub_vec(sub_vec_pos.size()) ;
    for ( typename std::vector<DataType>::iterator sv_it = sub_vec.begin(); sv_it != sub_vec.end() ; sv_it++, svp_it++ )
      *sv_it = main_vec[*svp_it];
    return sub_vec;
  };

template<typename ContType > 
void print_vec_of_conts( const std::vector<ContType>& vecvec, std::string name  ) {

   std::vector<int> maxs(vecvec.size());
   typename std::vector<ContType>::const_iterator vv_it = vecvec.begin();
   for ( std::vector<int>::iterator mx_it = maxs.begin(); mx_it!= maxs.end(); ++mx_it, ++vv_it )
     *mx_it =( vv_it->size() -1);

   std::vector<int> mins(vecvec.size(), 0);
   std::vector<int> fvec(vecvec.size(), 0);

   if (name != "" ){ 
     std::cout << name << " = "; std::cout.flush();
     std::string indent_s( name.size()+2, ' ');   
     bool first_block = true;
     do  {
       if ( first_block ) {
         first_block = false;
       } else {
         std::cout << indent_s ; std::cout.flush();
       }
       print_vector( get_1d_from_2d ( vecvec, fvec ) ); std::cout << std::endl;
     } while( fvec_cycle_skipper( fvec, maxs, mins) );

   
   } else { 
     do  {
       print_vector( get_1d_from_2d ( vecvec, fvec ) ); std::cout << std::endl;
     } while( fvec_cycle_skipper( fvec, maxs, mins) );

   } 
   return; 
}
};
#endif
