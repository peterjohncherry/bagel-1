#ifndef __SRC_PROP_PROPTOOL_DEBUGGING_UTILS_H
#define __SRC_PROP_PROPTOOL_DEBUGGING_UTILS_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <memory>

namespace Debugging_Utils {  

  template<typename DataType>
  void print_sizes( const std::vector<DataType> vec, std::string name ) {
 
    if ( name != "" ) { 
     std::cout << name << " sizes : " ; std::cout.flush();
    }
 
    std::cout << "[ " ; std::cout.flush();
    for ( const auto& elem : vec ) { 
      std::cout << elem.size() << " " ; std::cout.flush();
    }
    std::cout << "]" ; std::cout.flush();

  }
 
  template<typename DataType>
  void print_sizes( const std::vector<std::shared_ptr<DataType>> vec, std::string name ) {
 
    if ( name != "" ) { 
     std::cout << name << " sizes : " ; std::cout.flush();
    }
 
    std::cout << "[ " ;std::cout.flush();
    for ( const auto& elem : vec ) { 
      std::cout << elem->size() << " " ; std::cout.flush();
    }
    std::cout << "]"; std::cout.flush();
  }


  template<typename DataType>
  void print_names( const std::vector<DataType> vec, std::string name ) {
 
    if ( name != "" ) { 
     std::cout << name << " sizes : " ; std::cout.flush();
    }
 
    std::cout << "[ " ; std::cout.flush();
    for ( const auto& elem : vec ) { 
      std::cout << elem.name() << " " ; std::cout.flush();
    }
    std::cout << "]" ; std::cout.flush();

  }
 
  template<typename DataType>
  void print_names( const std::vector<std::shared_ptr<DataType>> vec, std::string name ) {
 
    if ( name != "" ) { 
     std::cout << name << " names : " ; std::cout.flush();
    }
 
    std::cout << "[ " ;std::cout.flush();
    for ( const auto& elem : vec ) { 
      std::cout << elem->name() << " " ; std::cout.flush();
    }
    std::cout << "]"; std::cout.flush();
  }
 

 
  template<typename DataType>
  std::shared_ptr<std::vector<DataType>> get_pos_from_arg2( const std::vector<int>& pos, const std::vector<DataType>& vec2 ) {

    std::shared_ptr<std::vector<DataType>> selected_elems = std::make_shared<std::vector<DataType>>(pos.size());
    std::vector<int>::const_iterator p_it = pos.begin(); 
    for ( typename std::vector<DataType>::iterator se_it = selected_elems->begin(); se_it != selected_elems->end() ; se_it++, p_it++ )
      *se_it = vec2[*p_it]; 
  
    return selected_elems;
  } 
  
 
  template<typename DataType>
  DataType sum_unique_ptr_elems( std::unique_ptr<DataType[]>& my_data, const size_t my_size ) {

    DataType* end_ptr = my_data.get()+my_size;
    DataType total = (DataType)(0.0);
    for ( DataType* ptr = my_data.get(); ptr != end_ptr; ptr++ )
      total += *ptr;
   
    return total;
  } 
 
  template<typename DataType>
  void print_array( DataType* ptr , std::vector<size_t> dimensions, std::string name = " "  )  { 
   
    if ( name != " " )  
      std::cout << std::endl << "============== " << name << " =============== " << std::endl;

    size_t row_length = dimensions.back(); 
    size_t block_size = dimensions[dimensions.size()-2];
    DataType* ptr_buff = ptr;
  
    int num_blocks = 1;
    if ( dimensions.size() > 2 )
      for (int xx = 0 ; xx != dimensions.size()-2 ; xx++ ) 
        num_blocks *= dimensions[xx];
   
    for ( int rr = 0; rr != num_blocks ; rr++ ) {
      for ( int qq = 0 ; qq != block_size ; qq++ ) {
        for ( int ii = 0; ii != row_length ; ptr_buff++, ii++ ){ 
          std::cout << *ptr_buff << " "; std::cout.flush();
        } 
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl; 
  }
 
  template<typename DataType>
  void print_1d_array( DataType* ptr , int array_length, std::string name ) { 
    std::cout << name << " = [ "; std::cout.flush();
    for ( int qq = 0; qq != array_length; qq++, ptr++ )
      std::cout << *ptr  << " ";
    std::cout << " ]" << std::endl;
  }
 
  template<typename DataType1, typename DataType2>
  void vector_to_num_in_base( const std::vector<DataType1>& vec , DataType2 base = 10.0 ) { 

   DataType1 num = (DataType1)(0.0);
   DataType2 basis_elem = (DataType2)(1.0);
   for ( typename std::vector<DataType1>::const_reverse_iterator v_rit = vec.crbegin(); v_rit != vec.crend(); v_rit++ ) {
     num += basis_elem* (*v_rit);
     basis_elem*=base;
   }
   return num;
  }


}
#endif
