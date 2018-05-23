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
 


}
#endif
