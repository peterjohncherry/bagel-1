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
    std::cout << "]" << std::endl;
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
    std::cout << "]" << std::endl;
  }
 


}
#endif
