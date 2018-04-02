#ifndef __SRC_PROP_PROPTOOL_algebra_utils_H
#define __SRC_PROP_PROPTOOL_algebra_utils_H
#include <iostream>     
#include <algorithm>    
#include <vector> 
#include <stdexcept>

namespace Algebra_Utils { 
  
  template<typename DataType>
  void transform_tens_vec( char transformation, std::vector<DataType>& invec  ) { 
     switch( transformation ) {
       case 'n' : return;
       case '0' : return;
       case 'h' : std::reverse( invec.begin(), invec.end() ); return;
       case 'H' : std::reverse( invec.begin(), invec.end() ); return;
       case 't' : std::reverse( invec.begin(), invec.end() ); return;
       case 'T' : std::reverse( invec.begin(), invec.end() ); return;
       default: 
           std::cout << "do not have transformation " << transformation << "implemented; please check the braket specification in the input file." << std::endl;
           throw std::logic_error( " Aborting!!" );
     } 
   }
  

  template<typename IterType>
  void transform_tens_vec( char transformation, IterType begin_it, IterType end_it  ) { 
     switch( transformation ) {
       case 'n' : return;
       case 'o' : return;
       case 'p' : return;
       case '0' : return;
       case '1' : return;
       case '2' : return;

       case 'h' : return;
       case 'H' : return;
       case 'i' : return;
       case 'I' : return;
       case 'j' : std::reverse( begin_it, end_it ); return;
       case 'J' : std::reverse( begin_it, end_it ); return;

       case 't' : std::reverse( begin_it, end_it ); return;
       case 'T' : std::reverse( begin_it, end_it ); return; 
       case 'u' : std::reverse( begin_it, end_it ); return; 
       case 'U' : std::reverse( begin_it, end_it ); return;
       case 'v' : return; 
       case 'V' : return;
       default: 
           std::cout << "do not have transformation " << transformation << "implemented; please check the braket specification in the input file." << std::endl;
           throw std::logic_error( " Aborting!!" );
     } 
  }
}

#endif
