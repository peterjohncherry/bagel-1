#ifndef __SRC_PROP_PROPTOOL_SYMMOPS_H
#define __SRC_PROP_PROPTOOL_SYMMOPS_H

#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include <map> 

namespace Symmetry_Operations {

      void hermitian( std::vector<std::string>& rngs); 
      void spin_flip( std::vector<std::string>& rngs);
      std::vector<std::string> f_1032( std::vector<std::string>& rngs); 
      std::vector<std::string> f_2143( std::vector<std::string>& rngs);  
      std::vector<std::string> f_2301( std::vector<std::string>& rngs);  
      std::vector<std::string> f_2103( std::vector<std::string>& rngs);  
      std::vector<std::string> f_3012( std::vector<std::string>& rngs);  
      std::vector<std::string> f_0321( std::vector<std::string>& rngs);  
      std::vector<std::string> f_1230( std::vector<std::string>& rngs);  

};

#endif
