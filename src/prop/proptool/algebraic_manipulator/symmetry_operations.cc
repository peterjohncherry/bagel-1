#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <algorithm>
using namespace std;

std::vector<std::string> 
Transformation_Hermitian::transform( const std::vector<std::string>& rngs) { 
  vector<std::string> rngs_new = rngs;
  std::reverse( rngs_new.begin(), rngs_new.end());
  return rngs_new;
}

std::vector<std::string> 
Transformation_Spinflip::transform( const std::vector<std::string>& rngs ) { 
  vector<string> rngs_new = rngs;
  for ( std::vector<std::string>::iterator rn_it = rngs_new.begin(); rn_it != rngs_new.end(); rn_it++ ) {
    if ( (*rn_it)[0] > 'Z' ) { 
      *rn_it = toupper((*rn_it)[0]);
    } else {                     
      *rn_it = tolower((*rn_it)[0]);
    } 
  }
  return rngs_new;
}

std::vector<std::string> 
Transformation_1032::transform( const std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[1], rngs[0] , rngs[3], rngs[2] };
  return rngs_new;
}

std::vector<std::string> 
Transformation_2103::transform( const std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[2], rngs[1] , rngs[0], rngs[3] };
  return rngs_new;
}

std::vector<std::string> 
Transformation_2301::transform( const std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[2], rngs[3] , rngs[0], rngs[1] };
  return rngs_new;
}

std::vector<std::string> 
Transformation_3012::transform( const std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[3], rngs[0] , rngs[1], rngs[2] };
  return rngs_new;
}

std::vector<std::string> 
Transformation_0321::transform( const std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[0], rngs[3] , rngs[2], rngs[1] };
  return rngs_new;
}

std::vector<std::string> 
Transformation_1230::transform( const std::vector<std::string>& rngs ) {
  vector<std::string> rngs_new = { rngs[1], rngs[2] , rngs[3], rngs[0] };
  return rngs_new;
}
