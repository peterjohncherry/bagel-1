#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <algorithm>
using namespace std;

void Symmetry_Operations::hermitian( std::vector<std::string>& rngs) { 
  std::reverse( rngs.begin(), rngs.end());
  return;
}

void Symmetry_Operations::spin_flip( std::vector<std::string>& rngs ) { 
  for ( std::vector<std::string>::iterator r_it = rngs.begin(); r_it != rngs.end(); r_it++ ) {
    if ( (*r_it)[0] > 'Z' ) { 
      *r_it = toupper((*r_it)[0]);
    } else {                     
      *r_it = tolower((*r_it)[0]);
    } 
  }
  return;
}
std::vector<std::string> 
Symmetry_Operations::f_1032( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[1], rngs[0] , rngs[3], rngs[2] };
  return rngs_new;
}

std::vector<std::string> 
Symmetry_Operations::f_2143( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[2], rngs[1] , rngs[4], rngs[3] };
  return rngs_new;
}

std::vector<std::string> 
Symmetry_Operations::f_2301( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[2], rngs[3] , rngs[0], rngs[1] };
  return rngs_new;
}

std::vector<std::string> 
Symmetry_Operations::f_2103( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[2], rngs[1] , rngs[0], rngs[3] };
  return rngs_new;
}

std::vector<std::string> 
Symmetry_Operations::f_3012( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[3], rngs[0] , rngs[1], rngs[2] };
  return rngs_new;
}

std::vector<std::string> 
Symmetry_Operations::f_0321( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[0], rngs[3] , rngs[2], rngs[1] };
  return rngs_new;
}
std::vector<std::string> 
Symmetry_Operations::f_1230( std::vector<std::string>& rngs ) { 
  vector<std::string> rngs_new = { rngs[1], rngs[2] , rngs[3], rngs[0] };
  return rngs_new;
}
