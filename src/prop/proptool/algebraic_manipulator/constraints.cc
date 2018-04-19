#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>
#include <iostream>
using namespace std;

bool Constraint_NotAllAct::apply_constraint( const std::vector<std::string>& rngs ){ 
  cout  << " Constraint_NotAllAct::apply_constraint" << endl;
  for ( vector<string>::const_iterator r_it = rngs.begin(); r_it != rngs.end() ; r_it++ ) 
    if ( *r_it != "a" && *r_it != "A" ) 
      return true;
  return false;
}
