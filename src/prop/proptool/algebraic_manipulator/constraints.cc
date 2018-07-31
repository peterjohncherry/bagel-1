#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>
#include <iostream>
#include <iterator>
#include <src/prop/proptool/proputils.h>


using namespace std;

#define __DEBUG_PROPTOOL_CONSTRAINTS

///////////////////////////////////////////////////////////////////////////////////////////
bool Constraint_NotAllAct::apply_constraint( const std::vector<std::string>& rngs ){ 
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
  cout  << " Constraint_NotAllAct::apply_constraint" << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////
 for ( vector<string>::const_iterator r_it = rngs.begin(); r_it != rngs.end() ; r_it++ )
    if ( *r_it != "a" && *r_it != "A" ) 
      return true;
  return false;
}
///////////////////////////////////////////////////////////////////////////////////////////
bool Constraint_Allowed_Blocks::apply_constraint( const std::vector<std::string>& rngs ){ 
#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
cout  << " Constraint_Allowed_Blocks::apply_constraint " << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////
  for ( const vector<std::string>& block  : allowed_blocks_ ){
    if ( block == rngs ){
      return true;
    }
    WickUtils::print_vector( rngs ); cout.flush(); cout << " != ";  WickUtils::print_vector( block ); cout << endl;
  }
#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
  WickUtils::print_vector( rngs, " does not satify constraints" ); cout << endl;
#endif
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//assumes normal ordering on input (Keep for now) 
bool Constraint_Spin_Neutral_Normal_Order::apply_constraint( const std::vector<std::string>& rngs){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
cout  << " Constraint_Spin_Neutral::apply_constraint" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////

  int spin_shift = 0;
  vector<string>::const_iterator r_it = rngs.begin();
  for ( int ii = 0; ii != rngs.size()/2 ; r_it++, ii++ ) 
    spin_shift =  (*r_it)[0] > 'Z'  ? spin_shift +1 : spin_shift -1 ;  

  for ( int ii = 0; ii != rngs.size()/2 ; r_it++, ii++ ) 
    spin_shift =  (*r_it)[0] < 'Z'  ? spin_shift +1 : spin_shift -1 ;  


#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
  cout << "spin_shift = " << spin_shift << " "; cout.flush(); WickUtils::print_vector(rngs, "rngs" ); cout << endl; 
#endif

  return ( spin_shift == 0 );
}

////////////////////////////////////////////////////////////////////////////////////////
//assumes normal ordering on input (Keep for now) 
bool Constraint_All_Same_Spin::apply_constraint( const std::vector<std::string>& rngs){ 
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CONSTRAINTS
cout  << " Constraint_All_Same_Spin::apply_constraint" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////

  vector<string>::const_iterator r_it = rngs.begin(); 
  if ( (*r_it)[0] > 'Z' ) {  
    ++r_it;
    for ( ; r_it != rngs.end() ; r_it++ ) 
      if ( ! ((*r_it)[0] > 'Z') )
        return false;

  } else {
    ++r_it;
    for ( ; r_it != rngs.end() ; r_it++ ) 
      if ( ! ((*r_it)[0] < 'Z') )
        return false;
  }
  return true;
}

//bool Constraint_Spin_Neutral::apply_constraint( const std::vector<std::string>& rngs, const std::vector<bool>& aops){ 
//  cout  << " Constraint_Spin_Neutral::apply_constraint" << endl;
//  vector<bool>::const_iterator aops_it = aops.begin();
//  for ( vector<string>::const_iterator r_it = rngs.begin(); r_it != rngs.end() ; r_it++, aops_it++ ) 
//     if ( *aops_it  )
//       spin_shift =  (*r_it)[0] > 'Z'  ? spin_shift -1 : spin_shift +1 ;  
//  return false;
//}
