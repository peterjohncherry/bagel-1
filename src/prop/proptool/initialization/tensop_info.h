#ifndef __SRC_PROPTOOL_TENSOP_INFO_INIT
#define __SRC_PROPTOOL_TENSOP_INFO_INIT
#include <src/global.h>
#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>

using namespace std;

namespace  TensOp_Info_Init {

/////////////////////////////////////////////////////////////////////////////////
//Build the operators here.
/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<TensOp::TensOp<DataType>> Initialize_Tensor_Op_Info( string op_name, shared_ptr<map<char, long unsigned int>> range_prime_map ) {
/////////////////////////////////////////////////////////////////////////////////
cout << "shared_ptr<TensOp::TensOp<DataType>>::Initialize_Tensor_Op_Info" << endl;


  static vector<string> free     = {"c", "C", "a", "A", "v", "V"};
  static vector<string> not_core = {"a", "A", "v", "V"};
  static vector<string> not_act  = {"c", "C", "v", "V"};
  static vector<string> not_virt = {"c", "C", "a", "A"};
  static vector<string> core     = {"c", "C"};
  static vector<string> act      = {"a", "A"};
  static vector<string> virt     = {"v", "V"};
 
  static vector<string> free_a     = {"c", "a", "v"};
  static vector<string> not_core_a = {"a", "v"};
  static vector<string> not_act_a  = {"c", "v"};
  static vector<string> not_virt_a = {"c", "a"};
  static vector<string> core_a     = {"c"};
  static vector<string> act_a      = {"a"};
  static vector<string> virt_a     = {"v"};
 
  static vector<string> free_b     = {"C", "A", "V"};
  static vector<string> not_core_b = {"A", "V"};
  static vector<string> not_act_b  = {"C", "V"};
  static vector<string> not_virt_b = {"C", "A"};
  static vector<string> core_b     = {"C"};
  static vector<string> act_b      = {"A"};
  static vector<string> virt_b     = {"V"};

  pair<double,double>                factor = make_pair(1.0, 0.0);
  shared_ptr<vector<string>>         idxs;
  shared_ptr<vector<bool>>           aops;
  shared_ptr<vector<vector<string>>> idx_ranges;
  string                             time_symm;
  vector<shared_ptr<Transformation>> symmfuncs = vector<shared_ptr<Transformation>>(0);
  vector<shared_ptr<Constraint>> constraints = vector<shared_ptr<Constraint>>(0);
  int state_dep;

  static shared_ptr<Transformation_Hermitian> hconj = make_shared<Transformation_Hermitian>( "hconj" );  
  static shared_ptr<Transformation_Spinflip>  spinflip = make_shared<Transformation_Spinflip>( "spinflip" );
  static shared_ptr<Transformation_1032>  perm_1032 = make_shared<Transformation_1032>( "1032" ); 
  static shared_ptr<Transformation_2301>  perm_2301 = make_shared<Transformation_2301>( "2301" );
  static shared_ptr<Transformation_2103>  perm_2103 = make_shared<Transformation_2103>( "2103" );
  static shared_ptr<Transformation_3012>  perm_3012 = make_shared<Transformation_3012>( "3012" );
  static shared_ptr<Transformation_0321>  perm_0321 = make_shared<Transformation_0321>( "0321" );
  static shared_ptr<Transformation_1230>  perm_1230 = make_shared<Transformation_1230>( "1230" );
  static shared_ptr<Transformation_ID>  identity = make_shared<Transformation_ID>( "Id" );

  static shared_ptr<Constraint_NotAllAct>  not_all_act = make_shared<Constraint_NotAllAct>();
  static shared_ptr<Constraint_Spin_Neutral_Normal_Order>  spin_neutral = make_shared<Constraint_Spin_Neutral_Normal_Order>();
  static shared_ptr<Constraint_All_Same_Spin>  all_same_spin = make_shared<Constraint_All_Same_Spin>();

  if ( op_name == "H" ) {  /* ---- H Tensor (2 electron Hamiltonian ----  */

   idxs = make_shared<vector<string>>(vector<string> { "H0", "H1", "H2", "H3" } );
   aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false } );//TODO check this ordering is correct
//   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, core, core });
//   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt_a, virt_a, core_a, core_a });
//   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt_a, virt_a, core_a, act_a });
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt_a, virt_a, act_a, core_a });

   //idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free, free, free, free });
//  symmfuncs = { hconj, perm_1032, perm_2301, perm_2103, perm_3012, perm_0321, perm_1230 };
//  constraints = { all_same_spin };
   time_symm = "none";
   state_dep = 0;
   factor = make_pair( -0.5, 0.0);

  } else if ( op_name == "S" ) {  /* ---- S Tensor ----  */

    idxs = make_shared<vector<string>>(vector<string>{"S0", "S1", "S2", "S3"}  );
    aops = make_shared<vector<bool>>  (vector<bool>  { true, true, false, false } );
    constraints = { not_all_act }; 
//    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { act_a, core_a, virt_a, virt_a } );
//    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { core_a, core_a, virt_a, virt_a } );
    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { act_a, core_a, virt_a, virt_a } );
    time_symm = "none";
    state_dep = 0;

 
  //} else if ( op_name == "Q" ) {  /*TEST*/
  //
  // idxs = make_shared<vector<string>>(vector<string> { "H0", "H1", "H2", "H3" } );
  // aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false } );//TODO check this ordering is correct
  // idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, virt, act, act });
  // constraints = { all_same_spin };
  // time_symm = "none";
  // state_dep = 0;
  //
  } else if ( op_name == "h" ) {  /* ---- h Tensor ( 1 electron Hamiltonian ) ----  */

   idxs = make_shared<vector<string>>(vector<string> {"h0", "h1"});
   aops = make_shared<vector<bool>>(vector<bool>  {true, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free });
   time_symm = "none";
   state_dep = 0;
  
  } else if ( op_name == "f" ) {  /* ---- state averaged fock operator ----  */

   idxs = make_shared<vector<string>>(vector<string> {"f0", "f1"});
   aops = make_shared<vector<bool>>(vector<bool>  {true, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free });
   time_symm = "none";
   state_dep = 0;

  } else if ( op_name == "T" ) {  /* ---- T Tensor ----  */

    idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
    aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} );
    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });
    //idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { virt_a, virt_a, act_a, act_a } );
    time_symm = "none";
    state_dep = 0;
   } else if ( op_name == "X" ) {

    idxs = make_shared<vector<string>>( vector<string> {"X3", "X2", "X1", "X0"} );
    aops = make_shared<vector<bool>>( vector<bool> { false, false, true, true } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt } );
    time_symm = "none";
    state_dep = 0;

  } else {
    
    throw runtime_error("Do not have in-built definition for operator \"" + op_name + "\", aborting !!!" ) ;

  }

  cout << "initializing tensop" << endl;
  shared_ptr<TensOp::TensOp<DataType>> new_tens =  make_shared<TensOp::TensOp<DataType>>( op_name, *idxs, *idx_ranges, *aops,
                                                                                          factor, symmfuncs, constraints, time_symm, state_dep, range_prime_map);
  return new_tens;
}
}
#endif
