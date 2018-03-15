#ifndef __SRC_PROPTOOL_TENSOP_INFO_INIT
#define __SRC_PROPTOOL_TENSOP_INFO_INIT
#include <src/global.h>
#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>

using namespace std;
using namespace Symmetry_Operations;
namespace  TensOp_Info_Init {

/////////////////////////////////////////////////////////////////////////////////
//Build the operators here.
/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<TensOp::TensOp<DataType>> Initialize_Tensor_Op_Info( string op_name, shared_ptr<map<char, long unsigned int>> range_prime_map ) {
/////////////////////////////////////////////////////////////////////////////////
cout << "shared_ptr<TensOp::TensOp<DataType>> System_Info<DataType>::System_Info::Initialize_Tensor_Op_Info" << endl;

  vector<string> free     = {"cor", "act", "vir"};
  vector<string> not_core = {"act", "vir"};
  vector<string> not_act  = {"cor", "vir"};
  vector<string> not_virt = {"cor", "act"};
  vector<string> core     = {"cor"};
  vector<string> act      = {"act"};
  vector<string> virt     = {"vir"};

  DataType                           factor;
  shared_ptr<vector<string>>         idxs;
  shared_ptr<vector<bool>>           aops;
  shared_ptr<vector<vector<string>>> idx_ranges;
  string                             time_symm;
  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> symmfuncs;
  vector<bool(*)(shared_ptr<vector<string>>)>  constraints;
  int state_dep;


  if ( op_name == "H" ) {  /* ---- H Tensor (2 electron Hamiltonian ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"H0", "H1", "H2", "H3"});
   aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false});//TODO check this ordering is correct
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free, free, free, free });
   time_symm = "none";
   symmfuncs = identity_only();
   constraints = {  &Symmetry_Operations::always_true };
   state_dep = 0;

  } else if ( op_name == "h" ) {  /* ---- h Tensor ( 1 electron Hamiltonian ) ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"h0", "h1"});
   aops = make_shared<vector<bool>>(vector<bool>  {true, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free });
   time_symm = "none";
   symmfuncs = identity_only(); //set_1el_symmfuncs();
   constraints = {  &Symmetry_Operations::always_true };
   state_dep = 0;
  
  } else if ( op_name == "Q" ) {  /* ---- test six index ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"Q0", "Q1", "Q2", "Q3", "Q4", "Q5" });
   aops = make_shared<vector<bool>>(vector<bool>  {true, true, true, false, false, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_core,  not_core, not_core, not_virt, not_virt, not_virt });
   time_symm = "none";
   symmfuncs = identity_only();
   constraints = {  &Symmetry_Operations::always_true };
   state_dep = 0;
 
  } else if ( op_name == "R" ) {  /* ---- test six index ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"R0", "R1", "R2", "R3", "R4", "R5" });
   aops = make_shared<vector<bool>>(vector<bool>  {true, true, true, false, false, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt,  virt, act, act, core, core });
   time_symm = "none";
   symmfuncs = identity_only();
   constraints = {  &Symmetry_Operations::always_true };
   state_dep = 0;

  } else if ( op_name == "f" ) {  /* ---- state averaged fock operator ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"f0", "f1"});
   aops = make_shared<vector<bool>>(vector<bool>  {true, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free });
   time_symm = "none";
   symmfuncs = set_1el_symmfuncs();
   constraints = {  &Symmetry_Operations::always_true };
   state_dep = 0;

  } else if ( op_name == "L" ) {  /* ---- L Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"L0", "L1", "L2", "L3"});
    aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, not_virt, not_virt });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &Symmetry_Operations::NotAllAct };
    state_dep = 2;
    
      
  } else if ( op_name == "M" ) {  /* ---- M Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"M0", "M1", "M2", "M3"});
    aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, core, core });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &Symmetry_Operations::NotAllAct };
    state_dep = 2;
 
  } else if ( op_name == "N" ) {  /* ---- N Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"N0", "N1", "N2", "N3"});
    aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act, act, act });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &Symmetry_Operations::always_true  };
    state_dep = 0;
    
  } else if ( op_name == "T" ) {  /* ---- T Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
    aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} );
    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = {  &Symmetry_Operations::NotAllAct };
    state_dep = 2;

  } else if ( op_name == "t" ) {  /* ---- T Tensor herm conj TODO  should find a better way fo dealing with this----  */
    cout << "getting t op " << endl;
    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string>{"t0", "t1", "t2", "t3"}  );
    aops = make_shared<vector<bool>>  (vector<bool>  {false, false, true, true } );
    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &Symmetry_Operations::NotAllAct };
    state_dep = 2;

  } else if ( op_name == "X" ) {

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>( vector<string> {"X0", "X1", "X2", "X3"} );
    aops = make_shared<vector<bool>>( vector<bool> { false, false, true, true } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt } );
    time_symm = "none";
    symmfuncs = identity_only();
   //constraints = {  &Symmetry_Operations::always_true };
    constraints = { &Symmetry_Operations::NotAllAct };
    state_dep = 0;

  } else if ( op_name == "x" ) {

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>( vector<string> {"X0", "X1"} );
    aops = make_shared<vector<bool>>( vector<bool> { false, true } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act } );
    time_symm = "none";
    symmfuncs = identity_only();
   //constraints = {  &Symmetry_Operations::always_true };
    constraints = { &Symmetry_Operations::NotAllAct };
    state_dep = 2;

   } else if ( op_name == "Z" ) { /* 2el test op */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>( vector<string> { "Z0", "Z1", "Z2", "Z3" } );
    aops = make_shared<vector<bool>>( vector<bool>  { false, false, true, true } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act, act, act } );
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = {  &Symmetry_Operations::always_true };
    state_dep = 0;

  } else {
    
    throw runtime_error("Do not have in-built definition for operator \"" + op_name + "\", aborting !!!" ) ;

  }

  shared_ptr<TensOp::TensOp<DataType>> new_tens =  make_shared<TensOp::TensOp<DataType>>( op_name, *idxs, *idx_ranges, *aops,
                                                                                          factor, symmfuncs, constraints, time_symm, state_dep, range_prime_map);
  new_tens->get_ctrs_tens_ranges();

  return new_tens;
}
}
#endif
