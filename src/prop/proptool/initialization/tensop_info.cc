#include <bagel_config.h>
#include  <src/prop/proptool/algebraic_manipulator/system_info.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//Build the operators here.
/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<TensOp::TensOp<DataType>> System_Info<DataType>::System_Info::Initialize_Tensor_Op_Info( string op_name ) {
/////////////////////////////////////////////////////////////////////////////////
cout << "shared_ptr<TensOp::TensOp<DataType>> System_Info<DataType>::System_Info::Initialize_Tensor_Op_Info" << endl;

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
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free });
   time_symm = "none";
   symmfuncs = identity_only();
   constraints = {  &System_Info<double>::System_Info::always_true };
   state_dep = 0;

  } else if ( op_name == "h" ) {  /* ---- h Tensor ( 1 electron Hamiltonian ) ----  */

   factor = (DataType) (1.0);
   idxs = make_shared<vector<string>>(vector<string> {"h0", "h1"});
   aops = make_shared<vector<bool>>(vector<bool>  {true, false});
   idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free });
   time_symm = "none";
   symmfuncs = set_1el_symmfuncs();
   constraints = {  &System_Info<double>::System_Info::always_true };
   state_dep = 0;
     
  } else if ( op_name == "L" ) {  /* ---- L Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"L0", "L1", "L2", "L3"});
    aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, not_virt, not_virt });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &System_Info<double>::System_Info::NotAllAct };
    state_dep = 2;
    
      
  } else if ( op_name == "M" ) {  /* ---- M Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"M0", "M1", "M2", "M3"});
    aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, core, core });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &System_Info<double>::System_Info::NotAllAct };
    state_dep = 2;
 
  } else if ( op_name == "N" ) {  /* ---- N Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string> {"N0", "N1", "N2", "N3"});
    aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false });
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, act, core });
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = { &System_Info<double>::System_Info::NotAllAct };
    state_dep = 2;
    
  } else if ( op_name == "T" ) {  /* ---- T Tensor ----  */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
    aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} );
    idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, act, not_virt });
    time_symm = "none";
    symmfuncs = set_2el_symmfuncs();
    constraints = {  &System_Info<double>::System_Info::NotAllAct };
    state_dep = 2;

  } else if ( op_name == "X" ) {

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>( vector<string> {"X0", "X1", "X2", "X3"} );
    aops = make_shared<vector<bool>>( vector<bool> { false, false, true, true } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, act, act } );
    time_symm = "none";
    symmfuncs = set_2el_symmfuncs();
    constraints = { &System_Info<double>::System_Info::NotAllAct };
    state_dep = 2;

   } else if ( op_name == "Z" ) { /* 2el test op */

    factor = (DataType) (1.0);
    idxs = make_shared<vector<string>>( vector<string> { "Z0", "Z1", "Z2", "Z3" } );
    aops = make_shared<vector<bool>>( vector<bool>  { true, true, false, false } );
    idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { virt, virt, act, core } );
    time_symm = "none";
    symmfuncs = identity_only();
    constraints = {  &System_Info<double>::System_Info::always_true };
    state_dep = 2;

  } else {
    
    throw runtime_error("Do not have in-built definition for operator \"" + op_name + "\", aborting !!!" ) ;

  }

  shared_ptr<TensOp::TensOp<DataType>> new_tens = Build_TensOp( op_name, idxs, aops, idx_ranges, symmfuncs, constraints, factor, time_symm, false, state_dep ) ;
  new_tens->get_ctrs_tens_ranges();

  return new_tens;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
