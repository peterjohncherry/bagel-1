#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include  <src/smith/wicktool/expression_info.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Info<DataType>::Expression_Info::Initialize_Tensor_Op_Info( string OpName ) {
/////////////////////////////////////////////////////////////////////////////////
 
  if ( T_map->find(OpName) != T_map->end() ) { 
 
   cout << "already initialized " << OpName << " info " << endl;
  
  } else { 

    if ( OpName == "S" ) { /* ---- H Tensor  ACTIVE ONLY FOR TESTING----  */
      cout << "setting info for Hact" << endl; 
      pair<double,double>                S_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         S_idxs = make_shared<vector<string>>(vector<string> {"S0", "S1", "S2", "S3"});
      shared_ptr<vector<bool>>           S_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
      shared_ptr<vector<vector<string>>> S_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act, act, act }); 
      shared_ptr<double>                 S_dummy_data;
      string                             S_TimeSymm = "none";
     
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> S_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  S_constraints = {  &Expression_Info<double>::Expression_Info::always_true };
     
      shared_ptr<TensOp<double>> STens = Build_TensOp("S", S_dummy_data, S_idxs, S_aops, S_idx_ranges, S_symmfuncs, S_constraints, S_factor, S_TimeSymm, false ) ;
      T_map->emplace("S", STens);
    
    } else if ( OpName == "U" ) {  /* ---- 4idx UNIT Tensor ACTIVE ONLY FOR TESTING----  */
    
      pair<double,double>                U_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         U_idxs = make_shared<vector<string>>(vector<string> {"U0", "U1", "U2", "U3"});
      shared_ptr<vector<bool>>           U_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
      shared_ptr<vector<vector<string>>> U_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> {  act, act, act, act }); 
      shared_ptr<double>                 U_dummy_data;
      string                             U_TimeSymm = "none";
     
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> U_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  U_constraints = {  &Expression_Info<double>::Expression_Info::always_true };
     
      shared_ptr<TensOp<double>> UTens = Build_TensOp( "U", U_dummy_data, U_idxs, U_aops, U_idx_ranges, U_symmfuncs, U_constraints, U_factor, U_TimeSymm, false ) ;
      T_map->emplace("U", UTens);
    
    } else if ( OpName == "R" ) {  /* ---- 6idx UNIT Tensor ACTIVE ONLY FOR TESTING----  */
    
      pair<double,double>                R_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         R_idxs = make_shared<vector<string>>(vector<string> {"R0", "R1", "R2", "R3","R4", "R5"} );
      shared_ptr<vector<bool>>           R_aops = make_shared<vector<bool>>(vector<bool>  {true, true, true, false, false, false} ); 
      shared_ptr<vector<vector<string>>> R_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> {  act, act, act, act, act, act }); 
      shared_ptr<double>                 R_dummy_data;
      string                             R_TimeRymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> R_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  R_constraints = {  &Expression_Info<double>::Expression_Info::always_true };
      
      shared_ptr<TensOp<double>> RTens = Build_TensOp("R", R_dummy_data, R_idxs, R_aops, R_idx_ranges, R_symmfuncs, R_constraints, R_factor, R_TimeRymm, false ) ;
      T_map->emplace("R", RTens);
    
    } else if ( OpName == "Q" ) {   /* ---- 8idx UNIT Tensor ACTIVE ONLY FOR TESTING----  */
    
      pair<double,double>                Q_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         Q_idxs = make_shared<vector<string>>(vector<string> {"Q0", "Q1", "Q2", "Q3","Q4", "Q5", "Q6", "Q7"});
      shared_ptr<vector<bool>>           Q_aops = make_shared<vector<bool>>(vector<bool>  {true, true, true, true, false, false, false, false}); 
      shared_ptr<vector<vector<string>>> Q_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> {  act, act, act, act, act, act, act, act }); 
      shared_ptr<double>                 Q_dummy_data;
      string                             Q_TimeSymm = "none";
   
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Q_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  Q_constraints = {  &Expression_Info<double>::Expression_Info::always_true };
   
      shared_ptr<TensOp<double>> QTens = Build_TensOp("Q", Q_dummy_data, Q_idxs, Q_aops, Q_idx_ranges, Q_symmfuncs, Q_constraints, Q_factor, Q_TimeSymm, false ) ;
      T_map->emplace("Q", QTens);
   
    
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
