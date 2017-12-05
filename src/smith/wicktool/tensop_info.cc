#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include  <src/smith/wicktool/system_info.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void System_Info<DataType>::System_Info::Initialize_Tensor_Op_Info( string OpName ) {
/////////////////////////////////////////////////////////////////////////////////
 
  if ( T_map->find(OpName) != T_map->end() ) { 
 
   cout << "already initialized " << OpName << " info " << endl;
  
  } else { 

    if ( OpName == "S" ) { /* ---- H Tensor  ACTIVE ONLY FOR TESTING----  */
      cout << "setting info for Hact" << endl; 
      pair<double,double>                S_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         S_idxs = make_shared<vector<string>>(vector<string> {"S0", "S1", "S2", "S3"});
      shared_ptr<vector<bool>>           S_aops = make_shared<vector<bool>>( vector<bool>  {true, true, false, false} ); 
      shared_ptr<vector<vector<string>>> S_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act, act, act }); 
      shared_ptr<double>                 S_dummy_data;
      string                             S_TimeSymm = "none";
     
//      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> S_symmfuncs = set_2el_symmfuncs();
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> S_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  S_constraints = {  &System_Info<double>::System_Info::always_true };
     
      shared_ptr<TensOp<double>> STens = Build_TensOp("S", S_dummy_data, S_idxs, S_aops, S_idx_ranges, S_symmfuncs, S_constraints, S_factor, S_TimeSymm, false ) ;
      T_map->emplace("S", STens);
 
   }else if ( OpName == "Z" ) { /* ---- H Tensor  ACTIVE ONLY FOR TEZTING----  */
      cout << "setting info for Hact" << endl; 
      pair<double,double>                Z_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         Z_idxs = make_shared<vector<string>>(vector<string> {"Z0", "Z1", "Z2", "Z3"});
      shared_ptr<vector<bool>>           Z_aops = make_shared<vector<bool>>( vector<bool>  {true, false, true, false} ); 
      shared_ptr<vector<vector<string>>> Z_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { act, act, act, act }); 
      shared_ptr<double>                 Z_dummy_data;
      string                             Z_TimeSymm = "none";
     
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Z_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  Z_constraints = {  &System_Info<double>::System_Info::always_true };
     
      shared_ptr<TensOp<double>> ZTens = Build_TensOp("Z", Z_dummy_data, Z_idxs, Z_aops, Z_idx_ranges, Z_symmfuncs, Z_constraints, Z_factor, Z_TimeSymm, false ) ;
      T_map->emplace("Z", ZTens);
   
    } else if ( OpName == "U" ) {  /* ---- 4idx UNIT Tensor ACTIVE ONLY FOR TESTING----  */
    
      pair<double,double>                U_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         U_idxs = make_shared<vector<string>>(vector<string> {"U0", "U1", "U2", "U3"});
      shared_ptr<vector<bool>>           U_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
      shared_ptr<vector<vector<string>>> U_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> {  act, act, act, act }); 
      shared_ptr<double>                 U_dummy_data;
      string                             U_TimeSymm = "none";
     
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> U_symmfuncs = identity_only();
      vector<bool(*)(shared_ptr<vector<string>>)>  U_constraints = {  &System_Info<double>::System_Info::always_true };
     
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
      vector<bool(*)(shared_ptr<vector<string>>)>  R_constraints = {  &System_Info<double>::System_Info::always_true };
      
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
      vector<bool(*)(shared_ptr<vector<string>>)>  Q_constraints = {  &System_Info<double>::System_Info::always_true };
   
      shared_ptr<TensOp<double>> QTens = Build_TensOp("Q", Q_dummy_data, Q_idxs, Q_aops, Q_idx_ranges, Q_symmfuncs, Q_constraints, Q_factor, Q_TimeSymm, false ) ;
      T_map->emplace("Q", QTens);
 
    } else if ( OpName == "P" ) {  /* ---- test tensor with general H2el shape ----  */

      pair<double,double>                P_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         P_idxs = make_shared<vector<string>>(vector<string> {"P0", "P1", "P2", "P3"});
      shared_ptr<vector<bool>>           P_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false});//TODO check this ordering is correct 
      shared_ptr<vector<vector<string>>> P_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
      shared_ptr<double>                 P_dummy_data;
      string                             P_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> P_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  P_constraints = {  &System_Info<double>::System_Info::always_true };
      
      shared_ptr<TensOp<double>> PTens = Build_TensOp("P", P_dummy_data, P_idxs, P_aops, P_idx_ranges, P_symmfuncs, P_constraints, P_factor, P_TimeSymm, false ) ;
      T_map->emplace("P", PTens);
  
    } else if ( OpName == "H" ) {  /* ---- H Tensor (2 electron Hamiltonian ----  */

      pair<double,double>                H_factor = make_pair(0.5,0.5);
      shared_ptr<vector<string>>         H_idxs = make_shared<vector<string>>(vector<string> {"H0", "H1", "H2", "H3"});
      shared_ptr<vector<bool>>           H_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false});//TODO check this ordering is correct 
      shared_ptr<vector<vector<string>>> H_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
      shared_ptr<double>                 H_dummy_data;
      string                             H_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> H_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  H_constraints = {  &System_Info<double>::System_Info::always_true };
      
      shared_ptr<TensOp<double>> HTens = Build_TensOp("H", H_dummy_data, H_idxs, H_aops, H_idx_ranges, H_symmfuncs, H_constraints, H_factor, H_TimeSymm, false ) ;
      T_map->emplace("H", HTens);
     } else if ( OpName == "h" ) {  /* ---- h Tensor ( 1 electron Hamiltonian ) ----  */

      pair<double,double>                h_factor = make_pair(1.0,1.0);
      shared_ptr<vector<string>>         h_idxs = make_shared<vector<string>>(vector<string> {"h0", "h1"});
      shared_ptr<vector<bool>>           h_aops = make_shared<vector<bool>>(vector<bool>  {true, false}); 
      shared_ptr<vector<vector<string>>> h_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free }); 
      shared_ptr<double>                 h_dummy_data;
      string                             h_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> h_symmfuncs = set_1el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  h_constraints = {  &System_Info<double>::System_Info::always_true };
      
      shared_ptr<TensOp<double>> hTens = Build_TensOp("h", h_dummy_data, h_idxs, h_aops, h_idx_ranges, h_symmfuncs, h_constraints, h_factor, h_TimeSymm, false ) ;
      T_map->emplace("h", hTens);
 
    
    } else if ( OpName == "T" ) {  /* ---- T Tensor ----  */

      pair<double,double>                 T_factor = make_pair(1.0,1.0);
      shared_ptr<vector<string>>          T_idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
      shared_ptr<vector<bool>>            T_aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} ); 
      shared_ptr<vector<vector<string>>>  T_idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });   
      shared_ptr<double>                  T_dummy_data ;
      string                              T_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> T_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)> T_constraints = {  &System_Info<double>::System_Info::NotAllAct };
      
      shared_ptr<TensOp<double>> TTens = Build_TensOp("T", T_dummy_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
      T_map->emplace("T", TTens);
       
    } else if ( OpName == "L" ) {  /* ---- L Tensor ----  */

      pair<double,double>                 L_factor = make_pair(1.0,1.0);
      shared_ptr<vector<string>>          L_idxs = make_shared<vector<string>>(vector<string> {"L0", "L1", "L2", "L3"});
      shared_ptr<vector<bool>>            L_aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true }); 
      shared_ptr<vector<vector<string>>>  L_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_virt, not_virt, not_core, not_core }); 
      shared_ptr<double>                  L_dummy_data;
      string                              L_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  L_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  L_constraints = { &System_Info<double>::System_Info::NotAllAct };
      
      shared_ptr<TensOp<double>> LTens = Build_TensOp("L", L_dummy_data, L_idxs, L_aops, L_idx_ranges, L_symmfuncs, L_constraints, L_factor, L_TimeSymm, false ) ;
      T_map->emplace("L", LTens);
     
    } else if ( OpName == "X" ) {  /* ---- 2el Excitation Op  ----  */

      pair<double,double>                 X_factor = make_pair(1.0,1.0);
      shared_ptr<vector<string>>          X_idxs = make_shared<vector<string>>( vector<string> {"X0", "X1", "X2", "X3"} );
      shared_ptr<vector<bool>>            X_aops = make_shared<vector<bool>>( vector<bool> { false, false, true, true } ); 
      shared_ptr<vector<vector<string>>>  X_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_virt, not_virt, not_core, not_core } ); 
      shared_ptr<double>                  X_dummy_data;
      string                              X_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  X_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  X_constraints = { &System_Info<double>::System_Info::NotAllAct };
      
      shared_ptr<TensOp<double>> XTens = Build_TensOp("X", X_dummy_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
      T_map->emplace("X", XTens);
      
    } else if ( OpName == "Y" ) {  /*---- 2el Excitation Op dag  ----  */

      pair<double,double>                 Y_factor = make_pair(1.0,1.0);
      shared_ptr<vector<string>>          Y_idxs = make_shared<vector<string>>(vector<string> {"Y0", "Y1", "Y2", "Y3"});
      shared_ptr<vector<bool>>            Y_aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false }); 
      shared_ptr<vector<vector<string>>>  Y_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt }); 
      shared_ptr<double>                  Y_dummy_data;
      string                              Y_TimeSymm = "none";
      
      vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  Y_symmfuncs = set_2el_symmfuncs();
      vector<bool(*)(shared_ptr<vector<string>>)>  Y_constraints = { &System_Info<double>::System_Info::NotAllAct };
      
      shared_ptr<TensOp<double>> YTens = Build_TensOp("Y", Y_dummy_data, Y_idxs, Y_aops, Y_idx_ranges, Y_symmfuncs, Y_constraints, Y_factor, Y_TimeSymm, false ) ;
      T_map->emplace("Y", YTens);
    }

  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
