#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include  <src/smith/wicktool/expression_info.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression_Info<DataType>::Expression_Info::Expression_Info( int nact,  int nele,  bool spinfree ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map          = make_shared< map <string, shared_ptr<TensOp<DataType>>>>();
  BraKet_map     = make_shared<  std::map <std::string, 
                                           std::shared_ptr<std::vector<std::shared_ptr< TensOp<DataType>>>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
  expression_map = make_shared< map <string, shared_ptr<Equation<DataType>>>>();
 
  nact_ = nact;
  nele_ = nele;
  nelea_ = nele;
  neleb_ = 0;
  spinfree_ = spinfree;

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression_Info<DataType>::Expression_Info::Expression_Info( int nact, int nelea, int neleb,  bool spinfree ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map          = make_shared< map <string, shared_ptr<TensOp<DataType>>>>();
  BraKet_map     = make_shared<  std::map <std::string, 
                                          std::shared_ptr<std::vector<std::shared_ptr< TensOp<DataType>>>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
  expression_map = make_shared< map <string, shared_ptr<Equation<DataType>>>>();
 
  nact_ = nact;
  nelea_ = nelea;
  neleb_ = neleb;
  nele_ = nelea + neleb;

  spinfree_ = spinfree;

  return;
}


/////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Info<DataType>::Expression_Info::Construct_Tensor_Ops() { 
/////////////////////////////////////////////////////////////////////////////////
 
  //spinfree orbital ranges
  vector<string> free     = {"cor", "act", "vir"};
  vector<string> not_core = {"act", "vir"};
  vector<string> not_act  = {"cor", "vir"};
  vector<string> not_virt = {"cor", "act"};
  vector<string> core = {"cor"};
  vector<string> act  = {"act"};
  vector<string> virt = {"vir"};
 
  /* ---- H Tensor ----  */
  pair<double,double>                H_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>         H_idxs = make_shared<vector<string>>(vector<string> {"H0", "H1", "H2", "H3"});
  shared_ptr<vector<bool>>           H_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
  shared_ptr<vector<vector<string>>> H_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
  shared_ptr<double>                 H_dummy_data;
  string                             H_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> H_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  H_constraints = { &always_true };

  auto HTens = Build_TensOp("H", H_dummy_data, H_idxs, H_aops, H_idx_ranges, H_symmfuncs, H_constraints, H_factor, H_TimeSymm, false ) ;
  T_map->emplace("H", HTens);

  
  /* ---- T Tensor ----  */
  pair<double,double>                 T_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          T_idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
  shared_ptr<vector<bool>>            T_aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} ); 
  shared_ptr<vector<vector<string>>>  T_idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });   
  shared_ptr<double>                  T_dummy_data ;
  string                              T_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> T_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)> T_constraints = { &NotAllAct };

  auto TTens = Build_TensOp("T", T_dummy_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
  T_map->emplace("T", TTens);


   /* ---- L Tensor ----  */
  pair<double,double>                 L_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          L_idxs = make_shared<vector<string>>(vector<string> {"L0", "L1", "L2", "L3"});
  shared_ptr<vector<bool>>            L_aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true }); 
  shared_ptr<vector<vector<string>>>  L_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> {  not_core, not_core, not_virt, not_virt }); 
  shared_ptr<double>                  L_dummy_data;
  string                              L_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  L_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  L_constraints = { &NotAllAct };

  shared_ptr<TensOp<double>> LTens = Build_TensOp("L", L_dummy_data, L_idxs, L_aops, L_idx_ranges, L_symmfuncs, L_constraints, L_factor, L_TimeSymm, false ) ;
  T_map->emplace("L", LTens);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<TensOp<DataType>>
Expression_Info<DataType>::Expression_Info::Build_TensOp( string op_name,
                                                          shared_ptr<DataType> tensor_data, 
                                                          shared_ptr<vector<string>> op_idxs,
                                                          shared_ptr<vector<bool>> op_aops, 
                                                          shared_ptr<vector<vector<string>>> op_idx_ranges,
                                                          vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                                                          vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                                                          pair<double,double> factor, 
                                                          string Tsymmetry,
                                                          bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Expression_Info<DataType>::Expression_Info::Build_TensOp" <<   endl;

  //NOTE: change to use proper factor
  int tmpfac = 1;
  shared_ptr<TensOp<DataType>>  New_Op = make_shared<TensOp<DataType>>(op_name, Symmetry_Funcs, Constraint_Funcs);

  New_Op->data = tensor_data;
  New_Op->initialize(*op_idxs, *op_idx_ranges, *op_aops, tmpfac, Tsymmetry);
  New_Op->get_ctrs_tens_ranges();

  return New_Op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Info<DataType>::Expression_Info::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string term_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Info::Expression_Info::Build_BraKet(shared_ptr<vector<string>> BraKet_names, string expression_name ) " << endl;
  
  shared_ptr<vector<shared_ptr<TensOp<double>>>> BraKet_Ops = make_shared<vector< shared_ptr<TensOp<double>> > >();

  for ( string name : *Op_names ) 
    BraKet_Ops->push_back(T_map->at(name));

  BraKet_map->emplace(term_name, BraKet_Ops);

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Info<DataType>::Expression_Info::Build_Expression(shared_ptr<vector<string>> BraKet_names, string expression_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Info::Expression_Info::Build_Expression(shared_ptr<vector<string>> BraKet_names, string expression_name ) " << endl;

  auto BraKet_List = make_shared<std::vector<std::shared_ptr<std::vector<std::shared_ptr<TensOp<DataType>>>>>>(); 

  for ( int ii = 0 ; ii != BraKet_names->size() ; ii++ ) 
    BraKet_List->at(ii) = BraKet_map->at(BraKet_names->at(ii)) ;

  shared_ptr<Equation<DataType>> new_expression = make_shared<Equation<DataType>>(BraKet_List);

  expression_map->emplace(expression_name, new_expression);

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
