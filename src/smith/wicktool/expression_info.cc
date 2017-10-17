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

  auto BraKet_List = make_shared<std::vector<std::shared_ptr<std::vector<std::shared_ptr<TensOp<DataType>>>>>>(BraKet_names->size()); 

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
