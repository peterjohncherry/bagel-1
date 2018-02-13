#include <src/prop/proptool/algebraic_manipulator/equation.h>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Equation_Base<DataType>::set_maps(  std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map,
                                         shared_ptr< map <string, shared_ptr< GammaInfo >>> gamma_info_map,
                                         shared_ptr< map <string, shared_ptr< vector<shared_ptr<CtrOp_base>>>>> ACompute_map,
                                         shared_ptr< map< string, shared_ptr< TensOp::TensOp<DataType>>>> T_map,
                                         shared_ptr< map< string, shared_ptr< MultiTensOp::MultiTensOp<DataType>>>> MT_map,
                                         shared_ptr< map< string, shared_ptr< CtrTensorPart<DataType>>>> CTP_map,     
                                         shared_ptr< map< string, shared_ptr< CtrMultiTensorPart<DataType> >>> CMTP_map){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void Equation_Value<DataType>::set_maps" << endl;

  expression_map_ = expression_map;
  gamma_info_map_ = gamma_info_map;
  ACompute_map_   = ACompute_map;
  T_map_          = T_map;
  MT_map_         = MT_map;
  CTP_map_        = CTP_map;   
  CMTP_map_       = CMTP_map;

  cout << "leaving void Equation_Value<DataType>::set_maps" << endl;
  return;
}

//////////////////////////////////////////////////////////////////////////
using namespace std;
template<typename DataType>
void Equation_Base<DataType>::generate_all_expressions() {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_Value<DataType>::generate_all_expressions() " << endl;  
  
  for ( auto& expr_info : *expression_term_map_ ){ 
    cout <<"expr_info.first = " << expr_info.first << endl;
    if ( expression_map_->find( expr_info.first ) == expression_map_->end() )
      expression_map_->emplace( expr_info.first, build_expression( expr_info.first ) );
  }

  return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Expression<DataType>> Equation_Base<DataType>::build_expression( string expression_name ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Value<DataType>::Build_Expression : " << expression_name << endl;
 
  shared_ptr<vector<pair<double, string>>> term_name_list = expression_term_map_->at(expression_name);
  cout << "term_name_list->at(0).second = " << term_name_list->at(0).second << endl;
  shared_ptr<vector<BraKet<DataType>>> bk_list = term_braket_map_->at( term_name_list->at(0).second ); //term_info.second);

  for ( int ii = 1; ii != term_name_list->size(); ii++ ){
    shared_ptr<vector<BraKet<DataType>>> term_bk_list = term_braket_map_->at( term_name_list->at(ii).second ); //term_info.second);
    for ( BraKet<DataType> bk :  *term_bk_list ) 
      bk_list->push_back(bk);
 //   bk_list->insert( bk_list->end(), term_bk_list->begin(), term_bk_list->end()); // TODO Why doesn't this work? 
  }

  return build_expression( bk_list );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Expression<DataType>> Equation_Base<DataType>::build_expression( shared_ptr<vector<BraKet<DataType>>> expr_bk_list ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Value<DataType>::Build_Expression" << endl;
  shared_ptr< vector<pair<string, DataType>> > braKet_name_list = make_shared<vector<pair< string, DataType >>>(0);

  // This is looping over states; op sparsity depends on states, should replace with term_info_map, and
  // have double loop, outer for ket state, inner for brastate
  for ( BraKet<DataType>& braket_info : *expr_bk_list ) {

    for (string op_name : braket_info.op_list_ ) { // TODO should loop  over states defined in term_info, 
      
      auto T_loc = T_map_->find(op_name);
      if( T_loc == T_map_->end() ){ 
        shared_ptr<TensOp::TensOp<DataType>> new_op = TensOp_Info_Init::Initialize_Tensor_Op_Info<DataType>( op_name );
        CTP_map_->insert( new_op->CTP_map_->begin(), new_op->CTP_map_->end());
        T_map_->emplace( op_name, new_op );
      }
    } 

    if( MT_map_->find( braket_info.multiop_name_ ) == MT_map_->end() ){

      vector<shared_ptr<TensOp::TensOp<DataType>>> SubOps(braket_info.op_list_.size());
      for (int ii = 0 ; ii != braket_info.op_list_.size() ; ii++ )  
        SubOps[ii] = T_map_->at(braket_info.op_list_[ii]); 

      shared_ptr<MultiTensOp::MultiTensOp<DataType>> multiop = make_shared<MultiTensOp::MultiTensOp<DataType>>( braket_info.multiop_name_, /*spinfree_ = */ true, SubOps );
      multiop->get_ctrs_tens_ranges();
      CTP_map_->insert( multiop->CTP_map()->begin(), multiop->CTP_map()->end());
      CMTP_map_->insert( multiop->CMTP_map()->begin(), multiop->CMTP_map()->end());
      MT_map_->emplace(braket_info.multiop_name_, multiop );
    } 

    braKet_name_list->push_back( make_pair( braket_info.name(), braket_info.factor_ ) );
  }
  
  return  make_shared<Expression<DataType>>( expr_bk_list, states_info_, MT_map_, CTP_map_, CMTP_map_, ACompute_map_, gamma_info_map_ );

}

//////////////////////////////////////////////////////////////////////////
template class Equation_Base<double>;
template class Equation_Value<double>;
/////////////////////////////////////////////////////////////////////////
