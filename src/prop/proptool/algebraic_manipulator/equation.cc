#include <src/prop/proptool/algebraic_manipulator/equation.h>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Equation_Base<DataType>::set_maps(  std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map,
                                         shared_ptr< map <string, shared_ptr< GammaInfo >>> gamma_info_map,
                                         shared_ptr< map <string, shared_ptr< vector<shared_ptr<CtrOp_base>>>>> ACompute_map,
                                         shared_ptr< map< string, shared_ptr< TensOp_Base >>> MT_map,
                                         shared_ptr< map< string, shared_ptr< CtrTensorPart_Base>>> CTP_map,
                                         shared_ptr< map< char, long unsigned int>> range_prime_map  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void Equation_Value<DataType>::set_maps" << endl;

  expression_map_ = expression_map;   if(!expression_map_ ) throw logic_error( "expression_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " ); 
  gamma_info_map_ = gamma_info_map;   if(!gamma_info_map_ ) throw logic_error( "gamma_info_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " );
  ACompute_map_   = ACompute_map;     if(!ACompute_map_   ) throw logic_error( "ACompute_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " );
  MT_map_         = MT_map;           if(!MT_map_         ) throw logic_error( "MT_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " );
  CTP_map_        = CTP_map;          if(!CTP_map_        ) throw logic_error( "CTP_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " );
  range_prime_map_  = range_prime_map; if(!range_prime_map_        ) throw logic_error( "range_prime_map_ is null when set in Equation_Base<DataType>::set_maps !! Aborting!! " );

  cout << "leaving void Equation_Value<DataType>::set_maps" << endl;
  return;
}

//////////////////////////////////////////////////////////////////////////
using namespace std;
template<typename DataType>
void Equation_Base<DataType>::generate_all_expressions() {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_Base<DataType>::generate_all_expressions() " << endl;  
  
  for ( auto& expr_info : *expression_term_map_ ){ 
    cout <<"expr_info.first = " << expr_info.first << endl;
    if ( expression_map_->find( expr_info.first ) == expression_map_->end() )
      add_expression(expr_info.first);
  }

  return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Equation_Base<DataType>::add_expression( string expression_name ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Base<DataType>::Build_Expression : " << expression_name << endl;
 
  shared_ptr<vector<pair<double, string>>> term_name_list = expression_term_map_->at(expression_name);
  cout << "term_name_list->at(0).second = " << term_name_list->at(0).second << endl;
  shared_ptr<vector<BraKet<DataType>>> bk_list = term_braket_map_->at( term_name_list->at(0).second ); //term_info.second);

  for ( int ii = 1; ii != term_name_list->size(); ii++ ){
    shared_ptr<vector<BraKet<DataType>>> term_bk_list = term_braket_map_->at( term_name_list->at(ii).second ); //term_info.second);
    for ( BraKet<DataType> bk :  *term_bk_list ) 
      bk_list->push_back(bk);
  }
  
  string expression_type = add_expression_info( bk_list ) ;

  cout << "making expression" << endl;
  if ( expression_type == "orbital_excitation_derivative"  ) {
    cout << "expression_type = " << expression_type <<  endl;
    shared_ptr<Expression_Orb_Exc_Deriv<DataType>>  new_exp = make_shared<Expression_Orb_Exc_Deriv<DataType>>( bk_list, states_info_, MT_map_, CTP_map_, ACompute_map_, gamma_info_map_, expression_type );
    new_exp->generate_algebraic_task_list();
    expression_map_->emplace( expression_name, new_exp);
  } else  if ( expression_type == "full"  ) {
    cout << "expression_type = " << expression_type <<  endl;
    shared_ptr<Expression_Full<DataType>>  new_exp = make_shared<Expression_Full<DataType>>( bk_list, states_info_, MT_map_, CTP_map_, ACompute_map_, gamma_info_map_, expression_type );
    new_exp->generate_algebraic_task_list();
    expression_map_->emplace( expression_name, new_exp );
  } else { 
    throw std::logic_error( "have not implemented expression type... Aborting!!" );  
  //  assert(false);
  } 
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string Equation_Base<DataType>::add_expression_info( shared_ptr<vector<BraKet<DataType>>> expr_bk_list ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Base<DataType>::Build_Expression bk input" << endl;
  shared_ptr< vector<pair<string, DataType>> > braKet_name_list = make_shared<vector<pair< string, DataType >>>(0);

  string expression_type = "full"; 
  // This is looping over states; op sparsity depends on states, should replace with term_info_map, and
  // have double loop, outer for ket state, inner for brastate
  for ( BraKet<DataType>& braket_info : *expr_bk_list ) {

    for (string op_name : braket_info.op_list_ ) { // TODO should loop  over states defined in term_info, 

      if (op_name == "X" ) 
        expression_type = "orbital_excitation_derivative";

      cout << MT_map_->size();
      auto T_loc = MT_map_->find(op_name);
      if( T_loc == MT_map_->end() ){ 
        shared_ptr<TensOp::TensOp<DataType>> new_op = TensOp_Info_Init::Initialize_Tensor_Op_Info<DataType>( op_name, range_prime_map_ );
        CTP_map_->insert( new_op->CTP_map()->begin(), new_op->CTP_map()->end());
        MT_map_->emplace( op_name, new_op );
      }
    }

    if( MT_map_->find( braket_info.multiop_name_ ) == MT_map_->end() ){
      cout << "could not find MT " << braket_info.multiop_name_ << " in map" <<endl;
      vector<shared_ptr<TensOp_Base>> SubOps(braket_info.op_list_.size());

      for (int ii = 0 ; ii != braket_info.op_list_.size() ; ii++ )  
        SubOps[ii] = MT_map_->at(braket_info.op_list_[ii]); 


      shared_ptr<MultiTensOp::MultiTensOp<DataType>> multiop = make_shared<MultiTensOp::MultiTensOp<DataType>>( braket_info.multiop_name_, /*spinfree_ = */ true, SubOps, range_prime_map_ );
      multiop->get_ctrs_tens_ranges();
      CTP_map_->insert( multiop->CTP_map()->begin(), multiop->CTP_map()->end());
      MT_map_->emplace(braket_info.multiop_name_, multiop );
    } 
    
    braKet_name_list->push_back( make_pair( braket_info.name(), braket_info.factor_ ) );
    cout << "Pushed " <<  braket_info.multiop_name_ << " back into braket_name_list" << endl;
  }
  
  cout << "Full CMTP map" << endl;
  for ( auto& elem  : *CTP_map_) { 
    cout <<  elem.first << endl;
  } 
  return expression_type;
}

//////////////////////////////////////////////////////////////////////////
template class Equation_Base<double>;
template class Equation_Value<double>;
/////////////////////////////////////////////////////////////////////////
