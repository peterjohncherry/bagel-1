#include <src/prop/proptool/algebraic_manipulator/equation.h>

//////////////////////////////////////////////////////////////////////////
using namespace std;
template<typename DataType>
void Equation_Value<DataType>::generate_all_expressions() {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_Value<DataType>::generate_all_expressions() " << endl;  
  
  for ( auto& expr_info : *expression_term_map_ ) 
    if ( expression_map_->find( expr_info.first ) != expression_map_->end() ) 
      generate_expression( expr_info.first ) ;

  return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Expression<DataType>> Equation_Value<DataType>::generate_expression( string expression_name ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Value<DataType>::Build_Expression : " << expression_name << endl;
 
  shared_ptr<vector<pair<double, string>>> term_name_list = expression_term_map_->at(expression_name);
  cout << "term_name_list->at(0).second = " << term_name_list->at(0).second << endl;
  shared_ptr<vector<BraKet<DataType>>> bk_list = term_braket_map_->at( term_name_list->at(0).second ); //term_info.second);
  cout << "got " << endl;

  for ( int ii = 1; ii != term_name_list->size(); ii++ ){
    cout << "ii = " << ii << endl;
    shared_ptr<vector<BraKet<DataType>>> term_bk_list = term_braket_map_->at( term_name_list->at(ii).second ); //term_info.second);
    for ( BraKet<DataType> bk :  *term_bk_list ) 
      bk_list->push_back(bk);
 //   bk_list->insert( bk_list->end(), term_bk_list->begin(), term_bk_list->end()); // TODO Why doesn't this work? 
  }
  cout << " build the expression " << endl;
  shared_ptr<Expression<DataType>> new_expression ; 
//  string expression_name_gen =  Build_Expression( bk_list );
//  cout << "   assert( " <<  expression_name_gen << " =?= " << expression_name  << " );" << endl;
//  assert(expression_name_gen == expression_name );

  return new_expression;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string Equation_Value<DataType>::build_expression( shared_ptr<vector<BraKet<DataType>>> expr_bk_list ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Value<DataType>::Build_Expression" << endl;
  shared_ptr< vector<pair<string, DataType>> > braKet_name_list = make_shared<vector<pair< string, DataType >>>(0);


  // This is looping over states; op sparsity depends on states, should replace with term_info_map, and
  // have double loop, outer for ket state, inner for brastate
  for ( BraKet<DataType>& braket_info : *expr_bk_list ) {

    for (string op_name : braket_info.op_list_ ) { // TODO should loop  over states defined in term_info, 
      
      auto T_loc = T_map_->find(op_name);
      if( T_loc == T_map_->end() ){ 
        shared_ptr<TensOp::TensOp<DataType>> new_op;// = System_Info<DataType>::Initialize_Tensor_Op_Info( op_name );
        T_map_->emplace( op_name, new_op );
        CTP_map_->insert( new_op->CTP_map_->begin(), new_op->CTP_map_->end());
      //TODO do state specific definition; just builds new range_map, should not have any sparsity yet ...
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
   
   // if ( braket_map_->find(braket_info.name()) == braket_map_->end() ) 
      //Set_BraKet_Ops( make_shared<vector<string>>(braket_info.op_list_), braket_info.name() ) ;
     

    braKet_name_list->push_back( make_pair( braket_info.name(), braket_info.factor_ ) );
  }

  string expression_name = "";
  for ( pair<string, DataType> name_fac_pair : *braKet_name_list ) {
    if (name_fac_pair.second != 0.0 ) 
      expression_name += "+ (" + to_string(name_fac_pair.second) + ")" + name_fac_pair.first;
  }
  cout << "new_expression_name = " << expression_name << endl;
  
  cout << "List of things in CMTP_map_ " << endl;
  for ( auto& elem : *CMTP_map_ ) 
    cout << elem.first << endl;   

  shared_ptr<Expression<DataType>> new_expression = make_shared<Expression<DataType>>( expr_bk_list, states_info_, MT_map_, CTP_map_, CMTP_map_, ACompute_map_, gamma_info_map_ ); 

  expression_map_->emplace( expression_name, new_expression );

  return expression_name;

}

//////////////////////////////////////////////////////////////////////////
template class Equation_Base<double>;
template class Equation_Value<double>;
/////////////////////////////////////////////////////////////////////////
