#include <src/prop/proptool/algebraic_manipulator/equation.h>
#include <src/prop/proptool/algebraic_manipulator/op_info.h>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Equation_Base<DataType>::set_maps(  std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map,
                                         shared_ptr< map <string, shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                         shared_ptr< map <string, shared_ptr< vector<shared_ptr<CtrOp_base>>>>> ACompute_map,
                                         shared_ptr< map< string, shared_ptr< TensOp_Base >>> MT_map,
                                         shared_ptr< map< string, shared_ptr< CtrTensorPart_Base>>> CTP_map,
                                         shared_ptr< map< char, long unsigned int>> range_prime_map  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef Equation_Base_DBG 
cout << "void Equation_Value<DataType>::set_maps" << endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
// TODO This is a very strange way of doing things, but I remember it is for a nasty reason. Either fix it
//      or find the reason and commend here....
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Equation_Base<DataType>::add_expression( string expression_name ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Base<DataType>::add_Expression : " << expression_name << " hey!" <<  endl;
 
  shared_ptr<vector<pair<DataType, string>>> term_name_list = expression_term_map_->at(expression_name);

  cout << "term_name_list->at(0).second = "; cout.flush(); cout << (*term_name_list)[0].second << endl;
  shared_ptr<vector<shared_ptr<BraKet_Base>>> bk_list = term_braket_map_->at( (*term_name_list)[0].second );
  for ( int ii = 1; ii != term_name_list->size(); ii++ ){
    shared_ptr<vector<shared_ptr<BraKet_Base>>> term_bk_list = term_braket_map_->at( term_name_list->at(ii).second );
    for ( shared_ptr<BraKet_Base> bk : *term_bk_list ) 
      bk_list->push_back(bk);
  }
  
  string expression_type = add_expression_info( bk_list ) ;

  cout << "making expression of Type : " << expression_type << endl;
  if ( expression_type == "orbital_excitation_derivative" || expression_type == "orb" ) {
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
// This is looping over states; op sparsity depends on states, should replace with term_info_map, and
// have double loop, outer for ket state, inner for brastate
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string Equation_Base<DataType>::add_expression_info( shared_ptr<vector<shared_ptr<BraKet_Base>>> expr_bk_list ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Base<DataType>::add_Expression_info  (bk input version)" << endl;

  cout << "expr_bk_list->size() = " <<  expr_bk_list->size() << endl;

  //TODO get these quantities from input 
  bool spinfree = false;
  string expression_type = expr_bk_list->front()->type_;
  cout << "expression_type =  " <<  expression_type << endl; 
  

  for ( shared_ptr<BraKet_Base>& braket_info : *expr_bk_list ) {
    if ( braket_info->type_ != expression_type )
      throw logic_error( " inconsistent term types in expression defintion : " + braket_info->type_ + " != " + expression_type ) ;

    shared_ptr<vector<shared_ptr<Op_Info>>> op_info_vec = braket_info->multiop_info_->op_info_vec();

    for (std::vector<shared_ptr<Op_Info>>::const_iterator oiv_it = op_info_vec->begin(); oiv_it != op_info_vec->end(); oiv_it++ ) {  
      
      auto T_loc = MT_map_->find( (*oiv_it)->op_name_ );
      if( T_loc == MT_map_->end() ){ 
        cout << "Adding " << (*oiv_it)->op_name_ << endl;

        shared_ptr<TensOp::TensOp<DataType>> new_op = TensOp_Info_Init::Initialize_Tensor_Op_Info<DataType>( (*oiv_it)->op_name_, range_prime_map_ );
        new_op->add_state_ids( *oiv_it );
        new_op->generate_uncontracted_ctps( (*oiv_it ) );  

        CTP_map_->insert( new_op->CTP_map()->begin(), new_op->CTP_map()->end());
        MT_map_->emplace( (*oiv_it)->op_name_, new_op );

      } else { 
        T_loc->second->add_state_ids( (*oiv_it ) );
        T_loc->second->generate_uncontracted_ctps( *oiv_it );  
      }

    }

    if( MT_map_->find( braket_info->multiop_info_->op_name_ ) == MT_map_->end() ){
      cout << "could not find MT " << braket_info->multiop_name_ << " in map" <<endl;
      vector<shared_ptr<TensOp_Base>> tensop_list(op_info_vec->size());
      vector<shared_ptr<TensOp_Base>>::iterator tl_it = tensop_list.begin();
      
      for ( vector<shared_ptr<Op_Info>>::const_iterator oiv_it = op_info_vec->begin(); oiv_it != op_info_vec->end();  oiv_it++, tl_it++ )
        *tl_it = MT_map_->at((*oiv_it)->op_name_); 

      shared_ptr<MultiTensOp::MultiTensOp<DataType>> multiop = make_shared<MultiTensOp::MultiTensOp<DataType>>( braket_info->multiop_info_->op_name_, spinfree, tensop_list, range_prime_map_ );
      MT_map_->emplace(braket_info->multiop_info_->op_name_, multiop );
    } 
    
  }
  
  return expression_type;
}
//////////////////////////////////////////////////////////////////////////
template class Equation_Base<double>;
template class Equation_Base<std::complex<double>>;
template class Equation_Value<double>;
template class Equation_Value<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
