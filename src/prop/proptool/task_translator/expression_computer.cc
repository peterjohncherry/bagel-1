#include <bagel_config.h>
#include <src/prop/proptool/task_translator/expression_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;

///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
Expression_Computer::Expression_Computer<DataType>::Expression_Computer( shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>>   gamma_computer,
                                                                         shared_ptr<map< string, shared_ptr<Expression<DataType>>>> expression_map,
                                                                         shared_ptr<map< string, shared_ptr<IndexRange>>>           range_conversion_map,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<DataType>>>>    tensop_data_map,
                                                                         shared_ptr<MOInt_Computer<DataType>> moint_computer ) : 
  gamma_computer_(gamma_computer), expression_map_(expression_map), range_conversion_map_(range_conversion_map), tensop_data_map_(tensop_data_map),
  moint_computer_(moint_computer) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
  scalar_results_map = make_shared<map< string, DataType >>(); //TODO dumb, check why not in header and fix  
  tensor_results_map_ = make_shared<map< string, shared_ptr<Tensor_<DataType>>>>(); 
}  
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression( string expression_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout <<  "Expression_Computer::Expression_Computer::evaluate_expression : name input : " << expression_name <<  endl;

  shared_ptr<Expression<DataType>> expr = expression_map_->at(expression_name);
  cout << "got " << expression_name << " info "<< endl;
  evaluate_expression(expr);
  return;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression( shared_ptr<Expression<DataType>> expression ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Computer::Expression_Computer<DataType>::evaluate_expression" << endl; 
  
  if ( expression->type_ == "full" ) {
    evaluate_expression_full( expression );
 
  } else if ( expression->type_ == "orb") {
    evaluate_expression_orb_exc_deriv( expression );
  
  } else { 
    cout << "expression type \"" << expression->type_ << "\" is not implemented" <<  endl; 
  }  
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression_orb_exc_deriv( shared_ptr<Expression<DataType>> expression ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout <<  "Expression_Computer::Expression_Computer::evaluate_expression_orb_exc_deriv : " << expression->name() << endl;

  auto TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer<DataType>>( expression->ACompute_map_, expression->CTP_map_, range_conversion_map_, tensop_data_map_,
                                                                                 moint_computer_ );
  cout << endl<< "-------------------- Contents of CTP_map ----------------" << endl;
  for ( auto& elem : *(expression->CTP_map_) ) 
    cout << elem.first << endl;

  cout << endl << endl;
 
  TensOp_Machine->get_tensor_data_blocks( expression->required_blocks_ );

  // define seperate generic name for each of the required blocks in the target tensor.
  // Copy the results from the contractions into these places in the data map.
  // Perhaps have seperate field in expression which is map of these blocks; can be removed when we a finished, and helps seperate updates from actual tensor

  for ( auto& exc_block_G_to_A_map_elem : *(expression->exc_block_G_to_A_map() ) ) {
 
    string target_block_name = exc_block_G_to_A_map_elem.first; 

    shared_ptr<Tensor_<DataType>> target_block_data;
    auto target_block_data_loc = tensor_results_map_->find( target_block_name ); 
    if ( target_block_data_loc == tensor_results_map_->end() ) {
      // MUST ALLOCATE
    } else {
      target_block_data = target_block_data_loc->second;
    }

    auto G_to_A_map = exc_block_G_to_A_map_elem.second; 
       
    for ( auto AG_contrib : *(expression->gamma_info_map_) ) {
    
      string gamma_name = AG_contrib.first; 
      cout << "gamma_name = " << gamma_name << endl;
       
      auto final_reorderings_map_loc = G_to_A_map->find( gamma_name );
      if ( (final_reorderings_map_loc != G_to_A_map->end()) &&  (final_reorderings_map_loc->second->size() != 0) ) {

        // loop over different final reorderings 
        for ( auto& final_reordering_elem : *(final_reorderings_map_loc->second) ) {     // Loop through different final reorderings.

          shared_ptr<AContribInfo_OrbExcDeriv<DataType>> a_intermediate_info = dynamic_pointer_cast<AContribInfo_OrbExcDeriv<DataType>>(final_reordering_elem.second);

          cout << endl << endl; 
          cout << " final_reordering_elem.second->name_ = " << final_reordering_elem.second->name_  << endl;

          auto gamma_contraction_pos_map = final_reordering_elem.second->gamma_pos_map();

          shared_ptr<Tensor_<DataType>> post_a_gamma_contraction_data;
          if ( gamma_name != "ID" ) {
            post_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( *(TensOp_Machine->Get_Bagel_IndexRanges( a_intermediate_info->post_gamma_contraction_ranges() )));
          } else {
            post_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
          }  
          post_a_gamma_contraction_data->allocate();
          post_a_gamma_contraction_data->zero(); 

          // loop over different gamma contraction positions         
          for ( auto a_contrib_map_elem : *(gamma_contraction_pos_map) ){ 

            vector<int> gamma_contraction_pos = a_contrib_map_elem.first;
            vector<int> a_contraction_pos(gamma_contraction_pos.size());
            iota( a_contraction_pos.begin(), a_contraction_pos.end(), 0 );

            pair<vector<int>, vector<int>> gamma_a_contractions = make_pair( gamma_contraction_pos, a_contraction_pos ); 
            print_vector( gamma_a_contractions.first, "gamma_contraction_pos" ); cout.flush(); print_vector( gamma_a_contractions.second , "     a_contraction_pos" ); cout << endl; 

            shared_ptr<Tensor_<DataType>> pre_a_gamma_contraction_data;
            if ( gamma_name != "ID" ) {
              pre_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( *(TensOp_Machine->Get_Bagel_IndexRanges(expression->gamma_info_map_->at(gamma_name)->id_ranges())) );
            } else {
              pre_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
            }  

            pre_a_gamma_contraction_data->allocate();
            pre_a_gamma_contraction_data->zero(); 
                 
            // loop over different a_contractoins
            for ( auto& a_contrib_elem : *(a_contrib_map_elem.second) ) {
                            
              string  a_contrib_name = a_contrib_elem.first;    
              cout << "a_contrib_name = " << a_contrib_name ; cout.flush(); 

              shared_ptr<AContribInfo_Base> a_contrib = a_contrib_elem.second;    

              if (check_AContrib_factors(*a_contrib))
                continue;
//              print_AContraction_list( expression->ACompute_map_->at(a_contrib_name), a_contrib_name );
              TensOp_Machine->Calculate_CTP( *a_contrib );
              
              for ( int qq = 0 ; qq != a_contrib->id_orders().size(); qq++){

                print_vector( a_contrib->id_order(qq) , "   id_reordering ") ; cout << endl;
 
                shared_ptr<Tensor_<DataType>> a_contrib_reordered = TensOp_Machine->reorder_block_Tensor( a_contrib_name, make_shared<vector<int>>(a_contrib->id_order(qq)) );
                pre_a_gamma_contraction_data->ax_plus_y( (DataType)(a_contrib->factor(qq).first), a_contrib_reordered );
              
              }
            }
            if ( gamma_name != "ID" ) {
      
              gamma_computer_->get_gamma( gamma_name );
              shared_ptr<Tensor_<DataType>> gamma_data = gamma_computer_->gamma_data(gamma_name); 
              shared_ptr<Tensor_<DataType>> tmp_result =   Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_different_tensors( gamma_data, pre_a_gamma_contraction_data, gamma_a_contractions  );
              post_a_gamma_contraction_data->ax_plus_y( (DataType)(1.0), tmp_result ) ;
      
            } else {
              cout << "I ignore identity based terms for now " << endl; 
            }

          } 
          print_vector( *(a_intermediate_info->post_contraction_reordering()) , "post_contraction_reordering" ); cout << endl;                
          shared_ptr<Tensor_<DataType>> reordered_tensor_block =  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( post_a_gamma_contraction_data, a_intermediate_info->post_contraction_reordering() );
            
          cout << "X" << endl;
          target_block_data->ax_plus_y((DataType)(1.0), reordered_tensor_block );
        }
      }
    }
  }                  
  cout << "leaving Expression_Computer::evaluate_expression_orb_exc_deriv" << endl;
  return;  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression_full( shared_ptr<Expression<DataType>> expression ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Computer::Expression_Computer::evaluate_expression_full : " << expression->name() << endl;

  string expression_name = expression->name();
  cout << endl << endl;
  cout << " --------- required_blocks ---------" << endl; 
  for ( std::shared_ptr<Range_Block_Info> block : *(expression->required_blocks_) ) 
    cout << block->name() << endl;
  cout << endl << endl;

  bool new_result = ( scalar_results_map->find( expression_name ) == scalar_results_map->end() ); 
  if ( !new_result )  
    cout << "WARNING : You have already calculated this expression....." << expression_name << " = " << scalar_results_map->at( expression_name ) << endl;

  auto TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer<DataType>>( expression->ACompute_map_, expression->CTP_map_, range_conversion_map_, tensop_data_map_,
                                                                                 moint_computer_ );
 
  TensOp_Machine->get_tensor_data_blocks( expression->required_blocks_ );

  DataType result = 0.0;
  map< string, DataType > g_result_map;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(expression->gamma_info_map_) ) {

    string gamma_name = AG_contrib.first; 
    cout << " gamma_name  = " << gamma_name << endl; 

    shared_ptr<Tensor_<DataType>> A_combined_data;
    // Build A_tensor to hold sums of different A-tensors.
    if ( gamma_name != "ID" ) {
      A_combined_data = make_shared<Tensor_<DataType>>( *(TensOp_Machine->Get_Bagel_IndexRanges(expression->gamma_info_map_->at(gamma_name)->id_ranges())) );
    } else {
      A_combined_data = make_shared<Tensor_<DataType>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
    }  
    A_combined_data->allocate();
    A_combined_data->zero(); 
 
    // Loop through A-tensors needed for this gamma
    auto A_contrib_loc = expression->G_to_A_map()->find( gamma_name );
    if ( (A_contrib_loc != expression->G_to_A_map()->end()) &&  (A_contrib_loc->second->size() != 0) ) {
      for ( auto&  A_contrib_map_elem : *A_contrib_loc->second ) {
      
        string  A_contrib_name = A_contrib_map_elem.first;    
        shared_ptr<AContribInfo_Base> A_contrib = A_contrib_map_elem.second;    

        if (check_AContrib_factors(*A_contrib))
          continue;
      
	print_AContraction_list( expression->ACompute_map_->at(A_contrib_name), A_contrib_name );
        TensOp_Machine->Calculate_CTP( *A_contrib );

        if ( gamma_name != "ID" ) {
          if ( tensop_data_map_->find(A_contrib_name) != tensop_data_map_->end() ) {

            for ( int qq = 0 ; qq != A_contrib->id_orders().size(); qq++){
              shared_ptr<Tensor_<DataType>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib->id_order(qq)) );
              A_combined_data->ax_plus_y( (DataType)(A_contrib->factor(qq).first), A_contrib_reordered );

	      cout << " A_contrib->factor(" << qq<<").first), tensop_data_map_->at(" << A_contrib_name << ")-norm() = ";
	      cout <<  A_contrib->factor(qq).first << ", " <<  tensop_data_map_->at(A_contrib_name)->norm() << endl;
              cout << "A_combined_data->norm() = "<<  A_combined_data->norm() << endl;
            }

        } else {

           cout << A_contrib_name << " not found in map; must be formed from direct product" << endl; //TODO is a way to avoid this, implement it
           for ( int qq = 0 ; qq != A_contrib->id_orders().size(); qq++){

              shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> CTP_vec = expression->CTP_map_->at(A_contrib_name)->CTP_vec() ;
              vector<string> sub_tensor_names(CTP_vec->size()); 

              for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ )  
                sub_tensor_names[rr] = CTP_vec->at(rr)->name();
              
              cout << "sub_tensor_names = [ " ; cout.flush(); 
              for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ ) { 
                cout << sub_tensor_names[rr]  << " " ; cout.flush(); 
              }  cout << "]" << endl;
              
              shared_ptr<Tensor_<DataType>> A_contrib_data = TensOp_Machine->direct_product_tensors( sub_tensor_names );//TODO fix so uses piecewise contraction where possible 
              tensop_data_map_->emplace( A_contrib_name, A_contrib_data );
              shared_ptr<Tensor_<DataType>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib->id_order(qq)) );
              A_combined_data->ax_plus_y( (DataType)(A_contrib->factor(qq).first), A_contrib_reordered ); //TODO replace with interface function in tensop_computer

              cout << " A_contrib->factor(" << qq<<").first), tensop_data_map_->at(" << A_contrib_name << ")-norm() = ";
	      cout <<  A_contrib->factor(qq).first << ", " <<  tensop_data_map_->at(A_contrib_name)->norm() << endl;
              cout << "A_combined_data->norm() = "<<  A_combined_data->norm() << endl;
            }
          }

        } else {
  
          if ( tensop_data_map_->find(A_contrib_name) == tensop_data_map_->end() ) {cout << A_contrib_name <<  " not yet in map, must form from direct product" << endl;

            shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> CTP_vec = expression->CTP_map_->at(A_contrib_name)->CTP_vec() ;
            vector<string> sub_tensor_names(CTP_vec->size()); 
            for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ )
              sub_tensor_names[rr] = CTP_vec->at(rr)->name();
            
            shared_ptr<Tensor_<DataType>> A_contrib_data = TensOp_Machine->direct_product_tensors( sub_tensor_names );
            for ( int qq = 0 ; qq != A_contrib->id_orders().size(); qq++){
              A_combined_data->ax_plus_y( (DataType)(A_contrib->factor(qq).first), A_contrib_data );

              cout << " A_contrib->factor(" << qq<<").first), tensop_data_map_->at(" << A_contrib_name << ")-norm() = ";
              cout <<  A_contrib->factor(qq).first << ", " <<  tensop_data_map_->at(A_contrib_name)->norm() << endl;
	      cout << "A_combined_data->norm() = "<<  A_combined_data->norm() << endl;

             }
          } else {

            for ( int qq = 0 ; qq != A_contrib->id_orders().size(); qq++){
              cout << " A_contrib->factor(" << qq<<").first), tensop_data_map_->at(" << A_contrib_name << ")->norm() = ";
              cout <<  A_contrib->factor(qq).first << ", " <<  tensop_data_map_->at(A_contrib_name)->norm() << endl;
              A_combined_data->ax_plus_y( (DataType)(A_contrib->factor(qq).first), tensop_data_map_->at(A_contrib_name) );
	      cout << "A_combined_data->norm() = "<<  A_combined_data->norm() << endl;
            }
          }
        }
  
        cout << "added " << A_contrib_name << endl; 
        cout << "=========================================================================================================" << endl << endl;
      }
     
      if ( gamma_name != "ID" ) {

        cout << "A13" << endl;
        gamma_computer_->get_gamma( gamma_name );
        cout << "gamma_computer_->gamma_data_map()->at(gamma_name)->norm() = " <<  gamma_computer_->gamma_data_map()->at(gamma_name)->norm() << endl;
        cout << "A_combined_data->norm() = "<<  A_combined_data->norm() << endl;
 
        DataType tmp_result = A_combined_data->dot_product( gamma_computer_->gamma_data_map()->at(gamma_name) );
        g_result_map.emplace(gamma_name, tmp_result) ;
        result += tmp_result;

        cout << "A14" << endl;
      } else {

        cout << "A15" << endl;
        DataType tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::sum_tensor_elems( A_combined_data ) ; cout << "tmp_result = " << tmp_result << endl;
        g_result_map.emplace(gamma_name, tmp_result) ;
        result += tmp_result ; 

        cout << "A16" << endl;
      }
    }
  }
  cout << endl << endl;

  cout << "Contributions from different gamma terms " << endl;
  for ( auto elem : g_result_map ) 
    cout << elem.first << " :  " << elem.second << endl; 

  cout << "=========================================================================================================" << endl;
  cout << "                             RESULT FOR " << expression_name << " = " << result <<  endl;
  cout << "=========================================================================================================" << endl;
  cout << endl; 

  if (new_result ) {
    scalar_results_map->emplace( expression_name, result );
  } else {
    scalar_results_map->at( expression_name ) = result;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::print_AContraction_list(shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list, string A_contrib_name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Expression_Computer::print_AContraction_list" << endl;  
  cout << "=========================================================================================================" << endl;
  cout << A_contrib_name << endl;
  cout << "=========================================================================================================" << endl;
  for (shared_ptr<CtrOp_base> ctr_op : *ACompute_list ){
    if ( ctr_op->ctr_type()[0] == 'd' ){
      cout << "[" << ctr_op->T1name() << " , " << ctr_op->T2name() << " , (";
      cout << ctr_op->T1_ctr_abs_pos() << "," <<  ctr_op->T2_ctr_abs_pos() << ")" << " , " << ctr_op->Tout_name() << " ] " ; cout << ctr_op->ctr_type() << endl;
  
    } else if (ctr_op->ctr_type()[0] == 's' ){
      cout << "[" << ctr_op->T1name() << " , " << ctr_op->T1name() << " , (";
      cout << ctr_op->ctr_abs_pos().first << "," <<  ctr_op->ctr_abs_pos().second << ")" << " , " << ctr_op->Tout_name() << " ] " ;   cout << ctr_op->ctr_type()  << endl;
    } else { 
      cout << ctr_op->ctr_type() << endl;
    }
    
  }
  cout << "=========================================================================================================" << endl;
 
  return;
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
bool Expression_Computer::Expression_Computer<DataType>::check_AContrib_factors(AContribInfo_Base& AC_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  bool  skip = false;
  for ( int qq = 0 ; qq != AC_info.id_orders().size(); qq++) { 
     if ( AC_info.factor(qq).first != 0 || AC_info.factor(qq).second !=0) { 
       print_vector( AC_info.id_order(qq), "AC_info.id_orders["+to_string(qq)+"]") ; cout << "has non-zero factor (" <<  AC_info.factor(qq).first <<","<< AC_info.factor(qq).second << ")" <<endl; 
       break;
     }
     if ( qq == AC_info.id_orders().size()-1 ) {
       print_vector( AC_info.id_order(qq), "AC_info.id_orders["+to_string(qq)+"]") ; cout << "has NO non-zero factors" <<endl; 
       skip =true;
     }
  } 
  return skip;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::set_gamma_maps( string expression_name, 
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<DataType>>>>    gamma_data_map,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<DataType>>>>    sigma_data_map,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<DataType>>>>    civec_data_map  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Computer::Expression_Computer<DataType>::set_gamma_maps" << endl;

  gamma_computer_->set_maps(range_conversion_map_, expression_map_->at(expression_name)->gamma_info_map_, gamma_data_map, sigma_data_map, civec_data_map );
  return;
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
