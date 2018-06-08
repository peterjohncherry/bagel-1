#include <bagel_config.h>
#include <src/prop/proptool/task_translator/expression_computer.h>
#include <src/prop/proptool/debugging_utils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;
using namespace Debugging_Utils;

#define __DEBUG_EXPRESSION_COMPUTER
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
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout <<  "Expression_Computer::Expression_Computer::evaluate_expression : name input : " << expression_name <<  endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Expression<DataType>> expr = expression_map_->at(expression_name);
  evaluate_expression(expr);
  return;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression( shared_ptr<Expression<DataType>> expression ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::Expression_Computer<DataType>::evaluate_expression" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
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
#ifdef __DEBUG_EXPRESSION_COMPUTER 
cout <<  "Expression_Computer::Expression_Computer::evaluate_expression_orb_exc_deriv : " << expression->name() << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  tensop_machine_ = make_shared<TensOp_Computer::TensOp_Computer<DataType>>( expression->ACompute_map_, expression->CTP_map_, range_conversion_map_, tensop_data_map_,
                                                                                 moint_computer_ );
 
  tensop_machine_->get_tensor_data_blocks( expression->required_blocks_ );

  // define seperate generic name for each of the required blocks in the target tensor.
  // Copy the results from the contractions into these places in the data map.
  // Perhaps have seperate field in expression which is map of these blocks; can be removed when we a finished, and helps seperate updates from actual tensor
  for ( auto& exc_block_G_to_A_map_elem : *(expression->exc_block_G_to_A_map() ) ) {
 
    string target_block_name = exc_block_G_to_A_map_elem.first; 

    shared_ptr<Tensor_<DataType>> target_block_data;
    auto target_block_data_loc = tensor_results_map_->find( target_block_name ); 
    if ( target_block_data_loc == tensor_results_map_->end() ) {
      shared_ptr<vector<string>> target_block_ranges = expression->CTP_map_->at(target_block_name)->full_id_ranges();
      target_block_data = make_shared<Tensor_<DataType>>( *(tensop_machine_->Get_Bagel_IndexRanges( target_block_ranges )));
      target_block_data->allocate();
      target_block_data->zero(); 
      tensor_results_map_->emplace( target_block_name, target_block_data );
    } else {
      target_block_data = target_block_data_loc->second;
    }

    auto G_to_A_map = exc_block_G_to_A_map_elem.second; 
       
    for ( auto AG_contrib : *(expression->gamma_info_map_) ) {
    
      string gamma_name = AG_contrib.first; 
       
      auto final_reorderings_map_loc = G_to_A_map->find( gamma_name );
      if ( (final_reorderings_map_loc != G_to_A_map->end()) &&  (final_reorderings_map_loc->second->size() != 0) ) {

        // loop over different final reorderings 
        for ( auto& final_reordering_elem : *(final_reorderings_map_loc->second) ) {     // Loop through different final reorderings.

          shared_ptr<AContribInfo_OrbExcDeriv<DataType>> a_intermediate_info = dynamic_pointer_cast<AContribInfo_OrbExcDeriv<DataType>>(final_reordering_elem.second);

          auto gamma_contraction_pos_map = final_reordering_elem.second->gamma_pos_map();

          shared_ptr<Tensor_<DataType>> post_a_gamma_contraction_data;
          if ( a_intermediate_info->post_gamma_contraction_ranges()->size() > 0 ) {
            post_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( *(tensop_machine_->Get_Bagel_IndexRanges( a_intermediate_info->post_gamma_contraction_ranges() )));
            post_a_gamma_contraction_data->allocate();
            post_a_gamma_contraction_data->zero(); 
          } else {
            throw logic_error (" Expression_Computer::evaluate_expression_orb_exc_deriv : " + expression->name() + " : Should never have fully contracted term; implies all active indexes in excitation operator! Aborting! " );
          }  

          // loop over different gamma contraction positions         
          for ( auto a_contrib_map_elem : *gamma_contraction_pos_map ){ 

            vector<int> gamma_contraction_pos = a_contrib_map_elem.first;
            vector<int> a_contraction_pos( gamma_contraction_pos.size() );
            iota( a_contraction_pos.begin(), a_contraction_pos.end(), 0 );

            pair<vector<int>, vector<int>> gamma_a_contractions = make_pair( gamma_contraction_pos, a_contraction_pos ); 

            shared_ptr<Tensor_<DataType>> pre_a_gamma_contraction_data;
                 
            //TODO get rid of this hack to find A_block_ranges 
            bool first_a_contrib = true; 
            for ( auto& a_contrib_elem : *(a_contrib_map_elem.second) ) {
                          
              string  a_contrib_name = a_contrib_elem.first;    

              shared_ptr<AContribInfo_Base> a_contrib = a_contrib_elem.second;    

              if ( first_a_contrib ) {  
                first_a_contrib = false;
                shared_ptr<vector<string>> buff = a_contrib->a_block_ranges();
                vector<int> aid_order = a_contrib->id_orders().front();
                vector<int>::iterator aio_it = aid_order.begin();
                shared_ptr<vector<string>> pagc_ranges = make_shared<vector<string>>(buff->size()) ;
                for ( vector<string>::iterator pagcr_it = pagc_ranges->begin(); pagcr_it != pagc_ranges->end() ; pagcr_it++, aio_it++ )
                  *pagcr_it = buff->at( *aio_it ) ;
               
                pre_a_gamma_contraction_data = make_shared<Tensor_<DataType>>( *(tensop_machine_->Get_Bagel_IndexRanges( pagc_ranges ) ) );
                pre_a_gamma_contraction_data->allocate();
                pre_a_gamma_contraction_data->zero(); 
              }
 
              if (check_AContrib_factors(*a_contrib))
                continue;

              tensop_machine_->Calculate_CTP( *a_contrib );

              for ( int qq = 0 ; qq != a_contrib->id_orders().size(); qq++){
                shared_ptr<Tensor_<DataType>> a_contrib_reordered = tensop_machine_->reorder_block_Tensor( a_contrib_name, a_contrib->id_order(qq));
                pre_a_gamma_contraction_data->ax_plus_y( (DataType)(a_contrib->factor(qq).first), a_contrib_reordered );
              }

              tensop_data_map_->emplace( "pre_a_gamma_contraction_data", pre_a_gamma_contraction_data );
            }

            if ( gamma_name != "ID" ) {

              gamma_computer_->get_gamma( gamma_name );
              shared_ptr<Tensor_<DataType>> gamma_data = gamma_computer_->gamma_data(gamma_name);

              shared_ptr<Tensor_<DataType>> tmp_result =   Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_different_tensors( gamma_data, pre_a_gamma_contraction_data,  gamma_a_contractions );
 
              assert( post_a_gamma_contraction_data->size_alloc() == tmp_result->size_alloc() );

              post_a_gamma_contraction_data->ax_plus_y( (DataType)(1.0), tmp_result );

            } else {
             assert ( post_a_gamma_contraction_data->size_alloc() == pre_a_gamma_contraction_data->size_alloc() ) ;
             post_a_gamma_contraction_data->ax_plus_y( (DataType)(1.0), pre_a_gamma_contraction_data ) ;
            }

          }

          if ( a_intermediate_info->post_contraction_reordering()->size() == target_block_data->rank() ) {
            shared_ptr<Tensor_<DataType>> reordered_tensor_block =  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( post_a_gamma_contraction_data, *(a_intermediate_info->post_contraction_reordering()) );
            assert (target_block_data->size_alloc() == reordered_tensor_block->size_alloc() ) ;

            target_block_data->ax_plus_y((DataType)(1.0), reordered_tensor_block );
            tensor_results_map_->at( target_block_name ) =  target_block_data ; 
            
          } else if ( a_intermediate_info->post_contraction_reordering()->size() != 0 ) { 

            Tensor_Arithmetic::Tensor_Arithmetic<DataType>::add_tensor_along_trace( target_block_data, post_a_gamma_contraction_data, *(a_intermediate_info->target_block_positions()), 1.0 );
            tensop_data_map_->at( target_block_name ) =  target_block_data ; 
            tensor_results_map_->at( target_block_name ) =  target_block_data ; 

          } else {
            throw logic_error ("Expression_computer:: should never end up with contribution to target term having rank 1 " ); 
          }
        }
      }
    }
  }                  

  cout << endl << "Contributions from different gamma terms " << endl;
  for ( auto elem : *tensor_results_map_ )
    cout << elem.first << "->norm() = " << elem.second->norm() << endl; 
  
 
  return;  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::evaluate_expression_full( shared_ptr<Expression<DataType>> expression ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout <<  "Expression_Computer::Expression_Computer::evaluate_expression_full : " << expression->name() << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  string expression_name = expression->name();


  bool new_result = ( scalar_results_map->find( expression_name ) == scalar_results_map->end() ); 
  if ( !new_result )  
    cout << "WARNING : You have already calculated this expression....." << expression_name << " = " << scalar_results_map->at( expression_name ) << endl;

  tensop_machine_ = make_shared<TensOp_Computer::TensOp_Computer<DataType>>( expression->ACompute_map_, expression->CTP_map_, range_conversion_map_, tensop_data_map_,
                                                                             moint_computer_ );
 
  test_trace_subtraction();

  tensop_machine_->get_tensor_data_blocks( expression->required_blocks_ );
  throw logic_error("die here for testing");
  DataType result = 0.0;
  map< string, DataType > g_result_map;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(expression->gamma_info_map_) ) {

    string gamma_name = AG_contrib.first; 
    
    tensop_machine_->build_tensor( "a_combined_data" , *(expression->gamma_info_map_->at(gamma_name)->id_ranges()) ); 
    shared_ptr<Tensor_<DataType>> A_combined_data = tensop_data_map_->at("a_combined_data");

    // Loop through A-tensors needed for this gamma
    auto A_contrib_loc = expression->G_to_A_map()->find( gamma_name );
    if ( (A_contrib_loc != expression->G_to_A_map()->end()) &&  (A_contrib_loc->second->size() != 0) ) {
      for ( auto& A_contrib_map_elem : *A_contrib_loc->second ) {
      
        string  A_contrib_name = A_contrib_map_elem.first;    
        shared_ptr<AContribInfo_Base> A_contrib = A_contrib_map_elem.second;    

        if (check_AContrib_factors(*A_contrib))
          continue;
      
	print_AContraction_list( expression->ACompute_map_->at(A_contrib_name), A_contrib_name );
        tensop_machine_->Calculate_CTP( *A_contrib );
           
        tensop_machine_->sum_different_orderings( "a_combined_data" , A_contrib_name, A_contrib->factors(), A_contrib->id_orders() );
  
        cout << "added " << A_contrib_name << endl; 
        cout << "=========================================================================================================" << endl << endl;
      }


      if ( gamma_name != "ID" ) {

        gamma_computer_->get_gamma( gamma_name );
        cout << "gamma_computer_->gamma_data_map()->at("+gamma_name+")->norm() = "; cout.flush(); cout << gamma_computer_->gamma_data_map()->at(gamma_name)->norm() << endl;
        cout << " result = " << result << endl;
        DataType tmp_result = A_combined_data->dot_product( gamma_computer_->gamma_data_map()->at(gamma_name) );
        cout << " tmp_result = " << tmp_result << endl;
        g_result_map.emplace(gamma_name, tmp_result) ;
        result += tmp_result;

      } else {

        DataType tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::sum_tensor_elems( A_combined_data ) ;
        g_result_map.emplace(gamma_name, tmp_result) ;
        result += tmp_result ; 

      }
    }
  }

  cout << endl << "Contributions from different gamma terms " << endl;
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

  { // TEST 
    auto bob =  gamma_computer_->gamma_data_map();
    shared_ptr<Tensor_<DataType>> rdm1;
    shared_ptr<Tensor_<DataType>> rdm2;
    auto v2_act =  tensop_data_map_->at("H_{00}_aaaa");
    for ( auto elem : *bob ) { 
      if ( elem.second->rank() == 2 ) {      
        cout << elem.first << "->norm() = " ; cout.flush(); cout << elem.second->norm() << endl;
        rdm1 = elem.second;
        tensop_data_map_->emplace("rdm1",  rdm1 );
      } else if ( elem.second->rank() == 4 ) { 
        cout << elem.first << "->norm() = " ; cout.flush(); cout << elem.second->norm() << endl;
        rdm2 = make_shared<Tensor_<DataType>>(*(elem.second));
        tensop_data_map_->emplace("rdm2",  rdm2 );
      }
    }
    vector<int> id_pos = {1, 2};
    cout << "rdm1->size_alloc() = " ; cout.flush() ; cout << rdm1->size_alloc() <<endl;
    cout << "rdm2->size_alloc() = " ; cout.flush() ; cout << rdm2->size_alloc() <<endl;
    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::add_tensor_along_trace( rdm2, rdm1, id_pos, -1.0 ); 
    cout << "rdm1->norm() = "; cout.flush(); cout << rdm1->norm() << endl;
    cout << "rdm2->norm() = "; cout.flush(); cout << rdm2->norm() << endl;
  } // END TEST 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::print_AContraction_list(shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list, string A_contrib_name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::print_AContraction_list" << endl;  
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

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
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::Expression_Computer<DataType>::check_AContrib_factors" << endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////
  bool  skip = false;
  for ( int qq = 0 ; qq != AC_info.id_orders().size(); qq++) { 
     if ( AC_info.factor(qq).first != 0 || AC_info.factor(qq).second !=0) { 
       
       print_vector( AC_info.id_order(qq), AC_info.name() +".id_orders["+to_string(qq)+"]") ;
       cout << "has non-zero factor (" <<  AC_info.factor(qq).first <<","<< AC_info.factor(qq).second << ")" <<endl; 
       break;
     }
     if ( qq == AC_info.id_orders().size()-1 ) {
       print_vector( AC_info.id_order(qq), AC_info.name() + ".id_orders["+to_string(qq)+"]") ;
       cout << "has NO non-zero factors" <<endl; 
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
#ifdef __DEBUG_EXPRESSION_COMPUTER 
cout << "Expression_Computer::Expression_Computer<DataType>::set_gamma_maps" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  gamma_computer_->set_maps(range_conversion_map_, expression_map_->at(expression_name)->gamma_info_map_, gamma_data_map, sigma_data_map, civec_data_map );
  return;
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::test_trace_subtraction(){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER 
cout << "Expression_Computer::Expression_Computer<DataType>::test_trace_substraction()" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////


  vector<size_t> test_range2 = { 6, 7 };
  vector<size_t> max_blocks2 = { 4, 4 }; 
  tensop_machine_->build_test_tensor( "test_tens2", test_range2, max_blocks2 );
  shared_ptr<Tensor_<double>> test_tens2 = tensop_data_map_->at("test_tens2");
  print_tensor_with_indexes( test_tens2, "test_tens2" ); cout << endl << endl; 

  vector<size_t> test_range4 = { 4, 6, 4, 7 };
  vector<size_t> max_blocks4 = { 2, 4, 2, 4 }; 
  tensop_machine_->build_test_tensor( "test_tens4_new", test_range4, max_blocks4 );
  shared_ptr<Tensor_<double>> test_tens4 = tensop_data_map_->at("test_tens4_new");
  print_tensor_with_indexes( test_tens4, "test_tens4_new" ); cout << endl << endl; 
  Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( test_tens4, 0.0);
   
  vector<int> id_pos = { 1, 3 };
  Tensor_Arithmetic::Tensor_Arithmetic<double>::add_tensor_along_trace( test_tens4, test_tens2, id_pos, -1.0 ); 
  print_tensor_with_indexes( test_tens4, "test_tens4_minus {1,3}"  );

  {
  vector<size_t> test_range4746 = { 4, 7, 4, 6 };
  vector<size_t> max_blocks2424 = { 2, 4, 2, 4 }; 
  tensop_machine_->build_test_tensor( "test_tens_4746", test_range4746, max_blocks2424 );
  shared_ptr<Tensor_<double>> test_tens4746 = tensop_data_map_->at("test_tens_4746");
  print_tensor_with_indexes( test_tens4746, "test_tens_4746" ); cout << endl << endl; 
  Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( test_tens4746, 0.0);
  
  vector<int> id_pos2 = { 3, 1 };
  Tensor_Arithmetic::Tensor_Arithmetic<double>::add_tensor_along_trace( test_tens4746, test_tens2, id_pos2, -1.0 ); 
  print_tensor_with_indexes( test_tens4746, "test_tens4_minus {3,1}"  );
  }
  return;
}
//Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( test_tens4, 0.0);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
////  {
  //
  //    tensop_machine_->build_test_tensor( "tester", test_range4, max_blocks4 );
  //    shared_ptr<Tensor_<double>> tester = tensop_data_map_->at("tester");
  //    print_tensor_with_indexes( tester, "tester pre reordering" ); cout << endl << endl;
  //    vector<int> reordering = { 0, 2, 1, 3 };
  //
  //    shared_ptr<Tensor_<double>> tester_reordered = tensop_machine_->reorder_block_Tensor( "tester", reordering );
  //    tensop_data_map_->emplace( "tester_reordered", tester_reordered );
  //    print_tensor_with_indexes( tester_reordered, "tester_reordered { 0, 2, 1, 3 } " ); cout << endl << endl;
  //
  //    shared_ptr<Tensor_<double>> tester_orig = tensop_machine_->reorder_block_Tensor( "tester_reordered", reordering );
  //    
  //    tensop_data_map_->emplace( "tester_orig", tester_orig );
  //    print_tensor_with_indexes( tester_orig, "tester_orig" ); cout << endl << endl;
  //
  //    tensop_machine_->build_test_tensor( "tester2", test_range2, max_blocks2 );
  //    shared_ptr<Tensor_<double>> tester2 = tensop_data_map_->at("tester2");
  //    print_tensor_with_indexes( tester2, "tester2 pre reordering" ); cout << endl << endl;
  //    vector<int> reordering2 = { 1, 0 };
  //
  //    shared_ptr<Tensor_<double>> tester2_reordered = tensop_machine_->reorder_block_Tensor( "tester2", reordering2 );
  //    tensop_data_map_->emplace( "tester2_reordered", tester2_reordered );
  //    print_tensor_with_indexes( tester2_reordered, "tester2_reordered {1, 0}" ); cout << endl << endl;
  //
  //  }
