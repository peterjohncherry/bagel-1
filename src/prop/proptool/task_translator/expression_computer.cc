#include <bagel_config.h>
#include <src/prop/proptool/task_translator/expression_computer.h>
#include <src/prop/proptool/debugging_utils.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/proputils.h>

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

  scalar_results_map   = make_shared<map<string, DataType >>(); //TODO dumb, check why not in header and fix  
  tensor_results_map_  = make_shared<map<string, shared_ptr<Tensor_<DataType>>>>(); 
  gamma_contrib_maps_ = make_shared<map<string, shared_ptr<map<string, DataType> >>>(); 
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

  tensop_machine_ = make_shared<TensOp_Computer::TensOp_Computer<DataType>>( expression->ACompute_map_, expression->CTP_map_, range_conversion_map_,
                                                                             tensop_data_map_, moint_computer_ );
  tensop_machine_->get_tensor_data_blocks( expression->required_blocks_ );

  // define seperate generic name for each of the required blocks in the target tensor.
  // Copy the results from the contractions into these places in the data map.
  // Perhaps have seperate field in expression which is map of these blocks; can be removed when we a finished, and helps seperate updates from actual tensor
  for ( auto& exc_block_G_to_A_map_elem : *(expression->exc_block_G_to_A_map() ) ) {
 
    string target_block_name = exc_block_G_to_A_map_elem.first; 
    shared_ptr<Tensor_<DataType>> target_block_data;
    {
    auto target_block_data_loc = tensor_results_map_->find( target_block_name ); 
    if ( target_block_data_loc == tensor_results_map_->end() ) {
      target_block_data = tensop_machine_->get_tensor( *(expression->CTP_map_->at(target_block_name)->full_id_ranges() ), true, 0.0 );
      tensor_results_map_->emplace( target_block_name, target_block_data );
    } else {
      target_block_data = target_block_data_loc->second;
    }
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
            post_a_gamma_contraction_data = tensop_machine_->get_tensor(  *(a_intermediate_info->post_gamma_contraction_ranges()), true, 0.0 );
          } else {
            throw logic_error (" Have fully contracted term; implies all active indexes in excitation operator! Aborting!!! " );
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
                          
              shared_ptr<AContribInfo_Base> a_contrib = a_contrib_elem.second;    

              if ( first_a_contrib ) {  
                first_a_contrib = false;
                shared_ptr<vector<string>> buff = a_contrib->a_block_ranges();
                vector<string> pagc_ranges = get_subvector( *(a_contrib->a_block_ranges()), a_contrib->id_orders().front() );
                tensop_machine_->build_tensor( "pre_a_gamma_contraction_data" , pagc_ranges, 0.0 );
              }
 
              if (!check_acontrib_factors(*a_contrib))
                continue;

              tensop_machine_->Calculate_CTP( *a_contrib );

              string  a_contrib_name = a_contrib_elem.first;    
              tensop_machine_->sum_different_orderings(  "pre_a_gamma_contraction_data", a_contrib_name, a_contrib->factors(), a_contrib->id_orders() ); 
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
    cout << elem.first << "->norm() = " << setprecision(13) << elem.second->norm() << endl; 
 
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

  tensop_machine_->get_tensor_data_blocks( expression->required_blocks_ );
  tensop_machine_->gamma_data_map_ = gamma_computer_->gamma_data_map(); 

  DataType result = 0.0;
  shared_ptr<map< string, DataType >> g_result_map = make_shared<map< string, DataType >>();

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(expression->gamma_info_map_) ) {

    string gamma_name = AG_contrib.first; 
    
    tensop_machine_->build_tensor( "a_combined_data" , *(expression->gamma_info_map_->at(gamma_name)->id_ranges()) ); 
    shared_ptr<Tensor_<DataType>> A_combined_data = tensop_data_map_->at("a_combined_data");

    // Loop through A-tensors needed for this gamma
    auto A_contrib_loc = expression->G_to_A_map()->find( gamma_name );
    if ( (A_contrib_loc != expression->G_to_A_map()->end()) &&  (A_contrib_loc->second->size() != 0) ) {
      for ( auto& A_contrib_map_elem : *A_contrib_loc->second ) {
      
        string  a_contrib_name = A_contrib_map_elem.first;    
        shared_ptr<AContribInfo_Base> a_contrib = A_contrib_map_elem.second;    

        if (check_acontrib_factors(*a_contrib)){
          tensop_machine_->Calculate_CTP( *a_contrib );
          tensop_machine_->add_acontrib_to_target( "a_combined_data", a_contrib );
        } else { 
          continue;
        }
      }

      if ( gamma_name != "ID" ) {
        gamma_computer_->get_gamma( gamma_name );
        DataType tmp_result = tensop_machine_->dot_arg1_with_gamma( "a_combined_data", gamma_name );

        auto grm_loc = g_result_map->find( gamma_name);
        if ( grm_loc == g_result_map->end() ) {
          g_result_map->emplace(gamma_name, tmp_result);
        } else { 
          grm_loc->second += tmp_result;
        }

        cout << "result +=  tmp_result : " << result << " += " << tmp_result ; cout.flush(); 
        result += tmp_result;
        cout << " --> result = " << result << endl;

      } else {
        DataType tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::sum_tensor_elems( A_combined_data ) ;
        auto grm_loc = g_result_map->find( gamma_name);
        if ( grm_loc == g_result_map->end() ) {
          g_result_map->emplace(gamma_name, tmp_result);
        } else { 
          grm_loc->second += tmp_result;
        }

        cout << "result +=  tmp_result : " << result << " += " << tmp_result ; cout.flush(); 
        result += tmp_result ; 
        cout << " --> result = " << result << endl;

      }
    }
  }

  if (new_result ) {
    scalar_results_map->emplace( expression_name, result );
    gamma_contrib_maps_->emplace( expression_name, g_result_map );
  } else {
    scalar_results_map->at( expression_name ) = result;
    gamma_contrib_maps_->at( expression_name ) = g_result_map;
  }
  print_scalar_result( expression_name, true );  
  return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void
Expression_Computer::Expression_Computer<DataType>::print_scalar_result( std::string expression_name,
                                                                         bool print_gamma_contribs    ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::print_scalar_result" << endl;  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr< map<string, DataType >> g_result_map = gamma_contrib_maps_->at(expression_name);
  cout << endl << "Contributions from different gamma terms " << endl;
  for ( auto elem : *g_result_map ) 
    cout << elem.first << " :  " << setprecision(13) << elem.second << endl; 

  cout << "=========================================================================================================" << endl;
  cout << "                             RESULT FOR " << expression_name << " = " << scalar_results_map->at(expression_name) <<  endl;
  cout << "=========================================================================================================" << endl;
  cout << endl; 
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::check_rdms() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::check_rdms" << endl;  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  auto bob =  gamma_computer_->gamma_data_map();
  shared_ptr<Tensor_<DataType>> rdm1;
  shared_ptr<Tensor_<DataType>> rdm2;
  auto v2_act =  tensop_data_map_->at("H_{00}_aaaa");
  cout << "v2_act->norm() = "<< v2_act->norm() << endl;
  bool got_rdm1 = false;
  bool got_rdm2 = false;

  for ( auto elem : *bob ) { 
    if ( elem.second->rank() == 2 ) {      
      rdm1 = elem.second;
      tensop_data_map_->emplace("rdm1",  rdm1 );
      got_rdm1 = true;
    } else if ( elem.second->rank() == 4 ) { 
      rdm2 = make_shared<Tensor_<DataType>>(*(elem.second));
      tensop_data_map_->emplace("rdm2",  rdm2 );
      got_rdm2 = true;
    }
  }

  if ( got_rdm1 && got_rdm2 ) {

    vector<int> id_pos = { 1, 2 };
    vector<int> reorder_vector = { 2, 3, 0, 1 };
    tensop_data_map_->emplace( "rdm1_test", make_shared<Tensor_<DataType>>(*(tensop_data_map_->at("rdm1"))) );
    tensop_data_map_->emplace( "rdm2_test", make_shared<Tensor_<DataType>>(*(tensop_data_map_->at("rdm2"))) );
    shared_ptr<Tensor_<DataType>> rdm2_klij =  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( tensop_data_map_->at("rdm2_test" ), reorder_vector  );

    tensop_data_map_->emplace( "rdm2_klij", rdm2_klij);

    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::add_tensor_along_trace( tensop_data_map_->at("rdm2_klij"), tensop_data_map_->at("rdm1_test"), id_pos, -1.0 );

    vector<string> free4 = { "free", "free", "free", "free" };
    tensop_data_map_->emplace( "v2_smith_order", moint_computer_->calculate_v2_smith( free4 ) );
    vector<string> act4 = { "a", "a", "a", "a" };
    tensop_machine_->get_sub_tensor( "v2_smith_order", "v2_smith_order_act",  act4 );

    cout << "tensop_data_map_->at(\"v2_smith_order_act\")->norm() = " << tensop_data_map_->at("v2_smith_order_act")->norm()  << endl;

    cout << " rdm2_klij->dot_product( tensop_data_map_->at(\"v2_smith_order_act\") ) = "; cout.flush();
    cout << rdm2_klij->dot_product( tensop_data_map_->at("v2_smith_order_act") ) << endl;

  }

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
bool Expression_Computer::Expression_Computer<DataType>::check_acontrib_factors(AContribInfo_Base& AC_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER
cout << "Expression_Computer::Expression_Computer<DataType>::check_acontrib_factors" << endl;
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
  return !skip;
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void Expression_Computer::Expression_Computer<DataType>::test_sum_reordered_tensor_list(){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_EXPRESSION_COMPUTER 
cout << "Expression_Computer::Expression_Computer<DataType>::test_sum_reordered_tensor_list()" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<size_t> test_ranges_64 = { 6, 4 };
  vector<size_t> test_ranges_46 = { 4, 6 };
  vector<size_t> test_ranges_444 = { 4, 4, 4 };
  vector<size_t> test_ranges_446 = { 4, 4, 6 };
  vector<size_t> test_ranges_644 = { 6, 4, 4 };
  vector<size_t> test_ranges_345 = { 3, 4, 5 };
  vector<size_t> test_ranges_543 = { 5, 4, 3 };
  vector<size_t> test_ranges_4444 = { 4, 4, 4, 4 };
  vector<size_t> test_ranges_4422 = { 4, 4, 2, 2 };
  vector<size_t> test_ranges_2244 = { 2, 2, 4, 4 };
  vector<size_t> test_ranges_3355 = { 3, 3, 5, 5 }; 
  vector<size_t> test_ranges_5533 = { 5, 5, 3, 3 }; 
  
  vector<size_t> test_ranges_6325 = { 6, 3, 2, 5 }; 
  vector<size_t> test_ranges_2536 = { 2, 5, 3, 6 }; 

  vector<size_t> max_blocks_46 = { 4, 6 };
  vector<size_t> max_blocks_22 = { 2, 2 };
  vector<size_t> max_blocks_24 = { 2, 4 };
  vector<size_t> max_blocks_42 = { 4, 2 };
  vector<size_t> max_blocks_222 = { 2, 2, 2 }; 
  vector<size_t> max_blocks_224 = { 2, 2, 4 }; 
  vector<size_t> max_blocks_422 = { 4, 2, 2 }; 
  vector<size_t> max_blocks_4444 = { 4, 4, 4, 4 }; 
  vector<size_t> max_blocks_4422 = { 4, 4, 2, 2 }; 
  vector<size_t> max_blocks_2244 = { 2, 2, 4, 4 }; 
  vector<size_t> max_blocks_3322 = { 3, 3, 2, 2 }; 
  vector<size_t> max_blocks_2233 = { 2, 2, 3, 3 }; 
  
  vector<size_t> max_blocks_1524 = { 1, 5, 2, 4 }; 
  vector<size_t> max_blocks_2451 = { 2, 4, 5, 1 }; 

  vector<int> ord_01 = { 0, 1 };
  vector<int> ord_10 = { 1, 0 };
  vector<int> ord_012 = { 0, 1, 2 };
  vector<int> ord_210 = { 2, 1, 0 };
  vector<int> ord_0123 = { 0, 1, 2, 3 };
  vector<int> ord_3210 = { 3, 2, 1, 0 };
  vector<int> ord_3201 = { 3, 2, 0, 1 };
  vector<int> ord_2301 = { 2, 3, 0, 1 };
  vector<int> ord_2310 = { 2, 3, 1, 0 };

  vector<size_t> target_ranges  = test_ranges_2536;
  vector<size_t> summand_ranges = test_ranges_6325;

  vector<size_t> target_max_blocks  = max_blocks_2451;
  vector<size_t> summand_max_blocks = max_blocks_1524;

  vector<vector<int>> summand_reorderings = { ord_2310 };
  vector<double> summand_factors = { -1.0 };

  tensop_machine_->build_test_tensor( "target", target_ranges, target_max_blocks );
  shared_ptr<Tensor_<double>> target = tensop_data_map_->at("target");
  Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( target, 0.0);
  print_tensor_with_indexes( target, "target" ); cout << endl << endl; 

  tensop_machine_->build_test_tensor( "summand", summand_ranges, summand_max_blocks );
  shared_ptr<Tensor_<double>> summand = tensop_data_map_->at("summand");
  print_tensor_with_indexes( summand, "summand" ); cout << endl << endl; 

  Tensor_Arithmetic::Tensor_Arithmetic<double>::add_list_of_reordered_tensors( target, summand, summand_reorderings, summand_factors );
  print_tensor_with_indexes( target, "post target" ); cout << endl << endl; 
 
  {

  vector<int> summand_reordering_inverse = get_ascending_order(summand_reorderings[0]);
  shared_ptr<Tensor_<double>> summand_reordered = tensop_machine_->reorder_block_Tensor( "summand", summand_reorderings[0] );
  print_tensor_with_indexes( summand_reordered, "summand_reordered" ); cout << endl << endl; 

  summand_reordered->ax_plus_y ( 1.0 , target ); 
  print_tensor_with_indexes( summand_reordered, "summand-post target" ); cout << endl << endl; 
  } 
  throw logic_error("dumb"); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
