#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/expression_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;

///////////////////////////////////////////////////////////////////////////////////////////////////////
template < class DataType >
Expression_Computer::Expression_Computer<DataType>::Expression_Computer( shared_ptr<const Dvec >                                  civectors_in, 
                                                                         shared_ptr<map< string, shared_ptr<Expression<double>>>> Expression_map_in,
                                                                         shared_ptr<map< string, shared_ptr<IndexRange>>>         range_conversion_map_in,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<double>>>>    TensOp_data_map_in,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<double>>>>    Gamma_data_map_in,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<double>>>>    Sigma_data_map_in,
                                                                         shared_ptr<map< string, shared_ptr<Tensor_<double>>>>    CIvec_data_map_in  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Computer::Expression_Computer::Expression_Computer " << endl;

  Expression_map = Expression_map_in;
  range_conversion_map =  range_conversion_map_in;
  TensOp_data_map      =  TensOp_data_map_in; 
  Gamma_data_map       =  Gamma_data_map_in; 
  Sigma_data_map       =  Sigma_data_map_in; 
  CIvec_data_map       =  CIvec_data_map_in; 
  civectors            =  civectors_in;                   

  scalar_results_map = make_shared<map< string, double >>(); 

}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template < class DataType >
void Expression_Computer::Expression_Computer<DataType>::Evaluate_Expression( string Expression_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Computer::Expression_Computer::Evaluate_Expression : " << Expression_name <<  endl;

 bool new_result = ( scalar_results_map->find( Expression_name ) == scalar_results_map->end() ); 
 if ( !new_result )  
   cout << "WARNING : You have already calculated this expression....." << Expression_name << " = " << scalar_results_map->at( Expression_name ) << endl;

  shared_ptr<Expression<double>> Expr = Expression_map->at(Expression_name); cout << "got " << Expression_name << " info "<< endl;

  shared_ptr<TensOp_Computer::TensOp_Computer> TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer>( Expr->ACompute_map, Expr->CTP_map, range_conversion_map, TensOp_data_map);
  cout << "Built_TensOp_Computer" << endl;


  B_Gamma_Computer::B_Gamma_Computer B_Gamma_Machine( civectors, range_conversion_map, Expr->GammaMap, Gamma_data_map, Sigma_data_map, CIvec_data_map );
  cout << "Built_B_Gamma_Computer" << endl;

  double result = 0.0;
  map< string, double > g_result_map;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(Expr->GammaMap) ) {

    string Gamma_name = AG_contrib.first;
    cout << " Gamma_name  = " << Gamma_name << endl; 

    print_vector( *(Expr->GammaMap->at(Gamma_name)->id_ranges), "id_ranges" );

    shared_ptr<Tensor_<double>> A_combined_data;
    // Build A_tensor to hold sums of different A-tensors
    if ( Gamma_name != "ID" ) {
      A_combined_data = make_shared<Tensor_<double>>( *(TensOp_Machine->Get_Bagel_IndexRanges(Expr->GammaMap->at(Gamma_name)->id_ranges)) );
    } else {
      A_combined_data = make_shared<Tensor_<double>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
    }  
    A_combined_data->allocate();
    A_combined_data->zero(); 
 
    cout << "Built an A-tensor to hold contributions" << endl;

    // Loop through A-tensors needed for this gamma
    auto  A_contrib_loc =  Expr->G_to_A_map->find(Gamma_name);
  
    if ( A_contrib_loc !=  Expr->G_to_A_map->end() ) {
   
      for ( auto A_contrib : *A_contrib_loc->second ) {
    
        if (check_AContrib_factors(A_contrib.second))
          continue;

	print_AContraction_list(Expr->ACompute_map->at(A_contrib.first), A_contrib.first);
        TensOp_Machine->Calculate_CTP(A_contrib.first);

        if ( Gamma_name != "ID" ) {
          for ( int qq = 0 ; qq != A_contrib.second.id_orders.size(); qq++){
            if ( TensOp_data_map->find(A_contrib.first) != TensOp_data_map->end() ) {
              cout << endl; Print_Tensor(TensOp_data_map->at(A_contrib.first), A_contrib.first); cout << endl << endl << endl;
            } else {
              cout << A_contrib.first << " is not done" << endl; 
              //TODO  fix this so it recognises when tensor cannot be found due to decomposition.
              //      Needs alteration here, probably in ctrtensop, and new tensor contraction routine (diff tensor taking vector, not pair, would do) 
            }
            shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib.first, make_shared<vector<int>>(A_contrib.second.id_order(qq)) );
            A_combined_data->ax_plus_y( (double)(A_contrib.second.factor(qq).first), A_contrib_reordered );
          }

        } else {

          for ( int qq = 0 ; qq != A_contrib.second.id_orders.size(); qq++)
            A_combined_data->ax_plus_y( (double)(A_contrib.second.factor(qq).first), TensOp_data_map->at(A_contrib.first) );

        }
  
        cout << "added " << A_contrib.first << endl; 
        cout << "=========================================================================================================" << endl << endl;
      }
      
      if ( Gamma_name != "ID" ) {

        B_Gamma_Machine.get_gamma( Gamma_name );
       
        double tmp_result = A_combined_data->dot_product( Gamma_data_map->at(Gamma_name) );
        cout << "A_combined_data->dot_product( Gamma_data_map->at(" << Gamma_name << ") ) = " << tmp_result  << endl;
        g_result_map.emplace(Gamma_name, tmp_result) ;
        result += tmp_result;

      } else {

        Print_Tensor( A_combined_data, " A_combined_data for 1D " ) ; cout << endl;

        double tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<double>::sum_tensor_elems( A_combined_data ) ;
        cout << "tmp_result = " << tmp_result << endl;
        g_result_map.emplace(Gamma_name, tmp_result) ;
        result += tmp_result ; 

      }
    }
  }
  cout << endl << endl;

  cout << "Contributions from different gamma terms " << endl;
  for ( auto elem : g_result_map ) 
    cout << elem.first << " :  " << elem.second << endl; 

  cout << "=========================================================================================================" << endl;
  cout << "                             RESULT FOR " << Expression_name << " = " << result <<  endl;
  cout << "=========================================================================================================" << endl;
  cout << endl; 

  if (new_result ) {
    scalar_results_map->emplace( Expression_name, result );
  } else {
    scalar_results_map->at( Expression_name ) = result;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < class DataType >
void Expression_Computer::Expression_Computer<DataType>::print_AContraction_list(shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list, string A_contrib_name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
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
template < class DataType >
bool Expression_Computer::Expression_Computer<DataType>::check_AContrib_factors(AContribInfo& AC_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  bool  skip = false;
  for ( int qq = 0 ; qq != AC_info.id_orders.size(); qq++) { 
     if ( AC_info.factor(qq).first != 0 || AC_info.factor(qq).second !=0) 
       break;
     if ( qq == AC_info.id_orders.size()-1 ) {
       skip =true;
     }
  } 
  return skip;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
 
#endif
