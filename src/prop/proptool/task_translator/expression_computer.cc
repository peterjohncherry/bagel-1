#include <bagel_config.h>
#include <src/prop/proptool/task_translator/expression_computer.h>

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

  shared_ptr<TensOp_Computer::TensOp_Computer<double>> TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer<double>>( Expr->ACompute_map_, Expr->CTP_map_, range_conversion_map, TensOp_data_map);
  B_Gamma_Computer::B_Gamma_Computer B_Gamma_Machine( civectors, range_conversion_map, Expr->GammaMap_, Gamma_data_map, Sigma_data_map, CIvec_data_map );

  double result = 0.0;
  map< string, double > g_result_map;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(Expr->GammaMap_) ) {

    string Gamma_name = AG_contrib.first;  cout << " Gamma_name  = " << Gamma_name << endl; 

    shared_ptr<Tensor_<double>> A_combined_data;
    // Build A_tensor to hold sums of different A-tensors.
    if ( Gamma_name != "ID" ) {
      A_combined_data = make_shared<Tensor_<double>>( *(TensOp_Machine->Get_Bagel_IndexRanges(Expr->GammaMap_->at(Gamma_name)->id_ranges())) );
    } else {
      A_combined_data = make_shared<Tensor_<double>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
    }  
    A_combined_data->allocate();
    A_combined_data->zero(); 
 
    // Loop through A-tensors needed for this gamma
    auto A_contrib_loc = Expr->G_to_A_map_->find( Gamma_name );
    if ( (A_contrib_loc != Expr->G_to_A_map_->end()) &&  (A_contrib_loc->second->size() != 0) ) {
  
      for ( auto  A_contrib_map_elem : *A_contrib_loc->second ) {
      
        string  A_contrib_name = A_contrib_map_elem.first;    
        AContribInfo A_contrib = A_contrib_map_elem.second;    

        if (check_AContrib_factors(A_contrib))
          continue;
      
	print_AContraction_list( Expr->ACompute_map_->at(A_contrib_name), A_contrib_name );
        TensOp_Machine->Calculate_CTP( A_contrib );

        if ( Gamma_name != "ID" ) {
       
          if ( TensOp_data_map->find(A_contrib_name) != TensOp_data_map->end() ) { cout << A_contrib_name << " found in map" << endl;

            for ( int qq = 0 ; qq != A_contrib.id_orders.size(); qq++){
              shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib.id_order(qq)) );
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_reordered );
            }

         } else { cout << A_contrib_name << " not found in map; must be formed from direct product" << endl; //TODO is a way to avoid this, implement it
            
           for ( int qq = 0 ; qq != A_contrib.id_orders.size(); qq++){

              shared_ptr<vector<shared_ptr<CtrTensorPart<double>>>> CTP_vec = Expr->CMTP_map_->at(A_contrib_name)->CTP_vec ;
              vector<string> sub_tensor_names(CTP_vec->size()); 

              for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ )
                sub_tensor_names[rr] = CTP_vec->at(rr)->myname();

              shared_ptr<Tensor_<double>> A_contrib_data = TensOp_Machine->direct_product_tensors( sub_tensor_names );//TODO fix so uses piecewise contraction where possible 
              TensOp_data_map->emplace( A_contrib_name, A_contrib_data );
              shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib.id_order(qq)) );
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_reordered );

            }
          }

        } else {
  
          if ( TensOp_data_map->find(A_contrib_name) == TensOp_data_map->end() ) {cout << A_contrib_name <<  " not yet in map, must form from direct product" << endl;

            shared_ptr<vector<shared_ptr<CtrTensorPart<double>>>> CTP_vec = Expr->CMTP_map_->at(A_contrib_name)->CTP_vec ;
            vector<string> sub_tensor_names(CTP_vec->size()); 
            for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ )
              sub_tensor_names[rr] = CTP_vec->at(rr)->myname();
            
            shared_ptr<Tensor_<double>> A_contrib_data = TensOp_Machine->direct_product_tensors( sub_tensor_names );
            for ( int qq = 0 ; qq != A_contrib.id_orders.size(); qq++)
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_data );

          } else {

            for ( int qq = 0 ; qq != A_contrib.id_orders.size(); qq++)
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), TensOp_data_map->at(A_contrib_name) );

          }

        }
  
        cout << "added " << A_contrib_name << endl; 
        cout << "=========================================================================================================" << endl << endl;
      }
     
      if ( Gamma_name != "ID" ) {

        B_Gamma_Machine.get_gamma( Gamma_name );
 
        double tmp_result = A_combined_data->dot_product( Gamma_data_map->at(Gamma_name) );
        g_result_map.emplace(Gamma_name, tmp_result) ;
        result += tmp_result;

      } else {

        double tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<double>::sum_tensor_elems( A_combined_data ) ; cout << "tmp_result = " << tmp_result << endl;
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
template < class DataType >
bool Expression_Computer::Expression_Computer<DataType>::check_AContrib_factors(AContribInfo& AC_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  bool  skip = false;
  for ( int qq = 0 ; qq != AC_info.id_orders.size(); qq++) { 
     if ( AC_info.factor(qq).first != 0 || AC_info.factor(qq).second !=0) { 
       print_vector( AC_info.id_orders[qq], "AC_info.id_orders["+to_string(qq)+"]") ; cout << "has non-zero factor (" <<  AC_info.factor(qq).first <<","<< AC_info.factor(qq).second << ")" <<endl; 
       break;
     }
     if ( qq == AC_info.id_orders.size()-1 ) {
       print_vector( AC_info.id_orders[qq], "AC_info.id_orders["+to_string(qq)+"]") ; cout << "has NO non-zero factors" <<endl; 
       skip =true;
     }
  } 
  return skip;
}


  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < class DataType >
void Expression_Computer::Expression_Computer<DataType>::Known_TensOp_Initializer( shared_ptr<TensOp::TensOp<DataType>> TensOp_info  ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << " Expression_Computer::Expression_Computer<DataType>::TensOp_Data_Initializer" << endl;
  
   shared_ptr<Tensor_<DataType>> TensOp_data;    
   string opname = TensOp_info->name();

   cout <<" retrieve the operator from Bagel and convert it to tensor format " << endl;
   if ( opname == "H" ) {
     // This must make use of dfdist etc, similar to moint.cc here. 

     // get_op_in_Tensor_format( );
     // Check the blocks for sparsity
     // Feed information about sparsity bakc to 
   } else {
     
       cout << " build one operator tensor format " << endl;

   }

   return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
