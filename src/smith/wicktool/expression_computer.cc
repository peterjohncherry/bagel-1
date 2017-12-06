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
  B_Gamma_Computer::B_Gamma_Computer B_Gamma_Machine( civectors, range_conversion_map, Expr->GammaMap, Gamma_data_map, Sigma_data_map, CIvec_data_map );

  double result = 0.0;
  map< string, double > g_result_map;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(Expr->GammaMap) ) {

    string Gamma_name = AG_contrib.first; 

    cout << " Gamma_name  = " << Gamma_name << endl; print_vector( *(Expr->GammaMap->at(Gamma_name)->id_ranges), "id_ranges" );

    shared_ptr<Tensor_<double>> A_combined_data;
    // Build A_tensor to hold sums of different A-tensors.
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
    if ( (A_contrib_loc != Expr->G_to_A_map->end()) &&  (A_contrib_loc->second->size() != 0) ) {
   
      for ( auto  A_contrib_map_elem : *A_contrib_loc->second ) {
      
        string  A_contrib_name  = A_contrib_map_elem.first;    
        AContribInfo A_contrib  = A_contrib_map_elem.second;    

        if (check_AContrib_factors(A_contrib))
          continue;
      
	print_AContraction_list(Expr->ACompute_map->at(A_contrib_name), A_contrib_name);
        TensOp_Machine->Calculate_CTP(A_contrib_name);

        if ( Gamma_name != "ID" ) {
          for ( int qq = 0 ; qq != A_contrib.id_orders.size(); qq++){
       
            if ( TensOp_data_map->find(A_contrib_name) != TensOp_data_map->end() ) { cout << A_contrib_name << " found in map" << endl;
           
              Print_Tensor(TensOp_data_map->at(A_contrib_name), A_contrib_name); cout << endl << endl << endl; 
//              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), TensOp_data_map->at(A_contrib_name) );
              print_vector(A_contrib.id_order(qq), "A_contrib_order("+to_string(qq)+")"); cout << endl;
              shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib.id_order(qq)) );
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_reordered );
            
            } else { cout << A_contrib_name << " not found in map; must be formed from direct product" << endl; //TODO there should be a way to avoid this
            
              cout << A_contrib_name << " Tensor is decomposed, do contraction in parts"  << endl;
              shared_ptr<vector<shared_ptr<CtrTensorPart<double>>>> CTP_vec = Expr->CMTP_map->at(A_contrib_name)->CTP_vec ;
              vector<string> sub_tensor_names(CTP_vec->size()); 
              for ( int rr = 0 ; rr != CTP_vec->size() ; rr++ )
                sub_tensor_names[rr] = CTP_vec->at(rr)->myname();

              shared_ptr<Tensor_<double>> A_contrib_data = TensOp_Machine->direct_product_tensors( sub_tensor_names );//TODO fix so uses piecewise contraction where possible 
              TensOp_data_map->emplace(A_contrib_name, A_contrib_data );
              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_data );
              print_vector(A_contrib.id_order(qq), "A_contrib_order("+to_string(qq)+")"); cout << endl;
//              shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib_name, make_shared<vector<int>>(A_contrib.id_order(qq)) );
//              A_combined_data->ax_plus_y( (double)(A_contrib.factor(qq).first), A_contrib_reordered );

            }
          }

        } else {
  
          cout << "A_contrib_name = " << A_contrib_name ; cout.flush();

          if ( TensOp_data_map->find(A_contrib_name) == TensOp_data_map->end() ) {cout << " not yet in map, must form from direct product" << endl;
            shared_ptr<vector<shared_ptr<CtrTensorPart<double>>>> CTP_vec = Expr->CMTP_map->at(A_contrib_name)->CTP_vec ;
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
        cout << "A_combined_data->dot_product( Gamma_data_map->at(" << Gamma_name << ") ) = " << tmp_result  << endl;
        g_result_map.emplace(Gamma_name, tmp_result) ;
        result += tmp_result;

        if ( Gamma_data_map->at(Gamma_name)->rank() == 2 ) {
          
           Gamma_data_map->emplace( "Gamma2", Gamma_data_map->at(Gamma_name) );

        } else if ( Gamma_data_map->at(Gamma_name)->rank() == 4 ) {
           

          shared_ptr<vector<int>> gorder = make_shared<vector<int>>( vector<int> { 0, 3, 1, 2, } ); 
          shared_ptr<Tensor_<double>> Gamma4_reord = 
          Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( Gamma_data_map->at(Gamma_name), gorder );
          vector<IndexRange> gids = Gamma4_reord->indexrange(); 

          vector<Index> g2_id_blocks = { gids[0].range(0), gids[1].range(0) }; 
          vector<Index> g4_id_blocks = { gids[0].range(0), gids[1].range(0), gids[2].range(0), gids[3].range(0) }; 

          unique_ptr<double[]> gamma2_data = Gamma_data_map->at( "Gamma2" )->get_block(g2_id_blocks);
          unique_ptr<double[]> gamma4_data = Gamma4_reord->get_block(g4_id_blocks); 

          int stride2 = g4_id_blocks[1].size()*g4_id_blocks[0].size();
          int stride3 = g4_id_blocks[2].size()*stride2;

          double* gamma4_data_ptr = gamma4_data.get();
          double* gamma2_data_ptr = gamma2_data.get();
          for ( int qq = 0; qq != g4_id_blocks[3].size(); qq++ ) {
            daxpy_( stride2, -1.0, gamma2_data_ptr, 1, gamma4_data_ptr, 1 );
            gamma4_data_ptr += stride3+stride2;
          }
          gamma4_data_ptr = gamma4_data.get();

          Gamma4_reord->put_block( gamma4_data, g4_id_blocks );
          
          shared_ptr<vector<int>> gorder_back = make_shared<vector<int>>( vector<int> {  0, 2, 3, 1} ); 

          shared_ptr<Tensor_<double>> rdm2_from_Gamma4 = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( Gamma4_reord, gorder_back );

          Print_Tensor( rdm2_from_Gamma4, "rdm2 from gamma4 wicktool" ); cout << endl << endl << endl;    
          
                  
          cout <<"---------------------------Energy_act_test----------------------------------" << endl;
          cout << "Smith_rdm2_dot_2el_op   = " << TensOp_data_map->at("Smith_rdm2")->dot_product(TensOp_data_map->at("S")) << endl;
          cout << "rdm2_from_G4_dot_2el_op = " << rdm2_from_Gamma4->dot_product(TensOp_data_map->at("S")) << endl;
          cout << "Smith_rdm2_dot_A_contrib   = " << TensOp_data_map->at("Smith_rdm2")->dot_product(A_combined_data) << endl;
          cout << "rdm2_from_G4_dot_A_contrib = " << rdm2_from_Gamma4->dot_product(A_combined_data) << endl;
          cout <<"----------------------------------------------------------------------------" << endl;

        }

      } else {

//        Print_Tensor( A_combined_data, " A_combined_data for 1D " ) ; cout << endl;

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
template class Expression_Computer::Expression_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
 
#endif
