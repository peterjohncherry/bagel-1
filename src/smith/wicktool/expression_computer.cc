
///////////////////////////////////////////////////////////////////////////////////////////////////////
Expression_Computer::Expression_Computer::Expression_Computer::(
///////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( scalar_results_map->find(Expression_name) != scalar_results_map->end() )  
    cout << "WARNING : You have already calculated this expression....." << Expression_name
    << " = " << scalar_results_map->at(Expression_name) << endl;
 
  shared_ptr<Expression<double>> Expr = Expr_Info->expression_map->at(Expression_name); 
  double result = 0.0;

  shared_ptr<TensOp_Computer::TensOp_Computer> TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer>(ref, Expr, range_conversion_map, TensOp_data_map);

  B_Gamma_Computer::B_Gamma_Computer B_Gamma_Machine( ref->ciwfn()->civectors(), range_conversion_map, Expr->GammaMap, Gamma_data_map, Sigma_data_map, CIvec_data_map );

  map<string , double > g_result_map;

///////////////////////////////////////////////////////////////////////////////////////////////////////
void Expression_Computer::Expression_Computer::Execute_Compute_List(string Expression_name) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "Expression_Computer::Expression_Computer::Execute_Compute_List(string expression_name ) " << endl;

  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(Expr->GammaMap) ) {
   
    string Gamma_name = AG_contrib.first;

    // Build A_tensor to hold sums of different A-tensors
    shared_ptr<Tensor_<double>> A_combined_data = make_shared<Tensor_<double>>( *(TensOp_Machine->Get_Bagel_IndexRanges(Expr->GammaMap->at(Gamma_name)->id_ranges)) );
    A_combined_data->allocate();
    A_combined_data->zero(); cout << " Gamma_name  = " << Gamma_name << endl;

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
            if (TensOp_data_map->find(A_contrib.first) != TensOp_data_map->end() ) {
              cout << endl; Print_Tensor(TensOp_data_map->at(A_contrib.first), A_contrib.first); cout << endl << endl << endl;
            } else {
              cout << A_contrib.first << " is not done" << endl;
            }
            shared_ptr<Tensor_<double>> A_contrib_reordered = TensOp_Machine->reorder_block_Tensor( A_contrib.first, make_shared<vector<int>>(A_contrib.second.id_order(qq)) );
            A_combined_data->ax_plus_y( (double)(A_contrib.second.factor(qq).first), A_contrib_reordered );
          }
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
        double tmp_result = Tensor_Arithmetic::Tensor_Arithmetic<double>::sum_tensor_elems( A_combined_data) ;
        g_result_map.emplace(Gamma_name, tmp_result) ;
        result += tmp_result ; 

      }
    }
  }
  cout << endl << endl;

  cout << "Contributions from different gamma terms " << endl;
  for ( auto elem : g_result_map ) 
    cout << elem.first << " :  " << elem.second << endl; 

  cout << endl; 
 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Expression_Computer::Expression_Computer::print_AContraction_list(shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list, string A_contrib_name ) {
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
bool Expression_Computer::Expression_Computer::check_AContrib_factors(AContribInfo& AC_info ) {
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
