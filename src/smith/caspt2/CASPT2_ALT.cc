//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/caspt2/CASPT2_ALT.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::TensOp_Computer; 
using namespace bagel::SMITH::Gamma_Computer; 
using namespace Tensor_Arithmetic_Utils;

////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in ) { 
////////////////////////////////////////////////////////////////////
  orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
  ref = orig_cpt2_in.info_;

  cc_ = ref->ciwfn()->civectors();
  det_ = ref->ciwfn()->civectors()->det(); //remove this ASAP

  T2_all     = orig_cpt2->t2all_;  
  lambda_all = orig_cpt2->lall_;   
  F_1el_all  = orig_cpt2->f1_;     cout << "F_1el_all->rms()" <<  F_1el_all->rms() << endl; 
  H_1el_all  = orig_cpt2->h1_;     cout << "H_1el_all->rms()" <<  H_1el_all->rms() << endl; 
  H_2el_all  = orig_cpt2->H_2el_;  cout << "H_2el_all->rms()" <<  H_2el_all->rms() << endl; 
 
  ncore   = ref->ncore();
  nclosed = ref->nclosed();
  nact    = ref->nact();
  nvirt   = ref->nvirt();
  maxtile = ref->maxtile();
  cimaxtile = 100000;

  TensOp_data_map    = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  Gamma_data_map     = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  GammaMap           = make_shared<map<string, shared_ptr<GammaInfo>>>();
  scalar_results_map = make_shared<map<string, double>>();

  CIvec_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    
  Sigma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    

  shared_ptr<vector<int>> states_of_interest = make_shared<vector<int>>( vector<int> { 0 } );
  set_target_info(states_of_interest) ;
  set_range_info(states_of_interest);
   
  Expr_Info = make_shared<System_Info<double>>(TargetsInfo, true);
  Expr_Info_map = Expr_Info->expression_map;

  /////////////////////////////////////////////////////////  
  orig_cpt2->set_rdm(0,0);
  Smith_rdm1 = orig_cpt2->rdm1_;
  Smith_rdm2 = orig_cpt2->rdm2_;
  Smith_rdm3 = orig_cpt2->rdm3_;
  Smith_rdm4 = orig_cpt2->rdm4_;
  cout << endl << endl << endl << "-------------------------------------------------------------------------------" << endl;
  Tensor_Arithmetic_Utils::Print_Tensor(Smith_rdm1, "Smith rdm1" );
  cout << endl << "-------------------------------------------------------------------------------" << endl;
  Tensor_Arithmetic_Utils::Print_Tensor(Smith_rdm2, "Smith rdm2" );
  cout << endl << endl << endl << "-------------------------------------------------------------------------------" << endl;
 
 shared_ptr<Tensor_<double>> Smith_rdm1_from_rdm2 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( Smith_rdm2 , make_pair(1,2));
  Tensor_Arithmetic_Utils::Print_Tensor( Smith_rdm1_from_rdm2, "Smith rdm1 from rdm2" );
  cout << endl << endl << endl << "-------------------------------------------------------------------------------" << endl;
  
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_range_info(shared_ptr<vector<int>> states_of_interest ) {
//////////////////////////////////////////////////////////////////////////////////////////////

  cout << maxtile << endl;
  int maxtile_buff = maxtile; 
  closed_rng  =  make_shared<IndexRange>(IndexRange(nclosed-ncore, maxtile, 0, ncore));
  active_rng  =  make_shared<IndexRange>(IndexRange(nact, min(10,maxtile_buff), closed_rng->nblock(), ncore + closed_rng->size()));
  virtual_rng =  make_shared<IndexRange>(IndexRange(nvirt, maxtile, closed_rng->nblock()+ active_rng->nblock(), ncore+closed_rng->size()+active_rng->size()));
  free_rng    = make_shared<IndexRange>(*closed_rng);
  free_rng->merge(*active_rng);
  free_rng->merge(*virtual_rng);

  not_closed_rng  =  make_shared<IndexRange>(*active_rng); not_closed_rng->merge(*virtual_rng);
  not_active_rng  =  make_shared<IndexRange>(*closed_rng); not_active_rng->merge(*virtual_rng);
  not_virtual_rng =  make_shared<IndexRange>(*closed_rng); not_virtual_rng->merge(*active_rng);

  range_conversion_map = make_shared<map<string, shared_ptr<IndexRange>>>();
  range_conversion_map->emplace("cor", closed_rng);//change the naming of the ranges from cor to clo... 
  range_conversion_map->emplace("act", active_rng);
  range_conversion_map->emplace("A", active_rng); //fudge, but leave for now
  range_conversion_map->emplace("vir", virtual_rng);
  range_conversion_map->emplace("free", free_rng);

  range_conversion_map->emplace("notcor", not_closed_rng);
  range_conversion_map->emplace("notact", not_active_rng);
  range_conversion_map->emplace("notvir", not_virtual_rng); 
 
  for ( int ii : *states_of_interest ) {

    range_conversion_map->emplace( get_civec_name( ii , cc_->data(ii)->det()->norb(), cc_->data(ii)->det()->nelea(), cc_->data(ii)->det()->neleb()),
                                                   make_shared<IndexRange>(cc_->data(ii)->det()->size(), cimaxtile ));  
    
    shared_ptr<IndexRange>  ci_index_ranges =  make_shared<IndexRange>(cc_->data(ii)->det()->size(), cimaxtile );
    cout << "cirngs = [ ";  for (auto irng : ci_index_ranges->range()) { cout << irng.size()  << " "; }; cout << "] " << endl;

  }    

  return;

}
//////////////////////////////////////////////////////////////////////////////////////////////
// defining this so we don't need to put determinants class into gamma generator etc.
// At the moment this is quite silly; the determinant space is the same for every space,
// however, this is not necessarily the case. This function will need to do something 
// different in the relativistic case.
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_target_info( shared_ptr<vector<int>> states_of_interest) {
//////////////////////////////////////////////////////////////////////////////////////////////
  TargetsInfo = make_shared<StatesInfo<double>> ( *states_of_interest ) ;
  
  for ( int state_num : *states_of_interest ) 
     TargetsInfo->add_state( cc_->data(state_num)->det()->nelea(), cc_->data(state_num)->det()->neleb(),
                             cc_->data(state_num)->det()->norb(), state_num ) ;
  
  return;
}

////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::solve() {
////////////////////////////////////////////////////////////////////
cout <<  " CASPT2_ALT::CASPT2_ALT::solve() " << endl;

  Construct_Tensor_Ops();
  
  Build_Compute_Lists();

  Execute_Compute_List("Hact_test");
  Execute_Compute_List("R6_test");
//  Execute_Compute_List("Q8_test");
//  Execute_Compute_List("U4_test");

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  
  // Get < M | W H | N  >

  // Input into linear equation solver to get T amplitudes
  
  // Should get \sum < M | W ( f- E_{L}+E_{s} ) T_{LN} | N>  in solver
 
  // Construct effective Hamiltonian
  
  // symmetrize and diagonalize effective Hmailtonian

}
/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Construct_Tensor_Ops() { 
/////////////////////////////////////////////////////////////////////////////////
cout << "CASPT2_ALT::CASPT2_ALT::Construct_Tensor_Ops() " << endl;

  // Setting_data for TensOps
 // Expr_Info->Initialize_Tensor_Op_Info( "U" );
 // shared_ptr<vector<IndexRange>> act_ranges_4 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> U_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges_4 );
 // TensOp_data_map->emplace( "U", U_Tens );

  Expr_Info->Initialize_Tensor_Op_Info( "R" );
  shared_ptr<vector<IndexRange>> act_ranges_6 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng, *active_rng, *active_rng } );
  shared_ptr<Tensor_<double>> R_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges_6);
  TensOp_data_map->emplace( "R", R_Tens );

  Expr_Info->Initialize_Tensor_Op_Info( "S" );
  vector<IndexRange> act_ranges = { *active_rng, *active_rng, *active_rng, *active_rng };
  shared_ptr<Tensor_<double>> Hact = Tensor_Arithmetic_Utils::get_sub_tensor( H_2el_all, act_ranges);
  TensOp_data_map->emplace( "S", Hact );
     
  TensOp_data_map->emplace("H" , H_2el_all);
  TensOp_data_map->emplace("h" , H_1el_all);
   
  //vector<IndexRange>                  Y_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng } ;
  //vector<IndexRange>                  X_bagel_ranges = { *not_virtual_rng, *not_virtual_rng, *not_closed_rng , *not_closed_rng} ;
  vector<IndexRange>                  T_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng} ;
  vector<IndexRange>                  L_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng} ;
  shared_ptr<Tensor_<double>> TTens_data =  make_shared<Tensor_<double>>(T_bagel_ranges); 
  TTens_data->allocate();
  TTens_data->zero();
  TensOp_data_map->emplace("T" , TTens_data);
   
  shared_ptr<Tensor_<double>> LTens_data =  make_shared<Tensor_<double>>(L_bagel_ranges); 
  LTens_data->allocate();
  LTens_data->zero();
  TensOp_data_map->emplace("L" , LTens_data);
     

/////////////////////////////////////////////////////////////////////////////////
  /*tests; keep for now so can check if anything starts to go wrong */
 // shared_ptr<vector<IndexRange>> act_ranges_1 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng } );
 // shared_ptr<Tensor_<double>> Vec_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges_1 );
 // cout << endl ; Print_Tensor( Vec_Tens, "Vec_Tens" ); cout << endl << endl << endl;
 // shared_ptr<Tensor_<double>> UV_contracted_0 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_tensor_with_vector( U_Tens, Vec_Tens, 0);
 // cout << endl ; Print_Tensor( UV_contracted_0, "UV_contracted_0" ); cout << endl << endl << endl;
 //
 // shared_ptr<vector<IndexRange>> act_ranges_2 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> Mat_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges_2 );
 // cout << endl ; Print_Tensor( Mat_Tens, "Mat_Tens" ); cout << endl << endl << endl;
 // shared_ptr<Tensor_<double>> MV_contracted_0 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_tensor_with_vector( Mat_Tens, Vec_Tens, 0);
 // cout << endl ; Print_Tensor( MV_contracted_0, "MV_contracted_0" ); cout << endl << endl << endl;
 //
 // shared_ptr<vector<IndexRange>> actvir_ranges_2 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *virtual_rng } );
 // shared_ptr<Tensor_<double>> Rect_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( actvir_ranges_2 );
 // cout << endl ; Print_Tensor( Rect_Tens, "Rect_Tens" ); cout << endl << endl << endl;
 // shared_ptr<Tensor_<double>> RectV_contracted_0 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_tensor_with_vector( Rect_Tens, Vec_Tens, 0);
 // cout << endl ; Print_Tensor( RectV_contracted_0, "RectV_contracted_0" ); cout << endl << endl << endl;
/////////////////////////////////////////////////////////////////////////////////
 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Build_Compute_Lists() { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
                                                                                     
  shared_ptr<vector<string>> op_list1 = make_shared<vector<string>>(vector<string> { "S" });
  Expr_Info->Set_BraKet_Ops( op_list1, "<I|H_act|J>" ) ;
  shared_ptr<vector<string>> BK_list_S = make_shared<vector<string>>(vector<string> { "<I|H_act|J>" });
  Expr_Info->Build_Expression( BK_list_S, "Hact_test") ;

  shared_ptr<vector<string>> op_list3 = make_shared<vector<string>>(vector<string> { "R" });
  Expr_Info->Set_BraKet_Ops( op_list3, "<I|R|J>" ) ;
  shared_ptr<vector<string>> BK_list_R = make_shared<vector<string>>(vector<string> { "<I|R|J>" });
  Expr_Info->Build_Expression( BK_list_R, "R6_test") ;

//  shared_ptr<vector<string>> op_list4 = make_shared<vector<string>>(vector<string> { "U" });
//  Expr_Info->Set_BraKet_Ops( op_list4, "<I|U|J>" ) ;
//  shared_ptr<vector<string>> BK_list_U = make_shared<vector<string>>(vector<string> { "<I|U|J>" });
//  Expr_Info->Build_Expression( BK_list_U, "U4_test") ;

  return ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(string Expression_name) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(string expression_name ) " << endl;

  if ( scalar_results_map->find(Expression_name) != scalar_results_map->end() )  
    cout << "WARNING : You have already calculated this expression....." << Expression_name
    << " = " << scalar_results_map->at(Expression_name) << endl;
 
  shared_ptr<Expression<double>> Expr = Expr_Info->expression_map->at(Expression_name); 
  double result = 0.0;

  shared_ptr<TensOp_Computer::TensOp_Computer> TensOp_Machine = make_shared<TensOp_Computer::TensOp_Computer>( Expr->ACompute_map, Expr->CTP_map, range_conversion_map, TensOp_data_map);

  B_Gamma_Computer::B_Gamma_Computer B_Gamma_Machine( ref->ciwfn()->civectors(), range_conversion_map, Expr->GammaMap, Gamma_data_map, Sigma_data_map, CIvec_data_map );

  map<string , double > g_result_map;

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

  double reference_result;
 
  if ( Expression_name == "Hact_test" ) {
    reference_result = Smith_rdm2->dot_product( TensOp_data_map->at("S"));
  } else if ( Expression_name == "U4_test" ) {
    reference_result = Smith_rdm2->dot_product( TensOp_data_map->at("U"));
  } else if ( Expression_name == "R6_test" ) {
    reference_result = Smith_rdm3->dot_product( TensOp_data_map->at("R"));
  } else if ( Expression_name == "Q8_test" ) {
    cout << "Smith_rdm4->norm() = " << Smith_rdm4->norm() << "     TensOp_data_map->at(\"Q\")->norm() = " << TensOp_data_map->at("Q")->norm() << endl;
    reference_result = Smith_rdm4->dot_product( TensOp_data_map->at("Q"));
  }
 

  cout << "==================================== RESULTS for "<< Expression_name << "===================" << endl << endl;
  cout << Expression_name << " = " << result << endl ;
  cout << "reference result  = " << reference_result << endl;
  cout << "difference = " << result-reference_result << endl << endl;
  cout << "=========================================================================================================" << endl << endl << endl;
  scalar_results_map->emplace( Expression_name, result );
  
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::print_AContraction_list(shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list, string A_contrib_name ) {
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
bool CASPT2_ALT::CASPT2_ALT::check_AContrib_factors(AContribInfo& AC_info ) {
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
#endif

// Hack to test excitation operators
//  shared_ptr<vector<IndexRange>> X_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { *not_virtual_rng, *not_virtual_rng, *not_closed_rng,  *not_closed_rng});
//  shared_ptr<vector<IndexRange>> Y_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { *not_closed_rng,  *not_closed_rng,  *not_virtual_rng, *not_virtual_rng});
//  shared_ptr<Tensor_<double>> XTens_data = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(X_ranges , 1.0 );
//  TensOp_data_map->emplace("X" , XTens_data);
//  shared_ptr<Tensor_<double>> YTens_data = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(Y_ranges , 1.0 );
//  TensOp_data_map->emplace("Y" , YTens_data);

