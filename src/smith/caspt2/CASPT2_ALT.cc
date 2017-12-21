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
  cout << " CASPT2_ALT::CASPT2_ALT::CASPT2_ALT" << endl;

  civectors = orig_cpt2_in.info_->ciwfn()->civectors();

  T2_all     = orig_cpt2_in.t2all_;  
  lambda_all = orig_cpt2_in.lall_;   
  F_1el_all  = orig_cpt2_in.f1_;     cout << "F_1el_all->rms()" <<  F_1el_all->rms() << endl; 
  H_1el_all  = orig_cpt2_in.h1_;     cout << "H_1el_all->rms()" <<  H_1el_all->rms() << endl; 
  H_2el_all  = orig_cpt2_in.H_2el_;  cout << "H_2el_all->rms()" <<  H_2el_all->rms() << endl; 
 
  ncore   = orig_cpt2_in.info_->ncore();
  nclosed = orig_cpt2_in.info_->nclosed();
  nact    = orig_cpt2_in.info_->nact();
  nvirt   = orig_cpt2_in.info_->nvirt();
  maxtile = orig_cpt2_in.info_->maxtile();
  cimaxtile = 100000;

  TensOp_data_map    = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  Gamma_data_map     = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  scalar_results_map = make_shared<map<string, double>>();
  range_conversion_map = make_shared<map<string, shared_ptr<IndexRange>>>();

  CIvec_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    
  Sigma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    

  shared_ptr<vector<int>> states_of_interest = make_shared<vector<int>>( vector<int> { 0 } );
  set_target_info(states_of_interest) ;
  set_range_info(states_of_interest);
   
  Sys_Info = make_shared<System_Info<double>>(TargetsInfo, true);
  Expression_map = Sys_Info->expression_map;
  Expression_Machine = make_shared<Expression_Computer::Expression_Computer<double>>(  civectors, Expression_map, range_conversion_map, TensOp_data_map, 
                                                                                       Gamma_data_map, Sigma_data_map, CIvec_data_map );

  ////////////////// FOR TESTING ///////////////////////////////////////  
  shared_ptr<CASPT2::CASPT2> orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
  orig_cpt2->set_rdm(0,0);
  Smith_rdm1 = orig_cpt2->rdm1_; TensOp_data_map->emplace( "Smith_rdm1", Smith_rdm1 );
  Smith_rdm2 = orig_cpt2->rdm2_; TensOp_data_map->emplace( "Smith_rdm2", Smith_rdm2 );
  Smith_rdm3 = orig_cpt2->rdm3_; TensOp_data_map->emplace( "Smith_rdm3", Smith_rdm3 );
  Smith_rdm4 = orig_cpt2->rdm4_; TensOp_data_map->emplace( "Smith_rdm4", Smith_rdm4 );

  shared_ptr<vector<int>> gorder = make_shared<vector<int>>( vector<int> {  1, 2, 0, 3} ); 
  shared_ptr<Tensor_<double>> rdm2_reord = 
  Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( Smith_rdm2, gorder );
 
  vector<Index> id_block_rdm2 = { active_rng->range(0), active_rng->range(0), active_rng->range(0),active_rng->range(0)}; 
  vector<Index> id_block_rdm1 = { active_rng->range(0), active_rng->range(0) }; 

  unique_ptr<double[]> rdm2_data = rdm2_reord->get_block( id_block_rdm2 );
  unique_ptr<double[]> rdm1_data = Smith_rdm1->get_block( id_block_rdm1 );
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_range_info(shared_ptr<vector<int>> states_of_interest ) {
//////////////////////////////////////////////////////////////////////////////////////////////
  cout << "CASPT2_ALT::CASPT2_ALT::set_range_info " << endl;

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

  PT2_ranges_ = { *not_closed_rng, *not_closed_rng, *not_virtual_rng, *not_virtual_rng};
  PT2_ranges_herm_conj_ = { *not_closed_rng, *not_closed_rng, *not_virtual_rng, *not_virtual_rng};

  range_conversion_map->emplace("cor", closed_rng);//change the naming of the ranges from cor to clo... 
  range_conversion_map->emplace("act", active_rng);
  range_conversion_map->emplace("A", active_rng); //fudge, but leave for now
  range_conversion_map->emplace("vir", virtual_rng);
  range_conversion_map->emplace("free", free_rng);

  range_conversion_map->emplace("notcor", not_closed_rng);
  range_conversion_map->emplace("notact", not_active_rng);
  range_conversion_map->emplace("notvir", not_virtual_rng); 
 
  for ( int ii : *states_of_interest ) {

    range_conversion_map->emplace( get_civec_name( ii , civectors->data(ii)->det()->norb(), civectors->data(ii)->det()->nelea(), civectors->data(ii)->det()->neleb()),
                                                   make_shared<IndexRange>(civectors->data(ii)->det()->size(), cimaxtile ));  
    
    shared_ptr<IndexRange>  ci_index_ranges =  make_shared<IndexRange>(civectors->data(ii)->det()->size(), cimaxtile );
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
     TargetsInfo->add_state( civectors->data(state_num)->det()->nelea(), civectors->data(state_num)->det()->neleb(),
                             civectors->data(state_num)->det()->norb(), state_num ) ;
  
  return;
}

////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::solve() {
////////////////////////////////////////////////////////////////////
cout <<  " CASPT2_ALT::CASPT2_ALT::solve() " << endl;
  {
  vector<string> op_list = { "H" };
  vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
  
  double factor = 0.0;
  // Building all necessary expressions 
  int  num_states = 1; 
  vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
  for ( int ii = 0 ; ii != num_states; ii++) {
    for ( int jj = 0 ; jj != num_states; jj++) {
  
      for ( pair<vector<string>,double> BK_info : BK_info_list ) {
        Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "H check" ));
        for ( string Op_name : BK_info.first )  
          Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
      }
  
      string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
  
      Expression_Machine->Evaluate_Expression( expression_name );
 
    }
  }
  }
  {
  vector<string> op_list = { "L", "F" };
  vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
  
  double factor = 0.0;
  // Building all necessary expressions 
  int  num_states = 1; 
  vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
  for ( int ii = 0 ; ii != num_states; ii++) {
    for ( int jj = 0 ; jj != num_states; jj++) {
  
      for ( pair<vector<string>,double> BK_info : BK_info_list ) {
        Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "wicktool HE ones test" ));
        for ( string Op_name : BK_info.first )  
          Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
      }
  
      string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
  
      Expression_Machine->Evaluate_Expression( expression_name );
 
    }
  }
  }

 return;
} 
/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data( string op_name, string bra_name , string ket_name )  { 
/////////////////////////////////////////////////////////////////////////////////
cout << "CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data() " << endl;

  //Setting_data for TensOps
  if ( op_name  == "P" ) { 

    shared_ptr<vector<IndexRange>> H2el_ranges = make_shared<vector<IndexRange>>( H_2el_all->indexrange() );
    shared_ptr<Tensor_<double>> P_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( H2el_ranges);
    TensOp_data_map->emplace( "P", P_Tens );

  } else if ( op_name  == "Smith_rdms" ) { 
   
    TensOp_data_map->emplace("Smith_rdm1", Smith_rdm1 );
    TensOp_data_map->emplace("Smith_rdm2", Smith_rdm2 );

  } else if ( op_name  == "S" ) { 

    vector<IndexRange> act4_ranges = { *active_rng, *active_rng, *active_rng, *active_rng};
    shared_ptr<Tensor_<double>> S_Tens = Tensor_Arithmetic_Utils::get_sub_tensor( H_2el_all, act4_ranges );
    TensOp_data_map->emplace( "S", S_Tens );
    cout <<"---------------------------Energy_act_test----------------------------------" << endl;
    cout << "Energy_act S = " << Smith_rdm2->dot_product(S_Tens) << endl;
    cout <<"----------------------------------------------------------------------------" << endl;
    
  } else if ( op_name  == "Z" ) {
 
    shared_ptr<vector<IndexRange>> act4_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { *free_rng, *free_rng, *free_rng, *free_rng} );
    shared_ptr<Tensor_<double>> Z_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( act4_ranges, 1.0 );
    TensOp_data_map->emplace( "Z", Z_Tens );

  } else if ( op_name  == "G" ) {
 
    shared_ptr<vector<IndexRange>> act2_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { *not_closed_rng, *not_virtual_rng } );
    shared_ptr<Tensor_<double>> G_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( act2_ranges, 2.0 );
    TensOp_data_map->emplace( "G", G_Tens );

  } else if ( op_name  == "F" ) {
 
    shared_ptr<vector<IndexRange>> free4_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { *free_rng, *free_rng, *free_rng, *free_rng} );
    shared_ptr<Tensor_<double>> F_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( free4_ranges, 1.0 );
    TensOp_data_map->emplace( "F", F_Tens );

  } else if ( op_name  == "Y" ) {
 
    shared_ptr<vector<IndexRange>> act8_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng, *active_rng, *active_rng, *active_rng, *active_rng} );
    shared_ptr<Tensor_<double>> Y_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( act8_ranges, 2.0 );
    TensOp_data_map->emplace( "Y", Y_Tens );

  } else if ( op_name  == "H" ) { 

    Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( H_2el_all , 0.5  );
    TensOp_data_map->emplace("H" , H_2el_all);

  } else if ( op_name  == "h" ) { 

    TensOp_data_map->emplace("h" , H_1el_all);

  } else if ( op_name  == "T" ) { 

    shared_ptr<Tensor_<double>> TTens_data =  make_shared<Tensor_<double>>(PT2_ranges_); 
    TTens_data->allocate();
    Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( TTens_data , 1.0  );
    TensOp_data_map->emplace("T" , TTens_data);

  } else if ( op_name  == "L" ) { 

    shared_ptr<Tensor_<double>> LTens_data =  make_shared<Tensor_<double>>(PT2_ranges_herm_conj_); 
    LTens_data->allocate();
    Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( LTens_data , 1.0  );
    TensOp_data_map->emplace("L" , LTens_data);
  }

  return;
}
#endif

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
 //  shared_ptr<vector<string>> op_list4 = make_shared<vector<string>>(vector<string> { "U" });
 //  Sys_Info->Set_BraKet_Ops( op_list4, "<I|U|J>" ) ;
 //  shared_ptr<vector<string>> BK_list_U = make_shared<vector<string>>(vector<string> { "<I|U|J>" });
 //  Sys_Info->Build_Expression( BK_list_U, "U4_test") ;
 //
 // Sys_Info->Initialize_Tensor_Op_Info( "U" );
 // shared_ptr<vector<IndexRange>> act_ranges_4 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> U_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges_4 );
 // TensOp_data_map->emplace( "U", U_Tens );
 //
 // Sys_Info->Initialize_Tensor_Op_Info( "S" );
 // vector<IndexRange> act_ranges = { *active_rng, *active_rng, *active_rng, *active_rng };
 // shared_ptr<Tensor_<double>> Hact = Tensor_Arithmetic_Utils::get_sub_tensor( H_2el_all, act_ranges);

 // <proj_jst|H|0_K> set to sall in ms-caspt2
  
  // Get < M | W H | N  >

  // Input into linear equation solver to get T amplitudes
  
  // Should get \sum < M | W ( f- E_{L}+E_{s} ) T_{LN} | N>  in solver
 
  // Construct effective Hamiltonian
  
  // symmetrize and diagonalize effective Hmailtonian
// 
// Sys_Info->Initialize_Tensor_Op_Info( "X" );
// vector<IndexRange> X_ranges = { *not_virtual_rng, *not_virtual_rng, *not_closed_rng , *not_closed_rng} ;
// shared_ptr<Tensor_<double>> XTens_data =  make_shared<Tensor_<double>>( X_ranges ); 
// XTens_data->allocate();
// Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( XTens_data , 1.0  );
// TensOp_data_map->emplace( "X", XTens_data);
// TensOp_data_map->emplace( "S", Hact );
//
//    shared_ptr<vector<IndexRange>> P_cvvc_ranges = make_shared<vector<IndexRange>>(vector<IndexRange>{*closed_rng , *virtual_rng, *virtual_rng, *closed_rng}) ;
//    shared_ptr<Tensor_<double>> P_cvvc = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( P_cvvc_ranges );
//    Print_Tensor( P_cvvc, "P_cvvc tensor" ) ; cout << endl << endl; 
//    shared_ptr<Tensor_<double>> P_Tens03 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( P_cvvc, make_pair(0,3) );
//    Print_Tensor( P_Tens03, "P_cvvc_03_contract" ) ; cout << endl << endl; 
//    shared_ptr<Tensor_<double>> P_Tens12 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( P_cvvc, make_pair(1,2) );
//    Print_Tensor( P_Tens12, "P_cvvc_12_contract" ) ; cout << endl << endl; 
 // shared_ptr<vector<IndexRange>> act_ranges2 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> Test_Tens2 = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges2 );
 //
 // Print_Tensor( Test_Tens2 , "TT2" ) ; cout << endl << endl;
 //
 // shared_ptr<vector<IndexRange>> act_ranges4 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> Test_Tens4 = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges4 );
 //
 // Print_Tensor( Test_Tens4 , "TT4" ) ; cout << endl << endl;
 // shared_ptr<vector<IndexRange>> act_ranges6 = make_shared<vector<IndexRange>>( vector<IndexRange> { *active_rng, *active_rng, *active_rng, *active_rng, *active_rng, *active_rng } );
 // shared_ptr<Tensor_<double>> Test_Tens6 = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( act_ranges6 );
 //
 // vector<int> vec0 = { 0 };
 // vector<int> vec1 = { 1 };
 // vector<int> vec3 = { 3 };
 // vector<int> vec01 = { 0,1 };
 // vector<int> vec23 = { 2,3 };
 // vector<int> vec32 = { 3,2 };
 // vector<int> vec0123 = { 0,1,2,3 };
 // 
 // pair<vector<int>,vector<int>> ctrs_todo_00( vec0, vec0 );
 // shared_ptr<Tensor_<double>> T2T2_00 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens2, Test_Tens2, ctrs_todo_00 );
 // Print_Tensor( T2T2_00 , "T2T2_00" ) ; cout << endl << endl;
 //  
 // pair<vector<int>,vector<int>> ctrs_todo_01( vec0, vec1 );
 // shared_ptr<Tensor_<double>> T2T2_01 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens2, Test_Tens2, ctrs_todo_01 );
 // Print_Tensor( T2T2_01 , "T2T2_01" ) ; cout << endl << endl;
 //  
 // pair<vector<int>,vector<int>> ctrs_todo_10( vec1, vec0 );
 // shared_ptr<Tensor_<double>> T2T2_10 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens2, Test_Tens2, ctrs_todo_10 );
 // Print_Tensor( T2T2_10 , "T2T2_10" ) ; cout << endl << endl;
 //  
 // pair<vector<int>,vector<int>> ctrs_todo_11( vec1, vec1 );
 // shared_ptr<Tensor_<double>> T2T2_11 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens2, Test_Tens2, ctrs_todo_11 );
 // Print_Tensor( T2T2_11 , "T2T2_11" ) ; cout << endl << endl;
 // 
 // pair<vector<int>,vector<int>> ctrs_todo_30( vec3, vec0 );
 // shared_ptr<Tensor_<double>> T4T2_30 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens4, Test_Tens2, ctrs_todo_30 );
 // Print_Tensor( T4T2_30 , "T4T2_30" ) ; cout << endl << endl;
 //
 // pair<vector<int>,vector<int>> ctrs_todo_3201( vec32, vec01 );
 // shared_ptr<Tensor_<double>> T4T2_3201 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens4, Test_Tens2, ctrs_todo_3201 );
 // Print_Tensor( T4T2_3201 , "T4T2_3201" ) ; cout << endl << endl;
 //
 // pair<vector<int>,vector<int>> ctrs_todo_0101( vec01, vec01 );
 // shared_ptr<Tensor_<double>> T4T4_0101 =   Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( Test_Tens4, Test_Tens4, ctrs_todo_0101 );
 // Print_Tensor( T4T4_0101 , "T4T4_0101" ) ; cout << endl << endl;
/////////////////////////////////////////////////////////////////////////////////
 

