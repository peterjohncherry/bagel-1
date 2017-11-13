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
using namespace bagel::SMITH::Equation_Computer; 
using namespace bagel::SMITH::Gamma_Computer; 

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

  CTP_map            = make_shared<map<string, shared_ptr<CtrTensorPart<double>>>>();
  Data_map           = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  Gamma_data_map     = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  GammaMap           = make_shared<map<string, shared_ptr<GammaInfo>>>();
  scalar_results_map = make_shared<map<string, double>>();

  CIvec_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    
  Sigma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    
  Determinants_map = make_shared<std::map<std::string, std::shared_ptr<const Determinants>>>(); 

  shared_ptr<vector<int>> states_of_interest = make_shared<vector<int>>( vector<int> { 0 } );
  set_target_info(states_of_interest) ;
  set_range_info(states_of_interest);
   
  Expr_Info = make_shared<Expression_Info<double>>(TargetsInfo, true);
  Expr_Info_map = Expr_Info->expression_map;

  /////////////////////////////////////////////////////////  
  orig_cpt2->set_rdm(0,0);
  Smith_rdm1 = orig_cpt2->rdm1_;
  Smith_rdm2 = orig_cpt2->rdm2_;
  Smith_rdm3 = orig_cpt2->rdm3_;
  Smith_rdm4 = orig_cpt2->rdm4_;
  
  cout << endl << "-------------------------------------------------------------------------------" << endl;
  Tensor_Arithmetic_Utils::Print_Tensor(Smith_rdm1, "Smith rdm1" );
  cout << endl << "-------------------------------------------------------------------------------" << endl;
  Tensor_Arithmetic_Utils::Print_Tensor(Smith_rdm2, "Smith rdm2" );
  cout << endl << "-------------------------------------------------------------------------------" << endl;
  shared_ptr<Tensor_<double>> Smith_rdm1_from_rdm2 = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( Smith_rdm2 , make_pair(2,3));
  
  Tensor_Arithmetic_Utils::Print_Tensor( Smith_rdm1_from_rdm2, "Smith rdm1 from rdm2" );
  cout << endl << "-------------------------------------------------------------------------------" << endl;
  pair<int,int> ctr_tmp = make_pair(2,3);
  shared_ptr<Tensor_<double>> Smith_rdm1_from_rdm2_new = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor_new( Smith_rdm2 , ctr_tmp);
  Tensor_Arithmetic_Utils::Print_Tensor( Smith_rdm1_from_rdm2_new, "Smith rdm1 from rdm2 new" );
  cout << endl << "-------------------------------------------------------------------------------" << endl;

  /////////////////////////////////////////////////////////  
  
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

  cout <<" cc_->data("<<ii<<")->det()->size() = " << cc_->data(ii)->det()->size() <<  endl;
  cout <<" maxtile = " << maxtile << endl;
  cout <<" cimaxtile = " << cimaxtile << endl;
  
  shared_ptr<IndexRange>  ci_index_ranges =  make_shared<IndexRange>(cc_->data(ii)->det()->size(), cimaxtile );
  cout <<" ci_index_ranges->nblock()       = " << ci_index_ranges->nblock()     <<  endl;
  cout <<" ci_index_ranges->size()         = " << ci_index_ranges->size()       <<  endl;
  cout <<" ci_index_ranges->range().size() = " << ci_index_ranges->range().size() <<  endl;
  cout <<" ci_index_ranges->range().size() = " << ci_index_ranges->range().size() <<  endl;
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

  Execute_Compute_List("test_case");

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
 
  //spinfree orbital ranges
  vector<string> free     = {"cor", "act", "vir"};
  vector<string> not_core = {"act", "vir"};
  vector<string> not_act  = {"cor", "vir"};
  vector<string> not_virt = {"cor", "act"};
  vector<string> core = {"cor"};
  vector<string> act  = {"act"};
  vector<string> virt = {"vir"};

  /* ---- 2el Excitation Op  ----  */
  pair<double,double>                 X_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          X_idxs = make_shared<vector<string>>(vector<string> {"X0", "X1", "X2", "X3"});
  shared_ptr<vector<bool>>            X_aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true }); 
  shared_ptr<vector<vector<string>>>  X_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_virt, not_virt, not_core, not_core }); 
  vector<IndexRange>                  X_bagel_ranges = { *not_virtual_rng, *not_virtual_rng, *not_closed_rng , *not_closed_rng,} ;
  shared_ptr<double>                  X_dummy_data;
  string                              X_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  X_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  X_constraints = { &Expression_Info<double>::Expression_Info::NotAllAct };

  shared_ptr<TensOp<double>> XTens = Expr_Info->Build_TensOp("X", X_dummy_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("X", XTens);

  /* ---- 2el Excitation Op dag  ----  */
  pair<double,double>                 Y_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          Y_idxs = make_shared<vector<string>>(vector<string> {"Y0", "Y1", "Y2", "Y3"});
  shared_ptr<vector<bool>>            Y_aops = make_shared<vector<bool>>(vector<bool>  { true, true, false, false }); 
  shared_ptr<vector<vector<string>>>  Y_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt }); 
  vector<IndexRange>                  Y_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng } ;
  shared_ptr<double>                  Y_dummy_data;
  string                              Y_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  Y_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  Y_constraints = { &Expression_Info<double>::Expression_Info::NotAllAct };

  shared_ptr<TensOp<double>> YTens = Expr_Info->Build_TensOp("Y", Y_dummy_data, Y_idxs, Y_aops, Y_idx_ranges, Y_symmfuncs, Y_constraints, Y_factor, Y_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("Y", YTens);


  /* ---- H Tensor ----  */
  pair<double,double>                H_factor = make_pair(0.5,0.5);
  shared_ptr<vector<string>>         H_idxs = make_shared<vector<string>>(vector<string> {"H0", "H1", "H2", "H3"});
  shared_ptr<vector<bool>>           H_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
  shared_ptr<vector<vector<string>>> H_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
  shared_ptr<double>                 H_dummy_data;
  string                             H_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> H_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  H_constraints = {  &Expression_Info<double>::Expression_Info::always_true };

  shared_ptr<TensOp<double>> HTens = Expr_Info->Build_TensOp("H", H_dummy_data, H_idxs, H_aops, H_idx_ranges, H_symmfuncs, H_constraints, H_factor, H_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("H", HTens);

  Data_map->emplace("H" , H_2el_all);

  /* ---- h Tensor ----  */
  pair<double,double>                h_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>         h_idxs = make_shared<vector<string>>(vector<string> {"h0", "h1"});
  shared_ptr<vector<bool>>           h_aops = make_shared<vector<bool>>(vector<bool>  {true, false}); 
  shared_ptr<vector<vector<string>>> h_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free }); 
  shared_ptr<double>                 h_dummy_data;
  string                             h_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> h_symmfuncs = Expr_Info->set_1el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  h_constraints = {  &Expression_Info<double>::Expression_Info::always_true };

  shared_ptr<TensOp<double>> hTens = Expr_Info->Build_TensOp("h", h_dummy_data, h_idxs, h_aops, h_idx_ranges, h_symmfuncs, h_constraints, h_factor, h_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("h", hTens);

  Data_map->emplace("h" , H_1el_all);
  
  /* ---- T Tensor ----  */
  pair<double,double>                 T_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          T_idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
  shared_ptr<vector<bool>>            T_aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} ); 
  shared_ptr<vector<vector<string>>>  T_idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });   
  vector<IndexRange>                  T_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng} ;
  shared_ptr<double>                  T_dummy_data ;
  string                              T_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> T_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)> T_constraints = {  &Expression_Info<double>::Expression_Info::NotAllAct };

  shared_ptr<TensOp<double>> TTens = Expr_Info->Build_TensOp("T", T_dummy_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("T", TTens);

  shared_ptr<Tensor_<double>> TTens_data =  make_shared<Tensor_<double>>(T_bagel_ranges); 
  TTens_data->allocate();
  TTens_data->zero();
  Data_map->emplace("T" , TTens_data);
 

   /* ---- L Tensor ----  */
  pair<double,double>                 L_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>          L_idxs = make_shared<vector<string>>(vector<string> {"L0", "L1", "L2", "L3"});
  shared_ptr<vector<bool>>            L_aops = make_shared<vector<bool>>(vector<bool>  { false, false, true, true }); 
  shared_ptr<vector<vector<string>>>  L_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { not_virt, not_virt, not_core, not_core }); 
  vector<IndexRange>                  L_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng} ;
  shared_ptr<double>                  L_dummy_data;
  string                              L_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  L_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  L_constraints = { &Expression_Info<double>::Expression_Info::NotAllAct };

  shared_ptr<TensOp<double>> LTens = Expr_Info->Build_TensOp("L", L_dummy_data, L_idxs, L_aops, L_idx_ranges, L_symmfuncs, L_constraints, L_factor, L_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("L", LTens);

  shared_ptr<Tensor_<double>> LTens_data =  make_shared<Tensor_<double>>(L_bagel_ranges); 
  LTens_data->allocate();
  LTens_data->zero();
  Data_map->emplace("L" , LTens_data);


  //Hack to test excitation operators
  shared_ptr<vector<IndexRange>> X_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { *not_virtual_rng, *not_virtual_rng, *not_closed_rng,  *not_closed_rng});
  shared_ptr<vector<IndexRange>> Y_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { *not_closed_rng,  *not_closed_rng,  *not_virtual_rng, *not_virtual_rng});
  
  shared_ptr<Tensor_<double>> XTens_data = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(X_ranges , 1.0 );
  Data_map->emplace("X" , XTens_data);
  shared_ptr<Tensor_<double>> YTens_data = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(Y_ranges , 1.0 );
  Data_map->emplace("Y" , YTens_data);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Build_Compute_Lists() { 
///////////////////////////////////////////////////////////////////////////////////////////////////////

//  shared_ptr<vector<string>> HT = make_shared<vector<string>>(vector<string> { "H" , "T" });
//  Expr_Info->Set_BraKet_Ops( HT, "HT" ) ;

//  shared_ptr<vector<string>> LT = make_shared<vector<string>>(vector<string> { "L" , "T" });
//  Expr_Info->Set_BraKet_Ops( LT, "LT" ) ;
 
//  shared_ptr<vector<string>> HX = make_shared<vector<string>>(vector<string> { "H" , "X" });
//  Expr_Info->Set_BraKet_Ops( HX, "HX" ) ;
 
//  shared_ptr<vector<string>> hX = make_shared<vector<string>>(vector<string> { "h" , "X" });
//  Expr_Info->Set_BraKet_Ops( hX, "hX" ) ;
 
  shared_ptr<vector<string>> YH = make_shared<vector<string>>(vector<string> { "Y", "H" });
  Expr_Info->Set_BraKet_Ops( YH, "YH" ) ;                                             
                                                                                     
  shared_ptr<vector<string>> Yh = make_shared<vector<string>>(vector<string> { "Y", "h" });
  Expr_Info->Set_BraKet_Ops( Yh, "Yh" ) ;

  shared_ptr<vector<string>> HamT = make_shared<vector<string>>(vector<string> { "YH" } ); 

  Expr_Info->Build_Expression(HamT, "test_case") ;

  return ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(string Expression_name) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(string expression_name ) " << endl;

  if ( scalar_results_map->find(Expression_name) != scalar_results_map->end() )  
    cout << "WARNING : You have already calculated this expression....." << Expression_name
    << " = " << scalar_results_map->at(Expression_name) << endl;
 
  shared_ptr<Equation<double>> Expr = Expr_Info->expression_map->at(Expression_name); 
  double result = 0.0;

  shared_ptr<Equation_Computer::Equation_Computer> Expr_computer = make_shared<Equation_Computer::Equation_Computer>(ref, Expr, range_conversion_map, Data_map);

  Gamma_Machine = make_shared<Gamma_Computer::Gamma_Computer>( Expr->GammaMap, CIvec_data_map, Sigma_data_map, Gamma_data_map, Determinants_map, range_conversion_map );
  
  Gamma_Machine->convert_civec_to_tensor( cc_->data(0), 0 );


  //Loop through gamma names in map, ultimately the order should be defined so as to be maximally efficient, but leave this for now.
  for ( auto AG_contrib : *(Expr->GammaMap) ) {
   
    string Gamma_name = AG_contrib.first;

    // Build A_tensor to hold sums of different A-tensors
    shared_ptr<Tensor_<double>> A_combined_data = make_shared<Tensor_<double>>( *(Expr_computer->Get_Bagel_IndexRanges(Expr->GammaMap->at(Gamma_name)->id_ranges)) );
    A_combined_data->allocate();
    A_combined_data->zero();
    cout << " Gamma_name  = " << Gamma_name << endl;

    // Loop through A-tensors needed for this gamma
    auto  A_contrib_loc =  Expr->G_to_A_map->find(Gamma_name) ;
  
    if ( A_contrib_loc !=  Expr->G_to_A_map->end() ) { 
   
      for ( auto A_contrib : *(Expr->G_to_A_map->at(Gamma_name))){
    
        if (check_AContrib_factors(A_contrib.second))
          continue;
      
         shared_ptr<vector<shared_ptr<Tensor_<double>>>> gamma_tensors = Expr_computer->get_gamma( Gamma_name );
         cout << Gamma_name << "   2 " << endl;
      
        //fudge, purging of A_contribs should happen in gamma_generator or Equation
        bool skip = false ;
             cout << Gamma_name << " 3 " << endl;
	print_AContraction_list(Expr->ACompute_map->at(A_contrib.first), A_contrib.first);
        cout << Gamma_name << " 4 " << endl;
        Expr_computer->Calculate_CTP(A_contrib.first);
        cout << Gamma_name << " 5 " << endl;
        for ( int qq = 0 ; qq != A_contrib.second.id_orders.size(); qq++){ 
          shared_ptr<Tensor_<double>> A_contrib_reordered = Expr_computer->reorder_block_Tensor( A_contrib.first, make_shared<vector<int>>(A_contrib.second.id_order(qq)) );
          A_combined_data->ax_plus_y( (double)(A_contrib.second.factor(qq).first), Data_map->at(A_contrib.first));
        }
        cout << "added " << A_contrib.first << endl; 
        cout << "=========================================================================================================" << endl << endl;
      }
      
      
      if ( Gamma_name != "ID" ) {
        //  double bob = A_combined_data->dot_product(gamma_tensors->at(ii));
        double bob = 1;
        result += bob; // A_combined_data->dot_product(gamma_tensors->at(ii));  
      } else { 
        result += A_combined_data->rms(); //really dumb, and wrong for negative values...
      }
    }
  }   
  cout << Expression_name << " = " << result << endl;
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

 // shared_ptr<vector<string>> free_ranges = make_shared<vector<string>>(vector<string> {"free", "free", "free", "free"});
 // shared_ptr<vector<string>> omega_ranges = make_shared<vector<string>>(vector<string> {"notcor", "notcor", "notvir", "notvir"});
 // shared_ptr<vector<string>> omega_ranges_dag = make_shared<vector<string>>(vector<string> {"notvir", "notvir", "notcor", "notcor" });
 //
 // shared_ptr<Tensor_<double>> All_ones_tens_free = Eqn_computer->get_uniform_Tensor(free_ranges, 1.0 );
 // shared_ptr<Tensor_<double>> All_twos_tens_omega = Eqn_computer->get_uniform_Tensor(omega_ranges, 2.0 );
 // shared_ptr<Tensor_<double>> All_twos_tens_omega_dag = Eqn_computer->get_uniform_Tensor(omega_ranges_dag, 2.0 );
 // 
 // Data_map->at("X") =  All_ones_tens_free ;
 // Data_map->at("T") =  All_twos_tens_omega ;
 // Data_map->at("L") =  All_twos_tens_omega_dag ;


#endif
