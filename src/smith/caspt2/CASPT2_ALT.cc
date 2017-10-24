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

////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in ) { 
////////////////////////////////////////////////////////////////////
  orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
  ref = orig_cpt2_in.info_;

  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  cc_ = ref->ciwfn()->civectors();
  det_ = ref->ciwfn()->civectors()->det();
 
  T2_all     = orig_cpt2->t2all_;
  lambda_all = orig_cpt2->lall_;
  H_1el_all  = orig_cpt2->f1_;
  H_2el_all  = orig_cpt2->H_2el_;
 
  const int max = ref->maxtile();
  closed_rng  =  make_shared<IndexRange>(IndexRange(ref->nclosed()-ref->ncore(), max, 0, ref->ncore()));
  active_rng  =  make_shared<IndexRange>(IndexRange(ref->nact(), min(10,max), closed_rng->nblock(), ref->ncore()+closed_rng->size()));
  virtual_rng =  make_shared<IndexRange>(IndexRange(ref->nvirt(), max, closed_rng->nblock()+active_rng->nblock(), ref->ncore()+closed_rng->size()+active_rng->size()));
  free_rng = make_shared<IndexRange>(*closed_rng); free_rng->merge(*active_rng); free_rng->merge(*virtual_rng);

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

  CTP_map = make_shared<map<string, shared_ptr<CtrTensorPart<double>>>>();
  CTP_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  gamma_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  scalar_results_map = make_shared<map<string, double>>();
  //Deriv_results_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();

  Expr_Info = make_shared<Expression_Info<double>>(norb_, nelea_, neleb_, true);
  Expr_Info_map = Expr_Info->expression_map;
  
  
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
  vector<IndexRange>                  X_bagel_ranges = {*not_closed_rng , *not_closed_rng, *not_virtual_rng, *not_virtual_rng} ;
  shared_ptr<double>                  X_dummy_data;
  string                              X_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>  X_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  X_constraints = { &Expression_Info<double>::Expression_Info::NotAllAct };

  shared_ptr<TensOp<double>> XTens = Expr_Info->Build_TensOp("X", X_dummy_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("X", XTens);

  /* ---- H Tensor ----  */
  pair<double,double>                H_factor = make_pair(1.0,1.0);
  shared_ptr<vector<string>>         H_idxs = make_shared<vector<string>>(vector<string> {"H0", "H1", "H2", "H3"});
  shared_ptr<vector<bool>>           H_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
  shared_ptr<vector<vector<string>>> H_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
  shared_ptr<double>                 H_dummy_data;
  string                             H_TimeSymm = "none";

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> H_symmfuncs = Expr_Info->set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  H_constraints = {  &Expression_Info<double>::Expression_Info::always_true };

  shared_ptr<TensOp<double>> HTens = Expr_Info->Build_TensOp("H", H_dummy_data, H_idxs, H_aops, H_idx_ranges, H_symmfuncs, H_constraints, H_factor, H_TimeSymm, false ) ;
  Expr_Info->T_map->emplace("H", HTens);

  CTP_data_map->emplace("H" , H_2el_all);

  
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
  CTP_data_map->emplace("T" , TTens_data);
 

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
  CTP_data_map->emplace("L" , LTens_data);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Build_Compute_Lists() { 
///////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<string>> HT = make_shared<vector<string>>(vector<string> { "H" , "T" });
  Expr_Info->Set_BraKet_Ops( HT, "HT" ) ;

  shared_ptr<vector<string>> LT = make_shared<vector<string>>(vector<string> { "L" , "T" });
  Expr_Info->Set_BraKet_Ops( LT, "LT" ) ;
 
  shared_ptr<vector<string>> HamT = make_shared<vector<string>>(vector<string> { "HT", "LT" } ); 

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

  shared_ptr<Equation_Computer::Equation_Computer> Expr_computer = make_shared<Equation_Computer::Equation_Computer>(ref, Expr, CTP_data_map, range_conversion_map );

  //Get Amap for each gamma
  vector<string> Gname_vec(Expr->G_to_A_map->size());
  {
    compare_string_length csl;
    int ii = 0 ; 
   
    // std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< std::unordered_map<std::string, std::vector<std::pair<std::vector<int> , std::pair<int,int>>> > >>>  G_to_A_map
    for ( auto G_to_A_map_it : *(Expr->G_to_A_map) ){
      Gname_vec[ii] = G_to_A_map_it.first;
      ii++;
    }
  std::sort(Gname_vec.begin(), Gname_vec.end(), csl); 
  }

  shared_ptr<vector<shared_ptr<Tensor_<double>>>> gamma_tensors = Expr_computer->get_gammas( 0, 0, Gname_vec[0] );
  cout << "returned from get_gammas in tact " << endl;
  
  cout << "Gamma names = [ " ; cout.flush();
  for ( int ii = 0 ; ii != Gname_vec.size(); ii++ ) {
   cout << Gname_vec[ii] << " " ; cout.flush();
  }
  cout << " ] " << endl;

  for ( int ii = 0 ; ii != Gname_vec.size(); ii++ ) {
  
    string Gamma_name = Gname_vec[ii];

    // Build A_tensor to hold sums of different A-tensors
    shared_ptr<Tensor_<double>> A_combined_data = make_shared<Tensor_<double>>( *(Expr_computer->Get_Bagel_IndexRanges(Expr->GammaMap->at(Gamma_name)->id_ranges)) );
    A_combined_data->allocate();

    // Loop through A-tensors needed for this gamma,
    for ( auto A_contrib : *(Expr->G_to_A_map->at(Gamma_name))){

      //fudge, purging of A_contribs should happen in gamma_generator or Equation
      bool skip = false ;
      for ( int qq = 0 ; qq != A_contrib.second.id_orders.size(); qq++) { 
         if ( A_contrib.second.factor(qq).first != 0 || A_contrib.second.factor(qq).second !=0) 
           break;
         if ( qq == A_contrib.second.id_orders.size()-1 ) {
           skip =true;
         }
      } 
      if (skip){
        cout << A_contrib.first << " does not contribute!!" << endl;
        continue;
      }
    
      cout << "=========================================================================================================" << endl;
      cout << A_contrib.first << endl;
      cout << "=========================================================================================================" << endl;
      for (shared_ptr<CtrOp_base> ctr_op : *(Expr->ACompute_map->at(A_contrib.first))){
        if ( ctr_op->ctr_type()[0] == 'd' ){
          cout << "[" << ctr_op->T1name() << " , " << ctr_op->T2name() << " , (";
          cout << ctr_op->T1_ctr_abs_pos() << "," <<  ctr_op->T2_ctr_abs_pos() << ")" << " , " << ctr_op->Tout_name() << " ] " ; cout << ctr_op->ctr_type() << endl;

        } else if (ctr_op->ctr_type()[0] == 's' ){
          cout << "[" << ctr_op->T1name() << " , " << ctr_op->T1name() << " , (";
          cout << ctr_op->ctr_abs_pos().first << "," <<  ctr_op->ctr_abs_pos().second << ")" << " , " << ctr_op->Tout_name() << " ] " ;   cout << ctr_op->ctr_type()  << endl;
        }
      }
      cout << "=========================================================================================================" << endl;
      Expr_computer->Calculate_CTP(A_contrib.first);
      for ( int qq = 0 ; qq != A_contrib.second.id_orders.size(); qq++){ 
        shared_ptr<Tensor_<double>> A_contrib_reordered = Expr_computer->reorder_block_Tensor( A_contrib.first, make_shared<vector<int>>(A_contrib.second.id_order(qq)) );
        A_combined_data->ax_plus_y( (double)(A_contrib.second.factor(qq).first), CTP_data_map->at(A_contrib.first));
      }
      cout << "added " << A_contrib.first << endl; 
      cout << "=========================================================================================================" << endl << endl;
    }
  
 
  if ( Gamma_name != "ID" ) {
    
    //shared_ptr<vector<int>> WickUtils::reorder_vector(vector<int>& neworder , const vector<int>& origvec ) {
    result += A_combined_data->dot_product(gamma_tensors->at(ii)); 
  } else { 
    result += A_combined_data->rms(); //really dumb, and wrong for negative results...
  }
  }   
  cout << Expression_name << " = " << result << endl;
  scalar_results_map->emplace(Expression_name, result);
  
  return;
}


 // shared_ptr<vector<string>> free_ranges = make_shared<vector<string>>(vector<string> {"free", "free", "free", "free"});
 // shared_ptr<vector<string>> omega_ranges = make_shared<vector<string>>(vector<string> {"notcor", "notcor", "notvir", "notvir"});
 // shared_ptr<vector<string>> omega_ranges_dag = make_shared<vector<string>>(vector<string> {"notvir", "notvir", "notcor", "notcor" });
 //
 // shared_ptr<Tensor_<double>> All_ones_tens_free = Eqn_computer->get_uniform_Tensor(free_ranges, 1.0 );
 // shared_ptr<Tensor_<double>> All_twos_tens_omega = Eqn_computer->get_uniform_Tensor(omega_ranges, 2.0 );
 // shared_ptr<Tensor_<double>> All_twos_tens_omega_dag = Eqn_computer->get_uniform_Tensor(omega_ranges_dag, 2.0 );
 // 
 // CTP_data_map->at("X") =  All_ones_tens_free ;
 // CTP_data_map->at("T") =  All_twos_tens_omega ;
 // CTP_data_map->at("L") =  All_twos_tens_omega_dag ;


#endif
