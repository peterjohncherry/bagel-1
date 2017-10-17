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
//  closed_rng  =  make_shared<IndexRange>(IndexRange(ref->nclosed()-ref->ncore(), max, 0, ref->ncore()));
//  active_rng  =  make_shared<IndexRange>(IndexRange(ref->nact(), min(10,max), closed_rng->nblock(), ref->ncore()+closed_rng->size()));
//  virtual_rng =  make_shared<IndexRange>(IndexRange(ref->nvirt(), max, closed_rng->nblock()+active_rng->nblock(), ref->ncore()+closed_rng->size()+active_rng->size()));
//  free_rng = make_shared<IndexRange>(*closed_rng); free_rng->merge(*active_rng); free_rng->merge(*virtual_rng);

  IndexRange closed  = orig_cpt2_in.closed_; 
  IndexRange active  = orig_cpt2_in.active_;  
  IndexRange virt    = orig_cpt2_in.virt_;  
  IndexRange free    = orig_cpt2_in.all_; 

  closed_rng  = make_shared<IndexRange>(closed); 
  active_rng  = make_shared<IndexRange>(active);  
  virtual_rng = make_shared<IndexRange>(virt);  
  free_rng    = make_shared<IndexRange>(free); 

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

  CTP_data_map->emplace("T" , T2_all[0]->at(0));
  CTP_data_map->emplace("L" , T2_all[0]->at(0));
  CTP_data_map->emplace("H" , H_2el_all);
  CTP_data_map->emplace("f" , H_1el_all);

  Expr_info = make_shared<Expression_Info<double>>(norb_, nelea_, neleb_, true);
  Expr_info->Construct_Tensor_Ops();
  
}

////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::solve() {
////////////////////////////////////////////////////////////////////
cout <<  " CASPT2_ALT::CASPT2_ALT::solve() " << endl;
  
  Build_Compute_Lists();

  Execute_Compute_List("test_case");

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  
  // Get < M | W H | N  >

  // Input into linear equation solver to get T amplitudes
  
  // Should get \sum < M | W ( f- E_{L}+E_{s} ) T_{LN} | N>  in solver
 
  // Construct effective Hamiltonian
  
  // symmetrize and diagonalize effective Hmailtonian

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Build_Compute_Lists() { 
///////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<string>> HT = make_shared<vector<string>>(vector<string> { "H" , "T" });
  Expr_info->Set_BraKet_Ops( HT, "HT" ) ;

  shared_ptr<vector<string>> LT = make_shared<vector<string>>(vector<string> { "L" , "T" });
  Expr_info->Set_BraKet_Ops( LT, "LT" ) ;
 
  shared_ptr<vector<string>> HamT = make_shared<vector<string>>(vector<string> { "HT", "LT" } ); 

  Expr_info->Build_Expression(HamT, "test_case") ;

  return ;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(string Expression_name) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "CASPT2_ALT::CASPT2_ALT::Execute_Compute_List(shared_ptr<vector<string>> BraKet_names, string expression_name ) " << endl;

  shared_ptr<Equation<double>> Expr = Expr_info->expression_map->at(Expression_name); 

  shared_ptr<Equation_Computer::Equation_Computer> Expr_computer = make_shared<Equation_Computer::Equation_Computer>(ref, Expr, CTP_data_map, range_conversion_map );

  //Get Amap for each gamma
  vector<string> Gname_vec(Expr->G_to_A_map->size());
  {
    compare_string_length csl;
    int ii = 0 ; 
    for ( auto G_to_A_map_it : *(Expr->G_to_A_map) ){
      Gname_vec[ii] = G_to_A_map_it.first;
      ii++;
    }
  std::sort(Gname_vec.begin(), Gname_vec.end(), csl); 
  }

  shared_ptr<vector<shared_ptr<Tensor_<double>>>> gamma_tensors = Expr_computer->get_gammas( 0, 0, Gname_vec[0] );
  cout << "returned from get_gammas in tact " << endl;

  for ( string Gamma_name : Gname_vec ) {
    // Build A_tensor to hold sums of different A-tensors
    shared_ptr<Tensor_<double>> A_combined_data = make_shared<Tensor_<double>>( *(Expr_computer->Get_Bagel_IndexRanges(Expr->GammaMap->at(Gamma_name)->id_ranges)) );
    A_combined_data->allocate();

    // Loop through A-tensors needed for this gamma,
    for ( auto A_contrib : *(Expr->G_to_A_map->at(Gamma_name))){

      //fudge, should be replaced.
      pair<int,int> A_factor = A_contrib.second;
      pair<int,int> zero_pair = make_pair(0,0);
      if (A_factor == zero_pair ) {
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
      A_combined_data->ax_plus_y( (double)(A_factor.first), CTP_data_map->at(A_contrib.first));
      cout << "added " << A_contrib.first << endl; 
      cout << "=========================================================================================================" << endl << endl;
    }
  }   

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
