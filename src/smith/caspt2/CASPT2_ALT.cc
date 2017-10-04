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
//using namespace bagel::SMITH::CASPT2_ALT_EQN_INFO;

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
 // all_gamma1 = make_shared<VecRDM<1>>();
 // all_gamma2 = make_shared<VecRDM<2>>();
 // all_gamma3 = make_shared<VecRDM<3>>();
 
  T2_all     = orig_cpt2->t2all_;
  lambda_all = orig_cpt2->lall_;
  H_1el_all  = orig_cpt2->f1_;
  H_2el_all  = orig_cpt2->H_2el_;
 
  const int max = ref->maxtile();
  closed_rng  =  make_shared<IndexRange>(IndexRange(ref->nclosed()-ref->ncore(), max, 0, ref->ncore()));
  active_rng  =  make_shared<IndexRange>(IndexRange(ref->nact(), min(10,max), closed_rng->nblock(), ref->ncore()+closed_rng->size()));
  virtual_rng =  make_shared<IndexRange>(IndexRange(ref->nvirt(), max, closed_rng->nblock()+active_rng->nblock(), ref->ncore()+closed_rng->size()+active_rng->size()));
  free_rng= closed_rng; free_rng->merge(*active_rng); free_rng->merge(*virtual_rng);

  range_conversion_map = make_shared<map<string, shared_ptr<IndexRange>>>();
  range_conversion_map->emplace("cor", closed_rng);//change the naming of the ranges from cor to clo... 
  range_conversion_map->emplace("act", active_rng);
  range_conversion_map->emplace("A", active_rng); //fudge, but leave for now
  range_conversion_map->emplace("vir", virtual_rng);
  range_conversion_map->emplace("free", free_rng);

  CTP_map = make_shared<map<string, shared_ptr<CtrTensorPart<double>>>>();
  CTP_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  gamma_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
}

/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::test() { 
/////////////////////////////////////////////////////////////////////////////////
  auto Eqn = make_shared<Equation<double>>();

  Eqn->Initialize();
  
  //////////////////Initialization section/////////////////////////
  ///////////This should be read in from an input file/////////////
  //spinfree orbital ranges
  vector<string> free     = {"cor", "act", "vir"};
  vector<string> not_core = {"act", "vir"};
  vector<string> not_act  = {"cor", "vir"};
  vector<string> not_virt = {"cor", "act"};
  vector<string> core = {"cor"};
  vector<string> act  = {"act"};
  vector<string> virt = {"vir"};
 
  ///////////////////////////////////X Tensor ////////////////////////////////////////
  string X_TimeSymm = "none";
  auto X_factor = make_pair(1.0,1.0);
  auto X_idxs = make_shared<vector<string>>(vector<string> {"X0", "X1", "X2", "X3"});
  auto X_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
  auto X_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free,free,free,free }); 
  auto X_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)>  X_constraints = { &always_true };
  vector<IndexRange> TEMP_X_ranges = {*free_rng, *free_rng, *free_rng, *free_rng };
  shared_ptr<Tensor_<double>> X_data = make_shared<Tensor_<double>>( TEMP_X_ranges);  
  shared_ptr<double> X_dummy_data;
  X_data->allocate();

  auto XTens = Eqn->Build_TensOp("X", X_dummy_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
  CTP_data_map->emplace("X", X_data );
  Eqn->T_map->emplace("X", XTens);
  ///////////////////////////////////////////////////// T Tensor /////////////////////////////////////////////////////////////////
  string T_TimeSymm = "none";
  auto T_factor = make_pair(1.0,1.0);
  auto T_idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
  auto T_aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} ); 
  auto T_idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });   
  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> T_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)> T_constraints = { &NotAllAct };
  vector<IndexRange> TEMP_T_ranges = {*free_rng, *free_rng, *free_rng, *free_rng };
  shared_ptr<Tensor_<double>> T_data = make_shared<Tensor_<double>>( TEMP_T_ranges );  
  shared_ptr<double> T_dummy_data;
  T_data->allocate();

  auto TTens = Eqn->Build_TensOp("T", T_dummy_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
  CTP_data_map->emplace("T", T_data );
  Eqn->T_map->emplace("T", TTens);
  ////////////////////////////////////// rdms (no derivs) ////////////////////////////////////////////////////////////////
  auto BraKet_Tensors1 = make_shared<vector< shared_ptr<TensOp<double>> > >( vector<shared_ptr<TensOp<double>>> { XTens,  TTens} );
  auto BraKet_Tensors2 = make_shared<vector< shared_ptr<TensOp<double>> > >( vector<shared_ptr<TensOp<double>>> { XTens,  TTens} );

  auto BraKet_List = make_shared<std::vector<std::shared_ptr<vector< shared_ptr<TensOp<double>> > >>>();
                         
  BraKet_List->push_back(BraKet_Tensors1);

  Eqn->equation_build(BraKet_List);
  CTP_map = Eqn->CTP_map;
  auto Eqn_computer = make_shared<Equation_Computer::Equation_Computer>(ref, Eqn, CTP_data_map, range_conversion_map );

  //Get Amap for each gamma
  vector<string> Gname_vec(Eqn->G_to_A_map->size());
  {
  compare_string_length csl;
  int ii = 0 ; 
  for ( auto G_to_A_map_it : *(Eqn->G_to_A_map) ){
    Gname_vec[ii] = G_to_A_map_it.first;
    ii++;
  }
  std::sort(Gname_vec.begin(), Gname_vec.end(), csl); 
  cout << "sorted gamma_names = [ "; cout.flush();
  for (auto gname_ : Gname_vec ) {cout << gname_ << " " ; cout.flush(); } cout << "]" << endl;  
  }
  for ( string Gamma_name : Gname_vec ) {
    auto gamma_tensors = Eqn_computer->get_gammas( 0, 0, Gamma_name );

    // Build A_tensor to hold sums of different A-tensors
    shared_ptr<Tensor_<double>> A_combined_data = make_shared<Tensor_<double>>( *(Eqn_computer->Get_Bagel_IndexRanges(Eqn->GammaMap->at(Gamma_name)->id_ranges)) );
    A_combined_data->allocate();

    //auto A_contrib_data = make_shared<Tensor_<double>();
    // Loop through A-tensors needed for this gamma,
    for ( auto A_contrib : *(Eqn->G_to_A_map->at(Gamma_name))){

      pair<int,int> A_factor = A_contrib.second;
      
      cout << "=========================================================================================================" << endl;
      cout << A_contrib.first << endl;
      cout << "=========================================================================================================" << endl;
      for (shared_ptr<CtrOp_base> ctr_op : *(Eqn->ACompute_map->at(A_contrib.first))){
        if ( ctr_op->ctr_type()[0] == 'd' ){
          cout << "[" << ctr_op->T1name() << " , " << ctr_op->T2name() << " , (";
          cout << ctr_op->T1_ctr_abs_pos() << "," <<  ctr_op->T2_ctr_abs_pos() << ")" << " , " << ctr_op->Tout_name() << " ] " ;
          cout << ctr_op->ctr_type() << endl;
        } else if (ctr_op->ctr_type()[0] == 's' ){
          cout << ctr_op->Tout_name() << endl;
          cout << "[" << ctr_op->T1name() << " , " << ctr_op->T1name() << " , (";
          cout << ctr_op->ctr_abs_pos().first << "," <<  ctr_op->ctr_abs_pos().second << ")" << " , " << ctr_op->Tout_name() << " ] " ;
          cout << ctr_op->ctr_type() << endl;
        }
      }

      cout << "=========================================================================================================" << endl;

      // Loop through compute list for this A-Tensor
      for (shared_ptr<CtrOp_base> ctr_op : *(Eqn->ACompute_map->at(A_contrib.first))){

        string CTP1_name = ctr_op->T1name();
        string CTPout_name = ctr_op->Tout_name();

        string CTP2_name; 
        pair<int,int> ctr_todo;
        cout << "CTP1_name = " << CTP1_name << "    last two =  " << CTP1_name.substr( CTP1_name.length() - 2 ) << endl;
        if (ctr_op->ctr_type()[0] == 's'){
          CTP2_name = ctr_op->T1name();
          ctr_todo = ctr_op->ctr_rel_pos();
        } else if (ctr_op->ctr_type()[0] == 'd'){
          CTP2_name = ctr_op->T2name();
          ctr_todo = make_pair(ctr_op->T1_ctr_rel_pos(), ctr_op->T2_ctr_rel_pos());
        }  

        // check if this is an uncontracted multitensor (0,0) && check if the data is in the map
        if( CTP_data_map->find(CTPout_name) == CTP_data_map->end() ) {
      
          if ( CTP1_name == CTPout_name){  cout << " : no contraction, fetch this tensor part" << endl; 
            shared_ptr<Tensor_<double>>  New_Tdata  =  Eqn_computer->get_block_Tensor(CTP1_name);
            CTP_data_map->emplace(CTP1_name, New_Tdata); 
      
          } else if ( ctr_op->ctr_type()[0] == 'd' ){ cout << " : contract different tensors" << endl; 
            shared_ptr<Tensor_<double>>  New_Tdata  =  Eqn_computer->contract_different_tensors( ctr_todo, CTP1_name, CTP2_name);

            CTP_data_map->emplace(CTPout_name, New_Tdata); 
          
          } else if ( ctr_op->ctr_type()[0] == 's' ) { cout << " : contract on same tensor" <<  endl; 
            shared_ptr<Tensor_<double>>  New_Tdata  =  Eqn_computer->contract_on_same_tensor( ctr_todo, CTP1_name); 
            CTP_data_map->emplace(CTPout_name, New_Tdata); 
          }
        } else { 
          cout << " did not satisfy any if .... " <<  endl;
        }
        cout << "A_contrib.first = " << A_contrib.first << endl;
        cout << "CTPout_name =  " << CTPout_name << endl;
      }
      cout << "added " << A_contrib.first << endl; 
      cout << "=========================================================================================================" << endl << endl;
      // add this contribution to the merged A_tensor
    }
  }   

  return;
}


#endif
