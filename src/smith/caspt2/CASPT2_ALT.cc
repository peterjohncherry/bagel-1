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

  CTP_map = make_shared<map<string, shared_ptr<CtrTensorPart<Tensor_<double>>>>>();
  CTP_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  gamma_data_map = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
}

/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::test() { 
/////////////////////////////////////////////////////////////////////////////////
  auto Eqn = make_shared<Equation<Tensor>>();

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
  X_data->allocate();

  auto XTens = Eqn->Build_TensOp("X", X_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
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
  T_data->allocate();

  auto TTens = Eqn->Build_TensOp("T", T_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
  CTP_data_map->emplace("T", T_data );
  Eqn->T_map->emplace("T", TTens);
  ////////////////////////////////////// rdms (no derivs) ////////////////////////////////////////////////////////////////
  auto BraKet_Tensors1 = make_shared<vector< shared_ptr<TensOp<Tensor_<double>>> > >( vector<shared_ptr<TensOp<Tensor_<double>>>> { XTens,  TTens} );
  auto BraKet_Tensors2 = make_shared<vector< shared_ptr<TensOp<Tensor_<double>>> > >( vector<shared_ptr<TensOp<Tensor_<double>>>> { XTens,  TTens} );

  auto BraKet_List = make_shared<std::vector<std::shared_ptr<vector< shared_ptr<TensOp<Tensor_<double>>> > >>>();
                         
  BraKet_List->push_back(BraKet_Tensors1);
  BraKet_List->push_back(BraKet_Tensors2);

  Eqn->equation_build(BraKet_List);
  CTP_map = Eqn->CTP_map;

  auto Eqn_computer = make_shared<Equation_Computer::Equation_Computer>(ref, Eqn, CTP_data_map, range_conversion_map );

  for (auto MM = 0 ; MM != nstate_ ; MM++){
    for (auto NN = 0 ; NN != nstate_ ; NN++){
       
      for ( auto gidxs_loc = Eqn->CMTP_Eqn_Compute_List->rbegin(); gidxs_loc != Eqn->CMTP_Eqn_Compute_List->rend(); ++gidxs_loc ){
          
        auto gamma_range = make_shared<vector<string>>(gidxs_loc->first);
        auto gamma_ranges = make_shared<vector<shared_ptr<vector<string>>>>(1, gamma_range);
        cout << "gamma_range->size() = " << gamma_range->size();
        
        for(int ii = 2 ; ii!= gamma_range->size(); ii+=2)
          gamma_ranges->push_back(make_shared<vector<string>>(gamma_range->begin(), gamma_range->end()-ii));  
        
        auto gamma_tensors =  Eqn_computer->get_gammas( MM, NN, gamma_ranges ) ;
        //  CTP_data_map->emplace("T", T2_all[MM]->at(NN) );
        cout << "got gammas" <<endl; 
        for ( auto A_contribs : *(Eqn->CMTP_Eqn_Compute_List->at(*gamma_range))){
          pair<int,int> ctr_factor = A_contribs.second;
          for (auto A_contrib : *A_contribs.first){
            for (auto ctr_op : *(Eqn->ACompute_map->at(A_contrib))){

              if ( CTP_data_map->find(get<3>(ctr_op)) == CTP_data_map->end() ){
                if ( get<0> (ctr_op) == get<3>(ctr_op)){ 
                  if ( CTP_map->at(get<0>(ctr_op))->id_ranges->size()  > Eqn->T_map->at(get<0>(ctr_op).substr(0,1))->idxs->size()){
                     cout << "uncontracted multi tensor, do not store" << endl;       
                  } else {
                    cout << get<0> (ctr_op)<< " -->  no contraction" << endl; 
                    shared_ptr<Tensor_<double>>  New_Tdata  =  Eqn_computer->get_block_Tensor(get<0>(ctr_op));
                    CTP_data_map->emplace(get<0>(ctr_op), New_Tdata); 
                  }
                
                } else if ( get<0> (ctr_op) != get<1>(ctr_op)){
                  cout << get<0> (ctr_op)<<  "!=" <<  get<1>(ctr_op)  << " --> contract different tensors" << endl; 
                  shared_ptr<Tensor_<double>>  New_Tdata ; // =  Eqn_computer->contract_same_tensor( get<2>(ctr_op), get<0>(ctr_op), get<1>(ctr_op));
                  CTP_data_map->emplace(get<3>(ctr_op), New_Tdata); 
                
                } else {
                  cout << get<0> (ctr_op)<<  "==" <<  get<1>(ctr_op)  << " --> contract on same tensor" <<  endl; 
                  shared_ptr<Tensor_<double>>  New_Tdata  =  Eqn_computer->contract_on_same_tensor( get<2>(ctr_op), get<0>(ctr_op)); 
                  CTP_data_map->emplace(get<3>(ctr_op), New_Tdata); 
                }
              } else {
                shared_ptr<Tensor_<double>> Merged_Aterms ; // add new A term; 
              } 
            }
          }
        }   
      }
    }
  }
  
  return;
}


//  CTP_data_map->emplace(get<3>(ctr_op),  Eqn_computer->contract_same_tensor( get<2>(ctr_op), Eqn->CTP_map->at(get<0>(ctr_op)), CTP_data_map->at(get<2>(ctr_op))); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//Builds data map for use by contraction routines
///////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> CASPT2_ALT::CASPT2_ALT::get_Tensor_data( string Tname ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(CTP_map->at(Tname)->id_ranges);
  
   shared_ptr<Tensor_<double>> Tblock = make_shared<Tensor_<double>>(*Bagel_id_ranges);  
 
   //unique_ptr<double[]> Tblock_data = CTP_data_map->at(to_string(Tname[0]))->get_block(Bagel_id_ranges);

   shared_ptr<Tensor_<double>> new_T = make_shared<Tensor_<double>>();

   return new_T;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>> CASPT2_ALT::CASPT2_ALT::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map->at(rng)));

  return ranges_Bagel;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> CASPT2_ALT::CASPT2_ALT::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}


#endif
