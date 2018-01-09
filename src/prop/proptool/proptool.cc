// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: pseudospin.cc
// Copyright (C) 2015 Toru Shiozaki
//
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

#include <src/prop/proptool/proptool.h>

using namespace std;
using namespace bagel;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PropTool::PropTool::PropTool(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r): 
                   idata_(idata), geom_(g), ref_(r), ciwfn_(ref_->ciwfn())  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //should read from idata
  
  
  nclosed_  = ref_->nclosed() ; 
  nact_     = ref_->nact(); 
  ncore_    = ciwfn_->ncore(); 
  nvirt_    = ref_->nvirt(); 
  nocc_     = nclosed_ + nact_; 

  civectors_ = ciwfn_->civectors();
  set_target_info() ;
  set_range_info();
   
  Sys_Info_ = make_shared<System_Info<double>>(TargetsInfo_, true);
  Expression_map_ = Sys_Info_->expression_map;
  Expression_Machine_ = make_shared<SMITH::Expression_Computer::Expression_Computer<double>>( civectors_, Expression_map_, range_conversion_map_, TensOp_data_map_, 
                                                                                              Gamma_data_map_, Sigma_data_map_, CIvec_data_map_ );

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int maxtile = 100000; // define this more rationally
  int cimaxtile = 10000;

  closed_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nclosed_-ncore_, maxtile, 0, ncore_));
  active_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nact_, min(10,maxtile), closed_rng_->nblock(), ncore_ + closed_rng_->size()));
  virtual_rng_ =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nvirt_, maxtile, closed_rng_->nblock()+ active_rng_->nblock(), ncore_+closed_rng_->size()+active_rng_->size()));
  free_rng_    = make_shared<SMITH::IndexRange>(*closed_rng_);
  free_rng_->merge(*active_rng_);
  free_rng_->merge(*virtual_rng_);

  not_closed_rng_  =  make_shared<SMITH::IndexRange>(*active_rng_); not_closed_rng_->merge(*virtual_rng_);
  not_active_rng_  =  make_shared<SMITH::IndexRange>(*closed_rng_); not_active_rng_->merge(*virtual_rng_);
  not_virtual_rng_ =  make_shared<SMITH::IndexRange>(*closed_rng_); not_virtual_rng_->merge(*active_rng_);

  range_conversion_map_->emplace("cor", closed_rng_); 
  range_conversion_map_->emplace("act", active_rng_);
  range_conversion_map_->emplace("vir", virtual_rng_);
  range_conversion_map_->emplace("free", free_rng_);
 
  range_conversion_map_->emplace("notcor", not_closed_rng_);
  range_conversion_map_->emplace("notact", not_active_rng_);
  range_conversion_map_->emplace("notvir", not_virtual_rng_); 
                      
  for ( int ii : target_states_ ) {
 
    range_conversion_map_->emplace( get_civec_name( ii , civectors_->data(ii)->det()->norb(), civectors_->data(ii)->det()->nelea(), civectors_->data(ii)->det()->neleb()),
                                                   make_shared<SMITH::IndexRange>(civectors_->data(ii)->det()->size(), cimaxtile ));  
    
    shared_ptr<SMITH::IndexRange>  ci_index_ranges =  make_shared<SMITH::IndexRange>(civectors_->data(ii)->det()->size(), cimaxtile );
    cout << "cirngs = [ ";  for (auto irng : ci_index_ranges->range()) { cout << irng.size()  << " "; }; cout << "] " << endl;
 
  }    

 return;
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_target_info() {
//////////////////////////////////////////////////////////////////////////////////////////////
  TargetsInfo_ = make_shared<StatesInfo<double>> ( target_states_ ) ;
  
  for ( int state_num : target_states_ ) 
     TargetsInfo_->add_state( civectors_->data(state_num)->det()->nelea(), civectors_->data(state_num)->det()->neleb(),
                             civectors_->data(state_num)->det()->norb(), state_num ) ;
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::calculate_term( vector<int>& target_states, vector<pair<vector<string>,double>>& BK_info_list, string term_type ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  map< pair<string,string>, shared_ptr<vector<Term_Info<double>>> > Term_info_map_;
  
  int num_states = target_states.size();

  // Building all necessary expressions 
  vector<vector<Term_Info<double>>> Term_info_list( target_states.size()*target_states.size() );
  for ( int ii : target_states ) {
    string Bra_name = TargetsInfo_->name(ii); // TODO change to vector for relativistic case           
    for ( int  jj : target_states ) {
      string Ket_name = TargetsInfo_->name(jj); // TODO change to vector for relativistic case           

      shared_ptr<vector<Term_Info<double>>> term_info_ss = make_shared<vector<Term_Info<double>>>(); 
      for ( pair<vector<string>,double> BK_info : BK_info_list ){
        term_info_ss->push_back(Term_Info<double>( BK_info, Bra_name, Ket_name, term_type ));
        for ( string Op_name : BK_info.first )  
          cout << "tmp" << endl;
//          Set_Tensor_Ops_Data( Op_name, Bra_name, Ket_name ); 
      }
  
      string expression_name = Sys_Info_->Build_Expression( *term_info_ss );
//      Expression_Machine->Evaluate_Expression( expression_name );
  
    }
  }

}
