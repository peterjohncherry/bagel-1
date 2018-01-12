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
cout << "PropTool::PropTool::PropTool" << endl;
  maxtile_   = idata->get<int>("maxtile", 10);
  cimaxtile_ = idata->get<int>("cimaxtile", (ciwfn_->civectors()->size() > 10000) ? 100 : 10);
  
 // op_name_list_ = idata->get<int>("op_info") ;
  

  //should read from idata
  nclosed_  = idata->get<int>("nclosed", ref_->nclosed()) ; cout << " nclosed_ = " <<  nclosed_  << endl;
  nact_     = idata->get<int>("nact", ref_->nact());        cout << " nact_    = " <<  nact_  << endl;
  ncore_    = idata->get<int>("ncore", ciwfn_->ncore());    cout << " ncore_   = " <<  ncore_  << endl;
  nvirt_    = idata->get<int>("nvirt", ref_->nvirt());      cout << " nvirt_   = " <<  nvirt_  << endl;
  nocc_     = nclosed_ + nact_; 

  vector< Term_Init > expression_init;
 
  auto expression_inp = idata->get_child("expression"); 
  for ( auto& bk_info_tree : *expression_inp ){

    Term_Init term_inp;

    for ( auto& bk_info : *bk_info_tree ) {
      auto op_list =  bk_info->get_child("ops"); 

      vector<string> op_names(0);  
      for (auto& op_name : *op_list)
        op_names.push_back(lexical_cast<string>(op_name->data()));

      term_inp.op_names.push_back(op_names);

      auto target_states_ptree =  bk_info->get_child("target_states"); 
      vector<int> target_states(0);
      for (auto& state_num : *target_states_ptree)
        target_states.push_back(lexical_cast<int>(state_num->data()));

      term_inp.target_states.push_back(target_states);
      term_inp.types.push_back(bk_info->get<string>("type"));
      term_inp.factors.push_back(bk_info->get<double>("factor"));
       
    }
    expression_init.push_back(term_inp);
  }

  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  
  range_conversion_map_ = make_shared<map<string, shared_ptr<SMITH::IndexRange>>>();

  civectors_ = ciwfn_->civectors();
  set_target_info() ;
  set_range_info();
   
  sys_info_ = make_shared<System_Info<double>>(targets_info_, true);
  expression_map_ = sys_info_->expression_map;
  expression_machine_ = make_shared<SMITH::Expression_Computer::Expression_Computer<double>>( civectors_, expression_map_, range_conversion_map_, tensop_data_map_, 
                                                                                              gamma_data_map_, sigma_data_map_, civec_data_map_ );

  

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_range_info" << endl;

  closed_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nclosed_-ncore_, maxtile_, 0, ncore_));
  active_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nact_, min((size_t)10,maxtile_), closed_rng_->nblock(), ncore_ + closed_rng_->size()));
  virtual_rng_ =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nvirt_, maxtile_, closed_rng_->nblock()+ active_rng_->nblock(), ncore_+closed_rng_->size()+active_rng_->size()));
  free_rng_    = make_shared<SMITH::IndexRange>(*closed_rng_);
  free_rng_->merge(*active_rng_);
  free_rng_->merge(*virtual_rng_);

  cout << " set first ranges " << endl;
  not_closed_rng_  =  make_shared<SMITH::IndexRange>(*active_rng_); not_closed_rng_->merge(*virtual_rng_);
  not_active_rng_  =  make_shared<SMITH::IndexRange>(*closed_rng_); not_active_rng_->merge(*virtual_rng_);
  not_virtual_rng_ =  make_shared<SMITH::IndexRange>(*closed_rng_); not_virtual_rng_->merge(*active_rng_);

  cout << " set seocnd  ranges " << endl;
  range_conversion_map_->emplace("cor", closed_rng_); 
  range_conversion_map_->emplace("act", active_rng_);
  range_conversion_map_->emplace("vir", virtual_rng_);
  range_conversion_map_->emplace("free", free_rng_);
 
  cout << " set third ranges " << endl;
  range_conversion_map_->emplace("notcor", not_closed_rng_);
  range_conversion_map_->emplace("notact", not_active_rng_);
  range_conversion_map_->emplace("notvir", not_virtual_rng_); 
                      
  cout << " set ci ranges " << endl;
  for ( int ii : target_states_ ) {

    cout << " ii = " << ii << endl; 
    range_conversion_map_->emplace( get_civec_name( ii , civectors_->data(ii)->det()->norb(), civectors_->data(ii)->det()->nelea(), civectors_->data(ii)->det()->neleb()),
                                                   make_shared<SMITH::IndexRange>(civectors_->data(ii)->det()->size(), cimaxtile_ ));  
    
    shared_ptr<SMITH::IndexRange>  ci_index_ranges =  make_shared<SMITH::IndexRange>(civectors_->data(ii)->det()->size(), cimaxtile_ );
    cout << "cirngs = [ ";  for (auto irng : ci_index_ranges->range()) { cout << irng.size()  << " "; }; cout << "] " << endl;
 
  }    
  cout << " set done with rangesd  ranges " << endl;

 return;
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_target_info() {
//////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_target_info" << endl;
  targets_info_ = make_shared<StatesInfo<double>> ( target_states_ ) ;
  
  for ( int state_num : target_states_ ) 
     targets_info_->add_state( civectors_->data(state_num)->det()->nelea(), civectors_->data(state_num)->det()->neleb(),
                               civectors_->data(state_num)->det()->norb(), state_num ) ;
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::build_expressions( vector<int>& target_states, vector<pair<vector<string>,double>>& BK_info_list, string term_type ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::calculate_term" << endl;

  map< pair<string,string>, shared_ptr<vector<Term_Info<double>>> > Term_info_map_;
  int num_states = target_states.size();
  vector<string> op_names_got(0);
  map<pair<string,string>,shared_ptr<vector<Term_Info<double>>> > term_map; 

  vector<string>  expression_list(0);
  // Building all necessary expressions 
  for ( int ii : target_states ) {
    string Bra_name = targets_info_->name(ii); // TODO change to vector for relativistic case           
    for ( int  jj : target_states ) {
      string Ket_name = targets_info_->name(jj); // TODO change to vector for relativistic case           

      shared_ptr<vector<Term_Info<double>>> term_info_ss = make_shared<vector<Term_Info<double>>>(); 
      for ( pair<vector<string>,double> BK_info : BK_info_list )
        term_info_ss->push_back(Term_Info<double>( BK_info, Bra_name, Ket_name, term_type ));
    
      term_map.emplace(make_pair(Bra_name, Ket_name) , term_info_ss ); 
      expression_list.push_back(sys_info_->Build_Expression( *term_info_ss ));
    }
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::build_op_tensors( vector<string>& expression_list ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::build_op_tensors" << endl;

  // Creating tensors from existing matrices; seperate loop as must run through all states first to make proper use of symmetry
  for (string expression_name : expression_list ) {
    for ( auto tensop_it : *(sys_info_->expression_map->at(expression_name)->T_map) ) {
      std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks = *(tensop_it.second->unique_range_blocks());
      for ( shared_ptr<const vector<string>> range_block : unique_range_blocks ) {
        shared_ptr<vector<SMITH::IndexRange>> range_block_bgl = convert_to_indexrange( range_block ); 
      }
    }
  }


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<SMITH::IndexRange>> PropTool::PropTool::convert_to_indexrange( shared_ptr<const vector<string>> range_block_str ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::convert_to_indexrange" << endl;

  shared_ptr<vector<SMITH::IndexRange>> range_block = make_shared<vector<SMITH::IndexRange>>( range_block_str->size() );
  for( int ii = 0; ii != range_block_str->size() ; ii++ ) 
    range_block->at(ii) = *(range_conversion_map_->at( range_block_str->at(ii) ));

  return range_block; 

}
