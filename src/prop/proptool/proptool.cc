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
#include <src/prop/proptool/moint_computer.h>

using namespace std;
using namespace bagel;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PropTool::PropTool::PropTool(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r): 
                   idata_(idata), geom_(g), ref_(r), ciwfn_(ref_->ciwfn()), civectors_(ciwfn_->civectors())  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::PropTool" << endl;

  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();

  range_conversion_map_ = make_shared<map<string, shared_ptr<SMITH::IndexRange>>>();

  // get user specified variables (e.g. ranges, constant factors) which may appear in term definitions
  get_expression_variables( idata->get_child("variables") );

  //Initializing range sizes either from idate or reference wfn 
  maxtile_   = idata->get<int>("maxtile", 10);
  cimaxtile_ = idata->get<int>("cimaxtile", (ciwfn_->civectors()->size() > 10000) ? 100 : 10);
  nclosed_  = idata->get<int>( "nclosed" , ref_->nclosed()); cout << " nclosed_ = " <<  nclosed_  << endl;
  nact_     = idata->get<int>( "nact"  , ref_->nact());      cout << " nact_    = " <<  nact_  << endl;
  ncore_    = idata->get<int>( "ncore" , ciwfn_->ncore());   cout << " ncore_   = " <<  ncore_  << endl;
  nvirt_    = idata->get<int>( "nvirt" , ref_->nvirt());     cout << " nvirt_   = " <<  nvirt_  << endl;
  nocc_     = nclosed_ + nact_; 
  nfrozenvirt_ = idata->get<int>( "nfrozenvirt", 0 );

  // leave for now
  block_diag_fock_ = false;
  gaunt_    = false;
  breit_    = false;
  set_ao_range_info();

  //Getting info about target expression (this includes which states are relevant)
  shared_ptr<vector<Term_Init<double>>> expression_init = get_expression_init( idata->get_child("expression") ); 
  set_target_state_info();
  set_ci_range_info();

  sys_info_ = make_shared<System_Info<double>>( targets_info_, true );
  shared_ptr< const PTree > ops_def_tree = idata->get_child_optional( "operators" ) ;
  if (ops_def_tree)
    get_new_ops_init( ops_def_tree ); 

  cout << " built user defined ops " << endl;

  expression_map_ = sys_info_->expression_map;
  expression_machine_ = make_shared<SMITH::Expression_Computer::Expression_Computer<double>>( civectors_, expression_map_, range_conversion_map_, tensop_data_map_, 
                                                                                              gamma_data_map_, sigma_data_map_, civec_data_map_ );

  cout << "expression_init->size()  = " << expression_init->size() << endl;
  for ( Term_Init<double> term_inp : *expression_init ) {
    shared_ptr<vector<string>> expr_list = build_expressions( term_inp );
    cout << "expr_list = [ " ; cout.flush();
    for ( string expr_name : *expr_list )
      cout << expr_name << " " ; cout.flush();
    cout << "]" <<  endl;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_expression_variables( shared_ptr<const PTree> variable_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto range_info_tree =  variable_def_tree->get_child_optional("ranges"); // can be mo or state.. 
  for ( auto& range_info : *range_info_tree ) {
   
    string range_name = range_info->get<string>( "name" ) ;
    auto range_vec_inp = range_info->get_child( "range_vec" ) ;
    vector<int> range_vec(0);
   
    for (auto& id : *range_vec_inp) 
      range_vec.push_back(lexical_cast<int>( id->data() ) );
   
    if ( inp_range_map_.find(range_name) != inp_range_map_.end() )
      throw runtime_error("Range \"" + range_name + "\" has been defined twice in input!!! ... aborting" ) ;
     
    inp_range_map_.emplace(range_name, range_vec); 

  }

  auto factor_info_tree =  variable_def_tree->get_child_optional("factors"); //
  for ( auto& factor_info : *factor_info_tree ) {
   
    string factor_name = factor_info->get<string>( "name" ) ;
    double factor_value = factor_info->get<double>( "value" ) ;
    if ( inp_factor_map_.find(factor_name) != inp_factor_map_.end() )
      throw runtime_error("Factor \"" + factor_name + "\" has been defined twice in input!!! ... aborting" ) ;
     
    inp_factor_map_.emplace(factor_name, factor_value ); 

  }

  cout << "USER DEFINED FACTORS " << endl;
  for ( auto elem : inp_factor_map_ ) 
    cout << elem.first << " = " << elem.second << endl; 

  cout << "USER DEFINED RANGES " << endl;
  for ( auto elem : inp_range_map_ ) {
    cout << elem.first << " = [ " ;
    for (int orb_num :  elem.second ) {
      cout << orb_num << " " ; cout.flush();
    }
    cout << "]" << endl;
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_new_ops_init( shared_ptr<const PTree> ops_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //this shouldn't be necessary, but the conversion doesn't seem to work otherwise....
  auto conv_to_bool = []( int aop ) { return aop == 1 ? true : false ; };

  // keep this so you don't have to get names correct, needs modification for spin case
  auto conv_to_range = []( string rng ) {
           if( rng[0] == 'c' ) return "cor";
           if( rng[0] == 'a' ) return "act";
           if( rng[0] == 'v' ) return "vir";
      };

  for ( auto& op_def_inp : *ops_def_tree ){
  
    auto idxs_ptree =  op_def_inp->get_child("idxs"); 
    vector<string> idxs(0);
    for (auto& idx : *idxs_ptree)
      idxs.push_back(lexical_cast<string>(idx->data()));
 
    auto ranges_ptree =  op_def_inp->get_child("ranges"); 
    vector<vector<string>> ranges(0);
    for (auto& idx_range : *ranges_ptree){
      vector<string> idx_range_vec(0);
      for (auto& ranges : *idx_range )
        idx_range_vec.push_back(conv_to_range(lexical_cast<string>(ranges->data())));
      ranges.push_back(idx_range_vec);
    }

    auto aops_ptree =  op_def_inp->get_child("aops"); 
    vector<bool> aops(0);
    for (auto& aop : *aops_ptree)
      aops.push_back(conv_to_bool(lexical_cast<int>(aop->data())));

    string                             op_name = op_def_inp->get<string>( "name" );
    double                             factor = op_def_inp->get<double>( "factor", 1.0 );
    shared_ptr<vector<string>>         idxs_ptr = make_shared<vector<string>>( idxs );
    shared_ptr<vector<bool>>           aops_ptr = make_shared<vector<bool>>( aops );
    shared_ptr<vector<vector<string>>> ranges_ptr = make_shared<vector<vector<string>>>( ranges );
    string                             TimeSymm = op_def_inp->get<string>( "TimeSymm", "none" );
    bool                               hconj = conv_to_bool(op_def_inp->get<int>( "HermConj", false ));
      
    vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> symmfuncs =  sys_info_->identity_only() ; // TODO define this by list of pairs of vectors
    vector<bool(*)(shared_ptr<vector<string>>)> constraints = { &System_Info<double>::System_Info::always_true };  // TODO define this by list of vectors 
    
    cout << "user defined op name : " << op_name << endl;
    shared_ptr<TensOp::TensOp<double>> new_op = sys_info_->Build_TensOp( op_name, idxs_ptr, aops_ptr, ranges_ptr, symmfuncs, constraints, factor, TimeSymm, hconj); 
    sys_info_->T_map->emplace( op_name, new_op );

  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr< vector<Term_Init<double> > > PropTool::PropTool::get_expression_init( shared_ptr<const PTree> expression_inp ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::get_expression_init" << endl;

  vector<int> ground = {0};

  vector< Term_Init<double>> expression_init; 
  for ( auto& term_info_tree : *expression_inp ){
  
    std::vector<std::vector<std::string>> term_init_op_names;
    std::vector<std::vector<int>> term_init_bra_states; 
    std::vector<std::vector<int>> term_init_ket_states; 
    std::vector<double> term_init_factors;
    std::vector<std::string> term_init_types; 
  
    for ( auto& term_info : *term_info_tree ) {
      auto op_list = term_info->get_child("ops"); 
  
      vector<string> op_names(0);  
      for (auto& op_name : *op_list)
        op_names.push_back(lexical_cast<string>(op_name->data()));
      
      vector<int> bra_states = inp_range_map_.at(term_info->get<string>("bra_states")); 
      vector<int> ket_states = inp_range_map_.at(term_info->get<string>("ket_states")); 
 
      bool sum_bra = ( term_info->get<string>( "Sum_over_bra", "false" ) == "false" ) ? false : true ;
      bool sum_ket = ( term_info->get<string>( "Sum_over_ket", "false" ) == "false" ) ? false : true ;
 
      term_init_op_names.push_back(op_names);
      term_init_bra_states.push_back(bra_states);
      term_init_ket_states.push_back(ket_states);
      term_init_types.push_back(term_info->get<string>("type"));
      term_init_factors.push_back(term_info->get<double>("factor"));
    }

    Term_Init<double> new_term = Term_Init<double>(term_init_op_names, term_init_bra_states, term_init_ket_states, term_init_factors, term_init_types);
    vector<vector<int>> tmp = { target_states_, new_term.all_states_ };
    target_states_ = new_term.merge_states( tmp );
  
    expression_init.push_back(new_term);
  }
  return make_shared<vector<Term_Init<double>>>(expression_init); 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_ao_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_ao_range_info" << endl;

  closed_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nclosed_-ncore_, maxtile_, 0, ncore_));
  active_rng_  =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nact_, min((size_t)10,maxtile_), closed_rng_->nblock(), ncore_ + closed_rng_->size()));
  virtual_rng_ =  make_shared<SMITH::IndexRange>(SMITH::IndexRange(nvirt_, maxtile_, closed_rng_->nblock()+ active_rng_->nblock(), ncore_+closed_rng_->size()+active_rng_->size()));
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
  
  return;                    
}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_ci_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_ci_range_info" << endl;

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
void PropTool::PropTool::set_target_state_info() {
//////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_target_info" << endl;
  targets_info_ = make_shared<StatesInfo<double>> ( target_states_ ) ;

  for ( int state_num : target_states_ ) 
     targets_info_->add_state( civectors_->data(state_num)->det()->nelea(), civectors_->data(state_num)->det()->neleb(),
                               civectors_->data(state_num)->det()->norb(), state_num );
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<string>> PropTool::PropTool::build_expressions( Term_Init<double>& term_inp ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::build_expressions" << endl;

  vector<string>  expression_list(0);
  // Loops through states, and if states contributes to that term, add to term info for building expression
  // restatement of Bra and Ket name is redundant, and redundancy means I've probably gone wrong somewhere....
  for ( int bra_num : term_inp.bra_states_merged_ ) {
    for ( int  ket_num : term_inp.ket_states_merged_ ) {

      shared_ptr<vector<BraKet<double>>> term_info = make_shared<vector<BraKet<double>>>();

      for ( int kk = 0; kk != term_inp.op_names_.size(); kk++ ){

        if ( std::find(term_inp.bra_states_[kk].begin(), term_inp.bra_states_[kk].end(), bra_num) != term_inp.bra_states_[kk].end() ) {

          if ( std::find(term_inp.ket_states_[kk].begin(), term_inp.ket_states_[kk].end(), ket_num) != term_inp.ket_states_[kk].end() ) {
            term_info->push_back(BraKet<double>( make_pair( term_inp.op_names_[kk], term_inp.factors_[kk] ), bra_num, ket_num, term_inp.types_[kk] ));
            expression_list.push_back( sys_info_->Build_Expression( *term_info ));
          }
        }
      }
    }
  }
  return make_shared<vector<string>>(expression_list);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::build_op_tensors( vector<string>& expression_list ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::build_op_tensors" << endl;

  // Creating tensors from existing matrices; seperate loop as must run through all states first to make proper use of symmetry
  for (string expression_name : expression_list ) {
    for ( auto tensop_it : *(sys_info_->T_map) ) {
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
