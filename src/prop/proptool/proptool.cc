#include <src/prop/proptool/proptool.h>

using namespace std;
using namespace bagel;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PropTool::PropTool::PropTool(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r): 
                   idata_(idata), geom_(g), ref_(r), ciwfn_(ref_->ciwfn()), civectors_(ciwfn_->civectors())  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::PropTool" << endl;


  // sort out how to determine datatype!!
  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();

  inp_factor_map_ = make_shared<map<string, double>>();
  inp_indexed_factor_map_ = make_shared<map<string, shared_ptr<vector<double>>>>(); // TODO sort this
  inp_range_map_ = make_shared<map<string, shared_ptr<vector<int>>>>();

  range_conversion_map_ = make_shared<map<string, shared_ptr<SMITH::IndexRange>>>();

  term_init_map_ = make_shared<map<string, shared_ptr<Term_Init>>>();
  expression_init_map_ = make_shared<map<string, shared_ptr<Expression_Init>>>();
  equation_init_map_ = make_shared<map<string, shared_ptr<Equation_Init_Base>>>();

  // get user specified variables (e.g. ranges, constant factors) which may appear in term definitions
  get_expression_variables( idata->get_child("variables") );

  //Initializing range sizes either from idate or reference wfn 
  maxtile_   = idata->get<int>("maxtile", 10);
  cimaxtile_ = idata->get<int>("cimaxtile", (ciwfn_->civectors()->size() > 10000) ? 100 : 10);

  const bool frozen = idata->get<bool>("frozen", true);
  ncore_ = idata->get<int>("ncore", (frozen ? ref_->geom()->num_count_ncore_only()/2 : 0));
  if (ncore_)
    cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;
  nfrozenvirt_ = idata->get<int>("nfrozenvirt", 0);
  if (nfrozenvirt_)
    cout << "    * freezing " << nfrozenvirt_ << " orbital" << (nfrozenvirt_^1 ? "s" : "") << " (virtual)" << endl;

  nclosed_  = idata->get<int>( "nclosed" , ref_->nclosed()); cout << " nclosed_ = " <<  nclosed_  << endl;
  nact_     = idata->get<int>( "nact"  , ref_->nact());      cout << " nact_    = " <<  nact_  << endl;
  nvirt_    = idata->get<int>( "nvirt" , ref_->nvirt());     cout << " nvirt_   = " <<  nvirt_  << endl;
  nocc_     = nclosed_ + nact_; 
  nfrozenvirt_ = idata->get<int>( "nfrozenvirt", 0 );

  // leave for now
  block_diag_fock_ = false;
  gaunt_    = false;
  breit_    = false;
  set_ao_range_info();

  cout << "getting mo integrals " <<  endl; 
  shared_ptr<MOInt_Init<double>> moint_init = make_shared<MOInt_Init<double>>( geom_, dynamic_pointer_cast<const Reference>(ref_), ncore_, nfrozenvirt_, block_diag_fock_ );
  shared_ptr<MOInt_Computer<double>> moint_comp = make_shared<MOInt_Computer<double>>( moint_init, range_conversion_map_ );

  vector<string> test_ranges4 = { "notcor", "notcor", "notvir", "notvir" }; 
  vector<string> test_ranges2 = { "free", "free" }; 
  {
  shared_ptr<SMITH::Tensor_<double>> v2_  =  moint_comp->get_v2( test_ranges4 ) ;
  cout << " old_coeffs  v2_->norm() = " << v2_->norm() << endl; 
  }
  shared_ptr<SMITH::Tensor_<double>> h1_  =  moint_comp->get_h1( test_ranges2, true ) ;
  shared_ptr<SMITH::Tensor_<double>> v2_  =  moint_comp->get_v2( test_ranges4 ) ;
  cout << " new_coeffs  v2_->norm() = " << v2_->norm() << endl; 

  set_target_state_info();
  set_ci_range_info();

  sys_info_ = make_shared<System_Info<double>>( targets_info_, true );

  cout << "initialized sys_info" << endl;
  shared_ptr< const PTree > ops_def_tree = idata->get_child_optional( "operators" ) ;
  if (ops_def_tree)
    get_new_ops_init( ops_def_tree ); 

  cout << "got ops_init" << endl;

  // Getting info about target expression (this includes which states are relevant)
  get_terms_init( idata->get_child( "terms" ) ); 
  cout << "got terms_init" << endl;

  // Get user specified variables (e.g. ranges, constant factors) which may appear in term definitions
  get_equations_init( idata->get_child( "equations" ) );
  cout << "got  equations_init" << endl;

  expression_map_ = sys_info_->expression_map;
  expression_machine_ = make_shared<SMITH::Expression_Computer::Expression_Computer<double>>( civectors_, expression_map_, range_conversion_map_, tensop_data_map_, 
                                                                                              gamma_data_map_, sigma_data_map_, civec_data_map_ );

  cout << "built expression machine" << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_expression_variables( shared_ptr<const PTree> variable_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << " PropTool::PropTool::get_expression_variables" << endl;

  //TODO don't do this...........
  inp_factor_map_->emplace("one", 1.0 ); 
  inp_factor_map_->emplace("1.0", 1.0 ); 
  inp_range_map_->emplace("none", make_shared<vector<int>>(0) ); 

  auto range_info_tree =  variable_def_tree->get_child_optional("ranges"); // can be mo or state.. 
  for ( auto& range_info : *range_info_tree ) {
   
    string range_name = range_info->get<string>( "name" ) ;
    auto range_vec_inp = range_info->get_child( "range_vec" ) ;
    shared_ptr<vector<int>> range_vec = make_shared<vector<int>>(0);
   
    for (auto& id : *range_vec_inp) 
      range_vec->push_back(lexical_cast<int>( id->data() ) );
   
    if ( inp_range_map_->find(range_name) != inp_range_map_->end() )
      throw runtime_error("Range \"" + range_name + "\" has been defined twice in input!!! ... aborting" ) ;
     
    inp_range_map_->emplace(range_name, range_vec); 

  }

  auto factor_info_tree = variable_def_tree->get_child_optional("factors"); //
  for ( auto& factor_info : *factor_info_tree ) {
   
    string factor_name = factor_info->get<string>( "name" ) ;
    bool factor_const = factor_info->get<bool>( "type" , true ) ;
    int  num_indexes = factor_info->get<int>( "num indexes" , 0 ) ; // TODO fix so we can deal with more than one
    
    if (factor_const ) {   
      if ( num_indexes == 0 ) { 

        double factor_value = factor_info->get<double>( "value" ) ;

        if ( inp_factor_map_->find(factor_name) != inp_factor_map_->end() )
          throw runtime_error("Factor \"" + factor_name + "\" has been defined twice in input!!! ... aborting" ) ;

        inp_factor_map_->emplace( factor_name, factor_value ); 

      } if ( num_indexes == 1 ) { 

        if ( inp_indexed_factor_map_->find(factor_name) != inp_indexed_factor_map_->end() )
          throw runtime_error("Factor range \"" + factor_name + "\" has been defined twice in input!!! ... aborting" ) ;

        auto factor_list_ptree =  factor_info->get_child("value"); 
        shared_ptr<vector<double>> factor_list = make_shared<vector<double>>(0);
        for (auto& factor_inp : *factor_list_ptree)
          factor_list->push_back(lexical_cast<double>(factor_inp->data()));

        inp_indexed_factor_map_->emplace( factor_name, factor_list ); 

      }
 
    } else {
      cout << "need to sort internally defined variables this (most urgently the eigenvalues of the Fock operator)" << endl;
    }
        
  }

  cout << "USER DEFINED FACTORS " << endl;
  for ( auto elem : *inp_factor_map_ ) 
    cout << elem.first << " = " << elem.second << endl; 

  cout << "USER DEFINED RANGES " << endl;
  for ( auto elem : *inp_range_map_ ) {
    cout << elem.first << " = [ " ;
    for (int orb_num :  *elem.second ) {
      cout << orb_num << " " ; cout.flush();
    }
    cout << "]" << endl;
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_new_ops_init( shared_ptr<const PTree> ops_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "void PropTool::PropTool::get_new_ops_init" << endl; 

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

    int state_dep =  op_def_inp->get<int>("state dependence");

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
    shared_ptr<TensOp::TensOp<double>> new_op = sys_info_->Build_TensOp( op_name, idxs_ptr, aops_ptr, ranges_ptr, symmfuncs, constraints, factor, TimeSymm, hconj, state_dep); 
    sys_info_->T_map->emplace( op_name, new_op );

  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_equations_init( shared_ptr<const PTree> equation_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::get_equations_init" << endl;

  for ( auto& equation_inp : *equation_def_tree ){
 
    string eqn_type = equation_inp->get<string>( "type" );

    if ( eqn_type == "LinearRM" ) { 
      get_equation_init_LinearRM( equation_inp ) ;

    } else {
      cout << "this equation type not implemented" << endl; 

    }
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_equation_init_LinearRM( shared_ptr<const PTree> equation_inp ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << " PropTool::PropTool::get_linear_equation_init_LinearRM" << endl;

  string eqn_name = equation_inp->get<string>( "name" );
  string eqn_target = equation_inp->get<string>( "target" );

  auto target_indices = make_shared<vector<string>>(0);
  auto ti_ptree = equation_inp->get_child("target indexes"); // Must solve for all "target" with these indices
  for (auto& si : *ti_ptree) 
    target_indices->push_back( lexical_cast<string>(si->data()));
  
  auto term_list = make_shared<vector<pair<string,shared_ptr<Term_Init>>>>();
  auto term_idrange_map_list = make_shared<vector<shared_ptr<map<string,pair<bool,string>>>>>();

  auto expression_def = equation_inp->get_child( "expression" );
  for (auto& term_info: *expression_def) {
  
    string term_name = term_info->get<string>( "term" );
    string term_factor = term_info->get<string>( "factor" );
    term_list->push_back(make_pair(term_factor, term_init_map_->at(term_name)));

    auto term_idrange_map = make_shared<map<string, pair<bool,string>>>();
    auto indexes_ptree =  term_info->get_child("indexes"); 
    for (auto& index_info : *indexes_ptree){
      string id_name = index_info->get<string>("name"); 
      string id_range = index_info->get<string>("range");
       
      bool   id_sum =  index_info->get<bool>("sum", false );
      term_idrange_map->emplace( id_name, make_pair(id_sum, id_range));
    }
    term_idrange_map->emplace( "none" , make_pair( false, "none" ) ); // TODO sort a better way of dealing with this case 
    term_idrange_map_list->push_back( term_idrange_map );     
    
  }
  auto master_expression = make_shared<Expression_Init>( term_list, term_idrange_map_list ); 
  auto eqn = make_shared<Equation_Init_LinearRM<double>>( eqn_name,  "LinearRM", master_expression, inp_range_map_, eqn_target, target_indices,
                                                          inp_factor_map_ );
  eqn->initialize_expression();

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_terms_init( shared_ptr<const PTree> term_inp_list ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::get_term_init" << endl;
   
  for ( auto& term_inp : *term_inp_list ) { 

    shared_ptr<map<string,int>> id_val_map = make_shared<map<string,int>>();

    auto get_id_ptr = [&id_val_map]( string id ){ 
       auto map_iter = id_val_map->find(id); 
       int* id_ptr;
       if ( map_iter == id_val_map->end() ){ 
         int new_id_val = 0; 
         id_val_map->emplace(id, new_id_val);
         id_ptr = &(id_val_map->at(id)) ; 
       } else { 
         id_ptr = &(map_iter->second); 
       }          
       return id_ptr;
    };

    string term_name = term_inp->get<string>( "name" );
    string term_type = term_inp->get<string>( "type", "full" );

    shared_ptr<vector<string>> braket_factors = make_shared<vector<string>>();
    auto braket_list = make_shared<vector<BraKet_Init>>(); 
    auto brakets_list_inp = term_inp->get_child("brakets"); 
    for ( auto& braket_inp : *brakets_list_inp ){
      auto bk_ops = make_shared<vector<Op_Init>>();

      auto bk_ops_inp = braket_inp->get_child("ops"); 
      for ( auto& op_def : *bk_ops_inp ){
        string opname = op_def->get<string>("name");
        vector<string> op_idxs(0);
        vector<int*> op_idxs_ptrs(0);

        auto ids_inp = op_def->get_child("ids");//note these are not orbital indexes
        for (auto& idx : *ids_inp ) { 
          op_idxs.push_back(lexical_cast<string>(idx->data()));
          op_idxs_ptrs.push_back( get_id_ptr(op_idxs.back()) );            
        }
 
        bk_ops->push_back(Op_Init( opname, op_idxs, make_shared<vector<int*>>(op_idxs_ptrs) ));
      }

      string bra_index = braket_inp->get<string>("bra");
      int* bra_index_ptr = get_id_ptr(bra_index);

      string ket_index = braket_inp->get<string>("ket");
      int* ket_index_ptr = get_id_ptr(ket_index);

      braket_factors->push_back(braket_inp->get<string>("factor", "one"));
            
      braket_list->push_back( BraKet_Init( bk_ops, bra_index, bra_index_ptr, ket_index, ket_index_ptr ));
    }
    auto new_term = make_shared<Term_Init>( term_name, term_type, braket_list, braket_factors, id_val_map );  

    term_init_map_->emplace( term_name, new_term );
  }

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_ao_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "PropTool::PropTool::set_ao_range_info" << endl;

  closed_rng_  = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nclosed_-ncore_, maxtile_, 0, ncore_));
  active_rng_  = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nact_, min((size_t)10,maxtile_), closed_rng_->nblock(), ncore_ + closed_rng_->size()));
  virtual_rng_ = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nvirt_, maxtile_, closed_rng_->nblock()+ active_rng_->nblock(), ncore_+closed_rng_->size()+active_rng_->size()));
  free_rng_    = make_shared<SMITH::IndexRange>(*closed_rng_);
  free_rng_->merge(*active_rng_);
  free_rng_->merge(*virtual_rng_);

  not_closed_rng_  = make_shared<SMITH::IndexRange>(*active_rng_); not_closed_rng_->merge(*virtual_rng_);
  not_active_rng_  = make_shared<SMITH::IndexRange>(*closed_rng_); not_active_rng_->merge(*virtual_rng_);
  not_virtual_rng_ = make_shared<SMITH::IndexRange>(*closed_rng_); not_virtual_rng_->merge(*active_rng_);

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
void PropTool::PropTool::build_op_tensors( vector<string>& expression_list ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "PropTool::PropTool::build_op_tensors" << endl;

  // Creating tensors from existing matrices; seperate loop as must run through all states first to make proper use of symmetry
  for (string expression_name : expression_list ) {
    for ( auto tensop_it : *(sys_info_->T_map) ) {
      vector< shared_ptr< const vector<string>>> unique_range_blocks = *(tensop_it.second->unique_range_blocks());
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
