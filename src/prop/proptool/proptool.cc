#include <src/prop/proptool/proptool.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/debugging_utils.h>

using namespace std;
using namespace bagel;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PropTool::PropTool::PropTool( shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r ): 
                              idata_(idata), geom_(g), ref_(r), ciwfn_(ref_->ciwfn()), civectors_(ciwfn_->civectors())  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::PropTool" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sort out how to determine datatype!!
  inp_factor_map_ = make_shared<map<string, double>>();
  inp_indexed_factor_map_ = make_shared<map<string, shared_ptr<vector<double>>>>(); // TODO sort this
  inp_range_map_ = make_shared<map<string, shared_ptr<vector<int>>>>();

  set_primes();
  range_conversion_map_ = make_shared<map<string, shared_ptr<SMITH::IndexRange>>>();
  range_prime_map_ = make_shared<map<char,long unsigned int>>();

  term_init_map_ = make_shared<map<string, shared_ptr<Term_Init>>>();
  expression_init_map_ = make_shared<map<string, shared_ptr<Expression_Init>>>();
  equation_init_map_ = make_shared<map<string, shared_ptr<Equation_Init_Base>>>();

  read_input_and_initialize(); 

  build_algebraic_task_lists( idata_->get<string>("equation interdependence", "share" ) ) ;

  execute_compute_lists();  

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::execute_compute_lists(){  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::execute_compute_lists()" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 for ( string& equation_name : equation_execution_list_ ) 
   system_computer_->build_equation_computer( equation_name );
 
 return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// eqn_dependence = "share" :  equations can be done in any order and share information
// eqn_dependence = "sequential" :  equations must be done in the order they appear in the input and share information
// eqn_dependence = "noshare" :  equations can be done in any order, but do not share information
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::build_algebraic_task_lists( string  eqn_interdependence ){  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::build_algebraic_task_lists()" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  if ( eqn_interdependence == "share"  ) { 
    for ( string& equation_name : equation_execution_list_ ) 
       sys_info_->construct_equation_task_list( equation_name ) ; 
  } else { 
     throw logic_error( "form of equation interdependence \"" + eqn_interdependence +"\" not implemented yet" ) ; 
  } 
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::read_input_and_initialize(){  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "void PropTool::PropTool::read_input_and_initialize()" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Get user specified variables (e.g. ranges, constant factors) which may appear in term definitions
  get_expression_variables( idata_->get_child("variables") );
 
  // gets information about the wavefunction from the refence
  get_wavefunction_info();

  // build system information object ( for algebraic task list construction)
  sys_info_ = make_shared<System_Info<double>>( targets_info_, true );
  sys_info_->range_prime_map_ = range_prime_map_;

  // build system computer (for computational task list construction/execution)
  auto moint_init = make_shared<MOInt_Init<double>>( geom_, dynamic_pointer_cast<const Reference>(ref_), ncore_,
                                                     nfrozenvirt_, block_diag_fock_ );

  auto moint_computer = make_shared<MOInt_Computer<double>>( moint_init, range_conversion_map_ );
  
  vector<string> free2 = { "free" , "free" };
  moint_computer->get_fock( free2, true );

  //TODO should build gamma_computer inside system_computer, like this due to DVec class dependence of B_Gamma_Computer 
  auto gamma_computer = make_shared<B_Gamma_Computer::B_Gamma_Computer<double>>(civectors_); 

  system_computer_ = make_shared<System_Computer::System_Computer<double>>(sys_info_, moint_computer, range_conversion_map_, gamma_computer );

  shared_ptr< const PTree > ops_def_tree = idata_->get_child_optional( "operators" ) ;
  if (ops_def_tree)
    get_new_ops_init( ops_def_tree ); 

  // Getting info about target expression (this includes which states are relevant)
  get_terms_init( idata_->get_child( "terms" ) ); 

  get_equations_init( idata_->get_child( "equations" ) );

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_wavefunction_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::get_wavefunction_info()" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Initializing range sizes either from idate or reference wfn 
  maxtile_   = idata_->get<int>("maxtile", 10);

  //cimaxtile_ = idata_->get<int>("cimaxtile", (ciwfn_->civectors()->size() > 10000) ? 100 : 10);
  cimaxtile_ = 100000; //TODO fix this so it uses the above statement, issue in b_gamma_computer means must use large cimaxblock for now

  const bool frozen = idata_->get<bool>("frozen", true);
  ncore_ = idata_->get<int>("ncore", (frozen ? ref_->geom()->num_count_ncore_only()/2 : 0));
  if (ncore_)
    cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;
  nfrozenvirt_ = idata_->get<int>("nfrozenvirt", 0);
  if (nfrozenvirt_)
    cout << "    * freezing " << nfrozenvirt_ << " orbital" << (nfrozenvirt_^1 ? "s" : "") << " (virtual)" << endl;

  nclosed_  = idata_->get<int>( "nclosed" , ref_->nclosed());
  nact_     = idata_->get<int>( "nact"  , ref_->nact());
  nvirt_    = idata_->get<int>( "nvirt" , ref_->nvirt());
  nocc_     = nclosed_ + nact_; 
  nfrozenvirt_ = idata_->get<int>( "nfrozenvirt", 0 );

  // TODO : should be determined from summation ranges in expression
  nstates_ = idata_->get<int>( "nstates" , ciwfn_->nstates() );
  target_states_ = vector<int>(nstates_); 
  std::iota(target_states_.begin(), target_states_.end(), 0 ); 
     
  // leave for now
  block_diag_fock_ = false;
  gaunt_    = false;
  breit_    = false;

  set_ao_range_info();
  //creates the ci_info from the reference wavefunction
  set_ci_range_info();

  set_target_state_info();

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_expression_variables( shared_ptr<const PTree> variable_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
 cout << " PropTool::PropTool::get_expression_variables" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //TODO don't do this...........
  inp_factor_map_->emplace("one", 1.0 ); 
  inp_factor_map_->emplace("1.0", 1.0 ); 
  inp_range_map_->emplace("none", make_shared<vector<int>>(1, -1) ); 

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
          factor_list->push_back( lexical_cast<double>( factor_inp->data() ) );

        inp_indexed_factor_map_->emplace( factor_name, factor_list ); 

      }
 
    } else {
      cout << "need to sort internally defined variables this (most urgently the eigenvalues of the Fock operator)" << endl;
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_new_ops_init( shared_ptr<const PTree> ops_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "void PropTool::PropTool::get_new_ops_init" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //this shouldn't be necessary, but the conversion doesn't seem to work otherwise....
  auto conv_to_bool = []( int aop ) { return aop == 1 ? true : false ; };

  // keep this so you don't have to get names correct, needs modification for spin case
  auto conv_to_range = []( string rng ) {
           if( rng[0] == 'c' ) return "c";
           if( rng[0] == 'a' ) return "a";
           if( rng[0] == 'v' ) return "v";
      };

  //TODO  add in user defined ranges (not just on the operators; actual range extents using orbital numbers
  //      create these new ranges, put them into the range conversion map, and assign them a prime.
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
      
    vector<shared_ptr<Transformation>> symmfuncs(0);
    vector<shared_ptr<Constraint>> constraints(0);
    
    cout << endl << " user defined operator : " << op_name << endl << endl;
    shared_ptr<TensOp::TensOp<double>> new_op = sys_info_->Build_TensOp(op_name, idxs_ptr, aops_ptr, ranges_ptr, symmfuncs, constraints, factor, TimeSymm, hconj, state_dep);
    sys_info_->MT_map()->emplace( op_name, new_op );

  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_equations_init( shared_ptr<const PTree> equation_def_tree ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::get_equations_init" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( auto& equation_inp : *equation_def_tree ){
 
    string eqn_type = equation_inp->get<string>( "type" );

    if ( eqn_type == "LinearRM" ) { 
      get_equation_init_LinearRM( equation_inp ) ;

    } else if ( eqn_type == "Value" ) { 
      get_equation_init_Value( equation_inp ) ;
    
    } else { 
      cout << "this equation type not implemented" << endl; 

    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_equation_init_Value( shared_ptr<const PTree> equation_inp ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << " PropTool::PropTool::get_equation_init_Value" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string eqn_name = equation_inp->get<string>( "name" );

  // if target indices is not specified just get a different value for every combination of non-summed indexes
  auto target_indices = make_shared<vector<string>>(0);
  auto ti_ptree = equation_inp->get_child_optional("target indexes"); // Must solve for all "target" with these indices
  if (ti_ptree)
    for (auto& si : *ti_ptree) 
      target_indices->push_back( lexical_cast<string>(si->data()));

  auto term_list = make_shared<vector<pair<string,shared_ptr<Term_Init>>>>();
  auto term_idrange_map_list = make_shared<vector<shared_ptr<map<string,pair<bool,string>>>>>();

  auto expression_def = equation_inp->get_child( "expression" );

  string expression_type = "full"; 
  for (auto& term_info: *expression_def) {
     
    string term_name   = term_info->get<string>( "term" );
    string term_type   = term_info->get<string>( "type" );
    string term_factor = term_info->get<string>( "factor" );
    
    if (term_type[0] == 'o' )
      expression_type = "orb_excitation_derivative"; 
    
    shared_ptr<Term_Init> new_term_init = term_init_map_->at(term_name);

    term_list->push_back(make_pair(term_factor, new_term_init));

    auto term_idrange_map = make_shared<map<string, pair<bool,string>>>();
    auto indexes_ptree =  term_info->get_child("indexes"); 
    for (auto& index_info : *indexes_ptree){
      string id_name = index_info->get<string>("name"); 
      string id_range = index_info->get<string>("range");
      bool id_sum =  index_info->get<bool>("sum", false );
      term_idrange_map->emplace( id_name, make_pair(id_sum, id_range));
    }

    term_idrange_map->emplace( "none" , make_pair( false, "none" ) ); // TODO sort a better way of dealing with this case 
    term_idrange_map_list->push_back( term_idrange_map );     
    
  }
  auto master_expression = make_shared<Expression_Init>( term_list, term_idrange_map_list, expression_type ); 
  auto eqn_init = make_shared<Equation_Init_Value<double>>( eqn_name, "Value", master_expression, inp_range_map_, target_indices, inp_factor_map_ );
  eqn_init->initialize_expressions();
 
  sys_info_->create_equation( eqn_name, "Value", eqn_init->term_braket_map_ , eqn_init->expression_term_map_ );

  equation_execution_list_.push_back(eqn_name); 

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_equation_init_LinearRM( shared_ptr<const PTree> equation_inp ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << " PropTool::PropTool::get_linear_equation_init_LinearRM" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string eqn_name = equation_inp->get<string>( "name" );
  string eqn_target = equation_inp->get<string>( "target" );
  auto target_indices = make_shared<vector<string>>(0);
  auto ti_ptree = equation_inp->get_child("target indexes"); // Must solve for all "target" with these indices

  for (auto& si : *ti_ptree)
    target_indices->push_back( lexical_cast<string>(si->data()));

  auto term_list = make_shared<vector<pair<string,shared_ptr<Term_Init>>>>();
  auto term_idrange_map_list = make_shared<vector<shared_ptr<map<string,pair<bool,string>>>>>();
  auto expressions_inp = equation_inp->get_child( "expression" );
  auto expression_init_list  = make_shared<vector<shared_ptr<Expression_Init>>>();
  string expression_type = "full";
  for ( auto& expression_def : *expressions_inp ) {
    for ( auto& term_info : *expression_def) {
  
      string term_name = term_info->get<string>( "term" );
      string term_factor = term_info->get<string>( "factor" );

      shared_ptr<Term_Init> term_init = term_init_map_->at(term_name);
      if ( term_init->orbital_projector_ )
        expression_type = "orb_excitation_derivative";

      term_list->push_back( make_pair(term_factor, term_init ) );
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
    //TODO clean up these classes, names are confusing; initialization is badly scrambled due to changes in how to deal with indexes and varying term
    //     types in the expression_computer. 
    expression_init_list->push_back(make_shared<Expression_Init>( term_list, term_idrange_map_list, expression_type )); 

  }
    
  auto eqn_init = make_shared<Equation_Init_LinearRM<double>>( eqn_name, "LinearRM", expression_init_list, inp_range_map_, eqn_target, target_indices,
                                                               inp_factor_map_ );
  eqn_init->initialize_all_terms();

  sys_info_->create_equation( eqn_name, "LinearRM", eqn_init->term_braket_map_, eqn_init->expression_term_map_,
                              eqn_init->term_braket_map_state_spec_, eqn_init->expression_term_map_state_spec_ );

  equation_execution_list_.push_back(eqn_name); 
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::get_terms_init( shared_ptr<const PTree> term_inp_list ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
  cout << "PropTool::PropTool::get_term_init" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
  int zero = 0;
  for ( auto& term_inp : *term_inp_list ) { 

    shared_ptr<map<string,int>> id_val_map = make_shared<map<string,int>>();

    auto get_id_ptr = [&id_val_map, &zero ]( string id ){ 
         id_val_map->emplace(id, zero);
         return &(id_val_map->at(id)); 
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
        shared_ptr<vector<int*>> op_idxs_ptrs= make_shared<vector<int*>>(0);

        auto ids_inp = op_def->get_child("ids");//note these are not orbital indexes
        for (auto& idx : *ids_inp ) { 
          op_idxs.push_back(lexical_cast<string>(idx->data()));
          id_val_map->emplace(op_idxs.back(), zero); 
          op_idxs_ptrs->push_back( &( id_val_map->at(op_idxs.back()) ) ); 
        }
 
        if ( op_def->get_child_optional( "transform" ) ) { //TODO sort get_optional properly; won't work now for some reason... 
          string trans = op_def->get<string>("transform");
          bk_ops->push_back(Op_Init( opname, op_idxs, op_idxs_ptrs, trans ));
        } else {
          bk_ops->push_back(Op_Init( opname, op_idxs, op_idxs_ptrs ));
        }
      }

      string bra_index = braket_inp->get<string>("bra");
      int* bra_index_ptr = get_id_ptr(bra_index);

      string ket_index = braket_inp->get<string>("ket");
      int* ket_index_ptr = get_id_ptr(ket_index);

      braket_factors->push_back(braket_inp->get<string>("factor", "one"));
            
      braket_list->push_back( BraKet_Init( bk_ops, bra_index, bra_index_ptr, ket_index, ket_index_ptr ));
    }
   
    shared_ptr<Term_Init>  new_term = make_shared<Term_Init>( term_name, term_type, braket_list, braket_factors, id_val_map); 
    if ( term_type == "orb" ) 
      new_term->proj_op_name_ = term_inp->get<string>( "target op" , "X" );

    term_init_map_->emplace( term_name, new_term );
  }

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_ao_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::set_ao_range_info" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //TODO should also read user defined ranges from input!
  closed_rng_  = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nclosed_-ncore_, maxtile_, 0, ncore_));
  active_rng_  = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nact_, min((size_t)10,maxtile_), closed_rng_->nblock(), ncore_ + closed_rng_->size()));
  virtual_rng_ = make_shared<SMITH::IndexRange>(SMITH::IndexRange(nvirt_, maxtile_, closed_rng_->nblock()+ active_rng_->nblock(), ncore_+closed_rng_->size()+active_rng_->size()));
  free_rng_    = make_shared<SMITH::IndexRange>(*closed_rng_);
  free_rng_->merge(*active_rng_);
  free_rng_->merge(*virtual_rng_);

  not_closed_rng_  = make_shared<SMITH::IndexRange>(*active_rng_); not_closed_rng_->merge(*virtual_rng_);
  not_active_rng_  = make_shared<SMITH::IndexRange>(*closed_rng_); not_active_rng_->merge(*virtual_rng_);
  not_virtual_rng_ = make_shared<SMITH::IndexRange>(*closed_rng_); not_virtual_rng_->merge(*active_rng_);

  range_conversion_map_->emplace("c", closed_rng_); 
  range_conversion_map_->emplace("a", active_rng_);
  range_conversion_map_->emplace("v", virtual_rng_);
  range_conversion_map_->emplace("free", free_rng_);
 
  range_conversion_map_->emplace("notcor", not_closed_rng_);
  range_conversion_map_->emplace("notact", not_active_rng_);
  range_conversion_map_->emplace("notvir", not_virtual_rng_); 


  // TODO  sort this out properly
  vector<char> range_list;
  bool spinfree = false;
  if ( !spinfree ) {
    range_list = vector<char>( { 'v', 'V', 'c', 'C', 'a', 'A' } ); 
  }  else  {
    range_list = vector<char>( { 'v', 'c', 'a' } ); 
  }
  for ( char elem : range_list ) {
    range_prime_map_->emplace( elem , range_primes_.back());
    range_primes_.pop_back();
  }

  return;                    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_ci_range_info() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::set_ci_range_info" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( int ii : target_states_ )
    range_conversion_map_->emplace( get_civec_name( ii , civectors_->data(ii)->det()->norb(), civectors_->data(ii)->det()->nelea(), civectors_->data(ii)->det()->neleb()),
                                                    make_shared<SMITH::IndexRange>(civectors_->data(ii)->det()->size(), cimaxtile_ ));  
    
 return;
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void  PropTool::PropTool::identify_degeneracies( const vector<double>& energies ) {
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
  cout << "void PropTool::PropTool::set_target_state_info()" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int state_num = 0;
  double energy_buff = energies.front();
  degenerate_states_ = vector<vector<int>>(0);
 
  if ( energies.size() > 1 ) {
    vector<int> degenerate_set(0); 
    for ( vector<double>::const_iterator e_it = (energies.begin()+1) ; e_it != energies.end(); e_it++, state_num++ ) {
      if ( ( abs(energy_buff - *e_it) > 0.0000001 ) || ( e_it == energies.end() ) ) {
        degenerate_states_.push_back(degenerate_set);
        degenerate_set = vector<int>(0); 
      } else {
        degenerate_states_.push_back( degenerate_set );
      }
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_target_state_info() {
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::set_target_info" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  targets_info_ = make_shared<StatesInfo<double>> ( target_states_ ) ;

  for ( int state_num : target_states_ ){ 
     shared_ptr<map<char,int>> elec_range_map = make_shared<map<char,int>>(); 
     elec_range_map->emplace('c', nclosed_ - ncore_);
     elec_range_map->emplace('C', nclosed_ - ncore_);
     elec_range_map->emplace('A', civectors_->data(state_num)->det()->nelea());
     elec_range_map->emplace('a', civectors_->data(state_num)->det()->neleb());
     elec_range_map->emplace('v', 0 );
     elec_range_map->emplace('V', 0 );
  
     shared_ptr<map<char,int>> hole_range_map = make_shared<map<char,int>>(); 
     hole_range_map->emplace('c', 0);
     hole_range_map->emplace('C', 0);
     hole_range_map->emplace('a', nact_ - civectors_->data(state_num)->det()->nelea() );
     hole_range_map->emplace('A', nact_ - civectors_->data(state_num)->det()->neleb() );
     hole_range_map->emplace('v', 100 ); // TODO this is almost certainly always OK, but should be set properly...
     hole_range_map->emplace('V', 100 );

     targets_info_->add_state( nact_, civectors_->data(state_num)->det()->nelea() + civectors_->data(state_num)->det()->neleb(), state_num,
                               elec_range_map, hole_range_map); 
  }
 
  for ( auto& elem : targets_info_->civec_info_map_ )
    elem.second->set_elec_hole_pnums( range_prime_map_ );
  
  targets_info_->range_prime_map_ = range_prime_map_;

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<SMITH::IndexRange>> PropTool::PropTool::convert_to_indexrange( shared_ptr<const vector<string>> range_block_str ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
  cout << "PropTool::PropTool::convert_to_indexrange" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<SMITH::IndexRange>> range_block = make_shared<vector<SMITH::IndexRange>>( range_block_str->size() );
  for( int ii = 0; ii != range_block_str->size() ; ii++ ) 
    range_block->at(ii) = *(range_conversion_map_->at( range_block_str->at(ii) ));

  return range_block; 

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::print_input_info() { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
cout << "PropTool::PropTool::print_input_info()" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "range_prime_map" << endl;  
  for ( auto& elem : *range_prime_map_ ) {
    string rng = "";
    rng.push_back(elem.first);
    cout << rng << " : " << elem.second << endl; 
  }
 
  for ( auto& elem : targets_info_->civec_info_map_ )
    cout << elem.second->name() << " elec_pnum = " << elem.second->elec_pnum() << " hole_pnum = " << elem.second->hole_pnum() << endl;  

  cout <<  endl << "=========== Range map ============ " << endl;
  for ( auto& elem : *range_conversion_map_ ){ 
    Debugging_Utils::print_sizes(  elem.second->range(),elem.first );
    cout << " total_size = " << elem.second->size() << endl;
  }

  cout << endl << "============ user defined factors =============== " << endl;
  for ( auto elem : *inp_factor_map_ ) 
    cout << elem.first << " = " << elem.second << endl; 

  cout << endl << "============ user defiend ranges ================ " << endl;
  for ( auto elem : *inp_range_map_ ) {
    cout << elem.first << " = [ " ;
    for (int orb_num :  *elem.second ) {
      cout << orb_num << " " ; cout.flush();
    }
    cout << "]" << endl;
  }

  cout << endl << "======equation_execution_list_====== " << endl;
  for ( auto elem : equation_execution_list_ ) 
    cout << elem << endl;

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PropTool::PropTool::set_primes() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_DRIVER
  cout << "PropTool::PropTool::set_primes" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // for encoding ranges and contractions
  range_primes_ =  { 1013, 1009, 997, 991, 983, 977, 971, 967, 953, 947, 941, 937, 929, 919, 911, 907, 887, 883, 881, 877,
                      863, 859, 857, 853, 839, 829, 827, 823, 821, 811, 809, 797, 787, 773, 769, 761, 757, 751, 743, 739, 733,
                      727, 719, 709, 701, 691, 683, 677, 673, 661, 659, 653, 647, 643, 641, 631, 619, 617, 613, 607, 601, 599,
                      593, 587, 577, 571, 569, 563, 557, 547, 541, 523, 521, 509, 503, 499, 491, 487, 479, 467, 463, 461, 457,
                      449, 443, 439, 433, 431, 421, 419, 409, 401, 397, 389, 383, 379, 373, 367, 359, 353, 349, 347, 337, 331,
                      317, 313, 311, 307, 293, 283, 281, 277, 271, 269, 263, 257, 251, 241, 239, 233, 229, 227, 223, 211, 199,
                      197, 193, 191, 181, 179, 173, 167, 163, 157, 151, 149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97,
                      89, 83, 79, 73, 71, 67, 61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2 } ;

  ctr_primes_ =  { 1013, 1009, 997, 991, 983, 977, 971, 967, 953, 947, 941, 937, 929, 919, 911, 907, 887, 883, 881, 877,
                    863, 859, 857, 853, 839, 829, 827, 823, 821, 811, 809, 797, 787, 773, 769, 761, 757, 751, 743, 739, 733,
                    727, 719, 709, 701, 691, 683, 677, 673, 661, 659, 653, 647, 643, 641, 631, 619, 617, 613, 607, 601, 599,
                    593, 587, 577, 571, 569, 563, 557, 547, 541, 523, 521, 509, 503, 499, 491, 487, 479, 467, 463, 461, 457,
                    449, 443, 439, 433, 431, 421, 419, 409, 401, 397, 389, 383, 379, 373, 367, 359, 353, 349, 347, 337, 331,
                    317, 313, 311, 307, 293, 283, 281, 277, 271, 269, 263, 257, 251, 241, 239, 233, 229, 227, 223, 211, 199,
                    197, 193, 191, 181, 179, 173, 167, 163, 157, 151, 149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97,
                    89, 83, 79, 73, 71, 67, 61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2 } ;


   // contract_number (0,4)(1,3)(7,8) =  primes_[ (0*1 +4*aops->size() ]   
   // range_number_ = just check through defined ranges, pop from back  
   // allowed_contraction_number = product of all contraction numbers (not as big as you'd think, and still faster than lookup);
   // vector<bool> allowed_ctrs = ( aops->size() * aops->size() );
   // bool allowed_contraction( ii, jj ) {  allowed_contraction_ctrs[ ii + jj*aops->size()];  } is simpler. 
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
