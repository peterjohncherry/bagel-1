#include <bagel_config.h>
#include  <src/prop/proptool/algebraic_manipulator/system_info.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Note, spinfree should tell us not just if the wavefunction is free, but whether or 
//not the perturbation being applied is spin independent
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
System_Info<DataType>::System_Info( shared_ptr<StatesInfo<DataType>> states_info , bool spinfree ) :
                                                 states_info_(states_info), spinfree_(spinfree) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  braket_map_       = make_shared< map <string, shared_ptr<vector<shared_ptr< TensOp::TensOp<DataType>>>>>>();
  term_braket_map_ = make_shared<map<string, shared_ptr<vector<BraKet<DataType>>>>>();

  expression_term_map_ = make_shared<map<string, shared_ptr<vector<pair<DataType, string>>>>>();
  expression_map       = make_shared< map <string, shared_ptr<Expression<DataType>>>>();
  
  equation_map_       = make_shared< map <string, shared_ptr<Equation_Base<DataType>>>>();

  MT_map_         = make_shared< map <string, shared_ptr<TensOp_Base>>>();
  CTP_map_        = make_shared< map <string, shared_ptr<CtrTensorPart_Base>>>();
  ACompute_map   = make_shared< map <string, shared_ptr<vector<shared_ptr<CtrOp_base>> >>>();
  Gamma_map      = make_shared< map <string, shared_ptr<GammaInfo> > >();

  if (spinfree_ /*TODO like this for testing; obviously must put back*/ || !spinfree_ ) { 
    cout << " setting spinfree ranges" << endl;
    free     = {"cor", "act", "vir"};
    not_core = {"act", "vir"};
    not_act  = {"cor", "vir"};
    not_virt = {"cor", "act"};
    core     = {"cor"};
    act      = {"act"};
    virt     = {"vir"};
  } else { 
    free     = {"corA", "actA", "virA", "corB", "actB", "virB"};
    not_core = {"actA", "virA", "actB", "virB"};
    not_act  = {"corA", "virA", "corB", "virB"};
    not_virt = {"corA", "actA", "corB", "actB"};
    core     = {"corA", "corB"};
    act      = {"actA", "actB"};
    virt     = {"virA", "virB"};
  }

  // for encoding ranges and contractions
  primes_ranges_ =  { 1013, 1009, 997, 991, 983, 977, 971, 967, 953, 947, 941, 937, 929, 919, 911, 907, 887, 883, 881, 877,
                        863, 859, 857, 853, 839, 829, 827, 823, 821, 811, 809, 797, 787, 773, 769, 761, 757, 751, 743, 739, 733,
                        727, 719, 709, 701, 691, 683, 677, 673, 661, 659, 653, 647, 643, 641, 631, 619, 617, 613, 607, 601, 599,
                        593, 587, 577, 571, 569, 563, 557, 547, 541, 523, 521, 509, 503, 499, 491, 487, 479, 467, 463, 461, 457,
                        449, 443, 439, 433, 431, 421, 419, 409, 401, 397, 389, 383, 379, 373, 367, 359, 353, 349, 347, 337, 331,
                        317, 313, 311, 307, 293, 283, 281, 277, 271, 269, 263, 257, 251, 241, 239, 233, 229, 227, 223, 211, 199,
                        197, 193, 191, 181, 179, 173, 167, 163, 157, 151, 149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97,
                        89, 83, 79, 73, 71, 67, 61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2 } ;

  primes_ctrs_ =  { 1013, 1009, 997, 991, 983, 977, 971, 967, 953, 947, 941, 937, 929, 919, 911, 907, 887, 883, 881, 877,
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

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::construct_equation_task_list( string equation_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Info<DataType>::System_Info::construct_equation_task_list : " << equation_name << endl;

  if (equation_map_->at( equation_name)->type() == "Value" ) { 
    equation_map_->at( equation_name)->generate_all_expressions();
  } else if (equation_map_->at( equation_name)->type() == "LinearRM" ) { 
    equation_map_->at( equation_name)->generate_state_specific_terms();
  } 
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<TensOp::TensOp<DataType>>
System_Info<DataType>::Build_TensOp( string op_name,
                                     shared_ptr<vector<string>> op_idxs,
                                     shared_ptr<vector<bool>> op_aops, 
                                     shared_ptr<vector<vector<string>>> op_idx_ranges,
                                     vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                                     vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                                     DataType factor, string Tsymmetry, bool hconj,  int state_dep ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "System_Info<DataType>::System_Info::Build_TensOp" <<   endl;

  //NOTE: change to use proper factor
  int tmpfac = 1;
  shared_ptr<TensOp::TensOp<DataType>> new_op = make_shared<TensOp::TensOp<DataType>>( op_name, *op_idxs, *op_idx_ranges, *op_aops,
                                                                                        tmpfac,  Symmetry_Funcs, Constraint_Funcs, Tsymmetry, state_dep);
  
  // change to be state specific
  cout << "getting  ctr tens ranges for New_Op : " << op_name << endl;
  new_op->get_ctrs_tens_ranges();
  CTP_map_->insert( new_op->CTP_map()->begin(), new_op->CTP_map()->end());
  cout << "got ctr tens ranges for new_op : " << op_name << endl;

  return new_op;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::create_equation( std::string name, std::string type, 
                                        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                                        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map, 
                                        std::shared_ptr<std::map<std::pair< std::string, std::vector<std::pair<std::string, int>>>, 
                                                                            std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_state_spec, 
                                        std::shared_ptr<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>, 
                                                                  std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_state_spec ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Info<DataType>::System_Info::create_equation" << endl; 

  term_braket_map_->insert(term_braket_map->begin(), term_braket_map->end()) ;
  expression_term_map_->insert(expression_term_map->begin(), expression_term_map->end()) ;
  
  shared_ptr<Equation_Base<DataType>>  new_eqn;
  if ( type == "Value" ) { 

    auto new_eqn_val = make_shared<Equation_Value<DataType>>( name, type, states_info_, term_braket_map, expression_term_map );
    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>( new_eqn_val );

    if (!new_eqn) { throw runtime_error("cast from Equation_Value to equation_base failed" ); }

    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map, MT_map_, CTP_map_ );
    equation_map_->emplace( name, new_eqn ); 

  } else if ( type == "LinearRM") { 

    shared_ptr<Equation_LinearRM<DataType>> new_eqn_lrm  =
    make_shared<Equation_LinearRM<DataType>>( name, type, states_info_, term_braket_map, expression_term_map,
                                              term_braket_map_state_spec, expression_term_map_state_spec );

    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>(new_eqn_lrm);

    if (!new_eqn) { throw runtime_error("cast from Equation_LinearRM to Equation_Base failed" ); }

    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map,  MT_map_, CTP_map_ );
    new_eqn_lrm->generate_state_specific_terms();
    equation_map_->emplace( name, new_eqn); 
  
  } else {  
    throw logic_error( "equation type \""+ type + "\" not implemented yet! Aborting!"); 

  }
  cout << "put equation " << name << " into map" <<endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::create_equation( std::string name, std::string type, 
                                        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                                        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Info<DataType>::System_Info::create_equation " << endl; 

  term_braket_map_->insert(term_braket_map->begin(), term_braket_map->end()) ;
  expression_term_map_->insert(expression_term_map->begin(), expression_term_map->end()) ;
  
  shared_ptr<Equation_Base<DataType>>  new_eqn;
  if ( type == "Value" ) { 
    shared_ptr<Equation_Value<DataType>> new_eqn_val  = make_shared<Equation_Value<DataType>> ( name, type, states_info_,  term_braket_map, expression_term_map );
    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>(new_eqn_val);
    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map, MT_map_, CTP_map_ );
    equation_map_->emplace( name, new_eqn); 

  } else if ( type == "LinearRM") { 
  
    throw logic_error( "Must provide state specific term map for doing linearRM!! Aborting!! " ) ; 
    
  } else {  
    throw logic_error( "equation type \""+ type + "\" not implemented yet! Aborting!"); 

  }

  return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
