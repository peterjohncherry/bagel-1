#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/system_info.h>
#include <src/prop/proptool/proputils.h>

using namespace std;

#define __DEBUG_PROPTOOL_SYSTEM_INFO
#define __DEBUG_PROPTOOL_EQUATION_BASE
#define __DEBUG_PROPTOOL_EXPRESSION_FULL
#define __DEBUG_PROPTOOL_BRAKET_BASE
#define __DEBUG_PROPTOOL_BRAKET_FULL
#define __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE
#define __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Note, spinfree should tell us not just if the wavefunction is free, but whether or 
//not the perturbation being applied is spin independent
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
System_Info<DataType>::System_Info( shared_ptr<StatesInfo<DataType>> states_info , bool spinfree ) :
                                                 states_info_(states_info), spinfree_(spinfree) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __DEBUG_PROPTOOL_SYSTEM_INFO
cout << "System_Info::System_Info" <<endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////

  braket_map_       = make_shared< map <string, shared_ptr<vector<shared_ptr< TensOp::TensOp<DataType>>>>>>();
  term_braket_map_ = make_shared<map<string, shared_ptr<vector<shared_ptr<BraKet_Base>>>>>();

  expression_term_map_ = make_shared<map<string, shared_ptr<vector<pair<DataType, string>>>>>();
  expression_map       = make_shared< map <string, shared_ptr<Expression<DataType>>>>();
  
  equation_map_       = make_shared< map <string, shared_ptr<Equation_Base<DataType>>>>();

  MT_map_        = make_shared< map <string, shared_ptr<TensOp_Base>>>();
  CTP_map_       = make_shared< map <string, shared_ptr<CtrTensorPart_Base>>>();
  ACompute_map   = make_shared< map <string, shared_ptr<vector<shared_ptr<CtrOp_base>>>>>();
  Gamma_map      = make_shared< map <string, shared_ptr<GammaInfo_Base>>>();

  if (spinfree_ ) { 
    free     = {"c", "a", "v"};
    not_core = {"a", "v"};
    not_act  = {"c", "v"};
    not_virt = {"c", "a"};
    core     = {"c"};
    act      = {"a"};
    virt     = {"v"};
  } else { 
    free     = {"c", "C", "a", "A", "v", "V" };
    not_core = {"a", "A", "v", "V"};
    not_act  = {"c", "C", "v", "V"};
    not_virt = {"c", "C", "a", "A"};
    core     = {"c", "C"};
    act      = {"a", "A"};
    virt     = {"v", "V"};
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::construct_equation_task_list( string equation_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __DEBUG_PROPTOOL_SYSTEM_INFO
  cout << "System_Info<DataType>::System_Info::construct_equation_task_list : " << equation_name << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                                     vector<shared_ptr<Transformation>>& symmfuncs,
                                     vector<shared_ptr<Constraint>>& constraints,
                                     DataType factor, string Tsymmetry, bool hconj,  int state_dep  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __DEBUG_PROPTOOL_SYSTEM_INFO
cout << "System_Info<DataType>::System_Info::Build_TensOp" <<   endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //NOTE: change to use proper factor ( the factor is the factor by which the Re/Im parts of the Tensop are multiplied, not the Re/Im parts of the factor itself) 
  pair<double,double> tmp_fac = make_pair(1.0 , 1.0 );

  shared_ptr<TensOp::TensOp<DataType>> new_op = make_shared<TensOp::TensOp<DataType>>( op_name, *op_idxs, *op_idx_ranges, *op_aops,
                                                                                        tmp_fac, symmfuncs, constraints, Tsymmetry, state_dep, range_prime_map_);
  
  // change to be expression specific
  CTP_map_->insert( new_op->CTP_map()->begin(), new_op->CTP_map()->end());

  return new_op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::create_equation( string name, string type, 
                                        shared_ptr<map<string, shared_ptr<vector<shared_ptr<BraKet_Base>>>>>  term_braket_map,
                                        shared_ptr<map<string, shared_ptr<vector<pair<DataType,string>>>>> expression_term_map, 
                                        shared_ptr<map<pair< string, vector<pair<string, int>>>, 
                                                                            shared_ptr<vector<shared_ptr<BraKet_Base>>>>> term_braket_map_state_spec, 
                                        shared_ptr<map< pair<string, vector<pair<string, int>>>, 
                                                                  shared_ptr<vector<pair<DataType, string>>>>> expression_term_map_state_spec ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __DEBUG_PROPTOOL_SYSTEM_INFO
cout << "System_Info<DataType>::System_Info::create_equation" << endl;
cout << "name = " << name  <<  endl;  cout << "type = " << type  <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  term_braket_map_->insert(term_braket_map->begin(), term_braket_map->end()) ;
  expression_term_map_->insert(expression_term_map->begin(), expression_term_map->end()) ;
  
  shared_ptr<Equation_Base<DataType>>  new_eqn;
  if ( type == "Value" ) { 

    auto new_eqn_val = make_shared<Equation_Value<DataType>>( name, type, states_info_, term_braket_map, expression_term_map );
    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>( new_eqn_val );

    if (!new_eqn) { throw runtime_error("cast from Equation_Value to equation_base failed" ); }

    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map, MT_map_, CTP_map_, range_prime_map_ );
    equation_map_->emplace( name, new_eqn ); 

  } else if ( type == "LinearRM") { 

    shared_ptr<Equation_LinearRM<DataType>> new_eqn_lrm  =
    make_shared<Equation_LinearRM<DataType>>( name, type, states_info_, term_braket_map, expression_term_map,
                                              term_braket_map_state_spec, expression_term_map_state_spec );

    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>(new_eqn_lrm);

    if (!new_eqn) { throw runtime_error("cast from Equation_LinearRM to Equation_Base failed" ); }

    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map,  MT_map_, CTP_map_, range_prime_map_ );
    new_eqn_lrm->generate_state_specific_terms();
    equation_map_->emplace( name, new_eqn); 
  
  } else {  
    throw logic_error( "equation type \""+ type + "\" not implemented yet! Aborting!"); 
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::create_equation( string name, string type, 
                                        shared_ptr<map<string, shared_ptr<vector<shared_ptr<BraKet_Base>>>>>  term_braket_map,
                                        shared_ptr<map<string, shared_ptr<vector<pair<DataType,string>>>>> expression_term_map ){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __DEBUG_PROPTOOL_SYSTEM_INFO
cout << "System_Info<DataType>::System_Info::create_equation : "; cout.flush();
cout << "name = " << name  <<  endl; cout << "type = " << type  <<  endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  term_braket_map_->insert(term_braket_map->begin(), term_braket_map->end()) ;
  expression_term_map_->insert(expression_term_map->begin(), expression_term_map->end()) ;
 
  shared_ptr<Equation_Base<DataType>>  new_eqn;
  if ( type == "Value" ) { 
    shared_ptr<Equation_Value<DataType>> new_eqn_val  = make_shared<Equation_Value<DataType>> ( name, type, states_info_,  term_braket_map, expression_term_map );
    new_eqn = dynamic_pointer_cast<Equation_Base<DataType>>(new_eqn_val);
    new_eqn->set_maps( expression_map, Gamma_map, ACompute_map, MT_map_, CTP_map_, range_prime_map_ );
    equation_map_->emplace( name, new_eqn); 

  } else if ( type == "LinearRM") { 
    throw logic_error( "Must provide state specific term map for doing linearRM!! Aborting!! " ) ; 
    
  } else {  
    throw logic_error( "equation type \""+ type + "\" not implemented yet! Aborting!"); 

  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
