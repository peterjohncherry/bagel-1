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

  T_map          = make_shared< map <string, shared_ptr<TensOp::TensOp<DataType>>>>();
  MT_map         = make_shared< map <string, shared_ptr<MultiTensOp::MultiTensOp<DataType>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
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

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
System_Info<DataType>::construct_equation_task_list( string equation_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Info<DataType>::System_Info::construct_equation_task_list : " << equation_name << endl;

  equation_map_->at( equation_name)->generate_all_expressions(); 
//  if ( eqn->type() == "Value" ) 
//    for ( auto& expression_info : *eqn->expression_term_map_ ) 
//      Build_Expression ( expression_info.first );


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
  CTP_map->insert( new_op->CTP_map()->begin(), new_op->CTP_map()->end());
  cout << "got ctr tens ranges for new_op : " << op_name << endl;

  return new_op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void System_Info<DataType>::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string BraKet_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "System_Info::System_Info::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string term_name ) " << endl;
  
  shared_ptr<vector<shared_ptr<TensOp::TensOp<double>>>> BraKet_Ops = make_shared<vector< shared_ptr<TensOp::TensOp<double>> > >( Op_names->size());

  vector<shared_ptr<TensOp::TensOp<double>>>::iterator BraKet_Ops_it = BraKet_Ops->begin();
  for ( string name : *Op_names ){  
     cout << "looking_for " << name << " ... " ; cout.flush();
    *BraKet_Ops_it++ = T_map->at(name);
     cout << " found it! " << endl;    
  }
  braket_map_->emplace(BraKet_name, BraKet_Ops);

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
    equation_map_->emplace( name, new_eqn); 

  } else { 
    throw logic_error( "equation type \""+ type + "\" not implemented yet! Aborting!"); 

  }

  return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
