#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include  <src/smith/wicktool/system_info.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Note, spinfree should tell us not just if the wavefunction is free, but whether or 
//not the perturbation being applied is spin independent
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
System_Info<DataType>::System_Info::System_Info( shared_ptr<StatesInfo<DataType>> target_states , bool spinfree ) :
                                                 target_states_(target_states), spinfree_(spinfree) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map          = make_shared< map <string, shared_ptr<TensOp::TensOp<DataType>>>>();
  MT_map         = make_shared< map <string, shared_ptr<MultiTensOp::MultiTensOp<DataType>>>>();
  BraKet_map     = make_shared< map <string, shared_ptr<vector<shared_ptr< TensOp::TensOp<DataType>>>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
  expression_map = make_shared< map <string, shared_ptr<Expression<DataType>>>>();
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
shared_ptr<TensOp::TensOp<DataType>>
System_Info<DataType>::System_Info::Build_TensOp( string op_name,
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
  CTP_map->insert( new_op->CTP_map->begin(), new_op->CTP_map->end());
  cout << "got ctr tens ranges for new_op : " << op_name << endl;

  return new_op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void System_Info<DataType>::System_Info::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string BraKet_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "System_Info::System_Info::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string term_name ) " << endl;
  
  shared_ptr<vector<shared_ptr<TensOp::TensOp<double>>>> BraKet_Ops = make_shared<vector< shared_ptr<TensOp::TensOp<double>> > >( Op_names->size());

  vector<shared_ptr<TensOp::TensOp<double>>>::iterator BraKet_Ops_it = BraKet_Ops->begin();
  for ( string name : *Op_names ){  
     cout << "looking_for " << name << " ... " ; cout.flush();
    *BraKet_Ops_it++ = T_map->at(name);
     cout << " found it! " << endl;    
  }
  BraKet_map->emplace(BraKet_name, BraKet_Ops);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string System_Info<DataType>::System_Info::Build_Expression( vector<BraKet<DataType>>& term_info_list  ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "System_Info::System_Info::Build_Expression" << endl;                                                                                     
  shared_ptr< vector<pair<string, DataType>> > BraKet_name_list = make_shared<vector<pair< string, DataType >>>(0);
  
  // This is looping over states; op sparsity depends on states, should replace with term_info_map, and
  // have double loop, outer for ket state, inner for brastate
  for ( BraKet<DataType> BraKet_info : term_info_list ) {

    for (string op_name : BraKet_info.op_list_ ) { // TODO should loop  over states defined in term_info, 
      
      auto T_loc = T_map->find(op_name);
      if( T_loc == T_map->end() ){ 
        shared_ptr<TensOp::TensOp<DataType>> new_op = Initialize_Tensor_Op_Info( op_name );
        T_map->emplace( op_name, new_op );
        CTP_map->insert( new_op->CTP_map->begin(), new_op->CTP_map->end());
      //TODO do state specific definition; just builds new range_map, should not have any sparsity yet ...
      }
    } 

    if( MT_map->find( BraKet_info.multiop_name_ ) == MT_map->end() ){

      vector<shared_ptr<TensOp::TensOp<DataType>>> SubOps(BraKet_info.op_list_.size());
      for (int ii = 0 ; ii != BraKet_info.op_list_.size() ; ii++ )  
        SubOps[ii] = T_map->at(BraKet_info.op_list_[ii]); 

      shared_ptr<MultiTensOp::MultiTensOp<DataType>> multiop = make_shared<MultiTensOp::MultiTensOp<DataType>>( BraKet_info.multiop_name_, spinfree_, SubOps );
      multiop->get_ctrs_tens_ranges();
      CTP_map->insert( multiop->CTP_map->begin(), multiop->CTP_map->end());
      CMTP_map->insert( multiop->CMTP_map->begin(), multiop->CMTP_map->end());
      MT_map->emplace(BraKet_info.multiop_name_, multiop );
    } 
   
    //TODO do state specific definition for MultiTens Op; just builds new range_map, should not have any sparsity yet ...
      
    //TODO requires state specific looping to get name
    string BraKet_name =  Get_BraKet_name( BraKet_info  ); 

    if ( BraKet_map->find(BraKet_name) == BraKet_map->end() ) 
      Set_BraKet_Ops( make_shared<vector<string>>(BraKet_info.op_list_), BraKet_name ) ;

    BraKet_name_list->push_back( make_pair( BraKet_name, BraKet_info.factor_ ) );
    cout << "new BraKet_name = " << BraKet_name << endl;
  }

  string expression_name = "";
  for ( pair<string, DataType> name_fac_pair : *BraKet_name_list ) {
    if (name_fac_pair.second != 0.0 ) 
      expression_name += "+ (" + to_string(name_fac_pair.second) + ")" + name_fac_pair.first;
  }
  cout << "new_expression_name = " << expression_name << endl;
    
  shared_ptr<Expression<DataType>> new_expression = make_shared<Expression<DataType>>( term_info_list, target_states_, MT_map, CTP_map, CMTP_map, ACompute_map, Gamma_map ); 

  expression_map->emplace( expression_name, new_expression );

  return expression_name;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string System_Info<DataType>::System_Info::Get_BraKet_name( BraKet<DataType>& BraKet_info  ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " System_Info<DataType>::System_Info::Get_BraKet_name " << endl;

  string BraKet_name = "";

  if ( BraKet_info.type_ == "ci_derivative" ) {
    BraKet_name += "c_{I}^{" + to_string(BraKet_info.bra_num_) + "} < "+to_string(BraKet_info.bra_num_)  + " | ";
  } else if (BraKet_info.type_ == "expectation" || BraKet_info.type_ == "full")  { 
    BraKet_name += "< "+ to_string(BraKet_info.bra_num_) +" | ";
  }
  
  BraKet_name += BraKet_info.multiop_name_;
   
  BraKet_name += " | " + to_string(BraKet_info.ket_num_) + " > " ;
  
  return BraKet_name;  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
