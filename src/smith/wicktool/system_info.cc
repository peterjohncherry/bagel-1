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
  BraKet_map     = make_shared<  map < string, shared_ptr<vector<shared_ptr< TensOp::TensOp<DataType>>>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
  expression_map = make_shared< map <string, shared_ptr<Expression<DataType>>>>();

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
                                                  DataType factor, 
                                                  string Tsymmetry,
                                                  bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "System_Info<DataType>::System_Info::Build_TensOp" <<   endl;

  //NOTE: change to use proper factor
  int tmpfac = 1;
  shared_ptr<TensOp::TensOp<DataType>>  New_Op = make_shared<TensOp::TensOp<DataType>>( op_name, *op_idxs, *op_idx_ranges, *op_aops,
                                                                                        tmpfac,  Symmetry_Funcs, Constraint_Funcs, Tsymmetry, target_states_);
  // change to be state specific
  // New_Op->get_ctrs_tens_ranges();

  return New_Op;
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
string System_Info<DataType>::System_Info::Build_Expression( vector<Term_Info<DataType>>& term_info_list  ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO Argument should stay vector of term_info, but current term_info is really braket_info,
//     change so each element of vec is term_info containing range of states.

cout << "System_Info::System_Info::Build_Expression" << endl;                                                                                     
  shared_ptr< vector<pair<string, DataType>> > BraKet_name_list = make_shared<vector<pair< string, DataType >>>(0);
  
  // This is looping over states; op sparsity depends on states, should replace with term_info_map, and
  // have double loop, outer for ket state, inner for brastate
  for ( Term_Info<DataType> BraKet_info : term_info_list ) {


    string MT_name = "";
    for (string op_name : BraKet_info.op_list ) { // TODO should loop  over states defined in term_info, 

      if( T_map->find(op_name) == T_map->end() ) 
        Initialize_Tensor_Op_Info( op_name );
      MT_name += op_name;

      //TODO do state specific definition; just builds new range_map, should not have any sparsity yet ...
      
    }

    if( MT_map->find( MT_name ) == MT_map->end() ){

      vector<shared_ptr<TensOp::TensOp<DataType>>> SubOps(BraKet_info.op_list.size());
      for (int ii = 0 ; ii != BraKet_info.op_list.size() ; ii++ )  
        SubOps[ii] = T_map->at(BraKet_info.op_list[ii]); 

      shared_ptr<MultiTensOp::MultiTensOp<DataType>> MultiOp = make_shared<MultiTensOp::MultiTensOp<DataType>>( MT_name, spinfree_, SubOps, target_states_ );
      MT_map->emplace(MT_name , MultiOp );
    } 
   
    //TODO do state specific definition for MultiTens Op; just builds new range_map, should not have any sparsity yet ...
      
    //TODO requires state specific looping to get name
    string BraKet_name =  Get_BraKet_name( BraKet_info  ); 

    if ( BraKet_map->find(BraKet_name) == BraKet_map->end() ) 
      Set_BraKet_Ops( make_shared<vector<string>>(BraKet_info.op_list), BraKet_name ) ;

    BraKet_name_list->push_back( make_pair( BraKet_name, BraKet_info.factor ) );
    cout << "new BraKet_name = " << BraKet_name << endl;
  }

  string expression_name = "";
  for ( pair<string, DataType> name_fac_pair : *BraKet_name_list ) {
    if (name_fac_pair.second != 0.0 ) 
      expression_name += "+ (" + to_string(name_fac_pair.second) + ")" + name_fac_pair.first;
  }
  cout << "new_expression_name = " << expression_name << endl;
    
  shared_ptr<Expression<DataType>> new_expression = make_shared<Expression<DataType>>(BraKet_name_list, BraKet_map, target_states_);
  
  expression_map->emplace( expression_name, new_expression );

  return expression_name;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string System_Info<DataType>::System_Info::Get_BraKet_name( Term_Info<DataType>& BraKet_info  ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " System_Info<DataType>::System_Info::Get_BraKet_name " << endl;

  string BraKet_name = "";

  if ( BraKet_info.type == "ci_derivative" ) {
    BraKet_name += "c_{I}^{" + BraKet_info.Bra_name + "} < "+BraKet_info.Bra_name  + " | ";
  } else if (BraKet_info.type == "expectation")  { 
    BraKet_name += "< "+ BraKet_info.Bra_name +" | ";
  }
  
  for ( string op : BraKet_info.op_list ) 
    BraKet_name += op;
   
  BraKet_name += " | " + BraKet_info.Ket_name + " > " ;
  
  return BraKet_name;  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
