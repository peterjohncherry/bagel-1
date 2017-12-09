#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include  <src/smith/wicktool/system_info.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Note, spinfree should tell us not just if the wavefunction is free, but whether or 
//not the perturbation being applied is spin independent
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
System_Info<DataType>::System_Info::System_Info( shared_ptr<StatesInfo<DataType>> TargetStates_in , bool spinfree ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map          = make_shared< map <string, shared_ptr<TensOp<DataType>>>>();
  BraKet_map     = make_shared<  std::map <std::string, 
                                          std::shared_ptr<std::vector<std::shared_ptr< TensOp<DataType>>>>>>();
  CTP_map        = make_shared< map <string, shared_ptr<CtrTensorPart<DataType>>>>();
  CMTP_map       = make_shared< map <string, shared_ptr<CtrMultiTensorPart<DataType>>>>();
  expression_map = make_shared< map <string, shared_ptr<Expression<DataType>>>>();

  TargetStates = TargetStates_in;

  spinfree_ = spinfree;

  if (spinfree || !spinfree /*TODO like this for testing; obviously must put back*/ ) { 
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
shared_ptr<TensOp<DataType>>
System_Info<DataType>::System_Info::Build_TensOp( string op_name,
                                                  shared_ptr<DataType> tensor_data, 
                                                  shared_ptr<vector<string>> op_idxs,
                                                  shared_ptr<vector<bool>> op_aops, 
                                                  shared_ptr<vector<vector<string>>> op_idx_ranges,
                                                  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                                                  vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                                                  pair<double,double> factor, 
                                                  string Tsymmetry,
                                                  bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "System_Info<DataType>::System_Info::Build_TensOp" <<   endl;

  //NOTE: change to use proper factor
  int tmpfac = 1;
  shared_ptr<TensOp<DataType>>  New_Op = make_shared<TensOp<DataType>>(op_name, Symmetry_Funcs, Constraint_Funcs,
                                                                       *op_idxs, *op_idx_ranges, *op_aops, tmpfac, Tsymmetry);
  New_Op->get_ctrs_tens_ranges();

  return New_Op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void System_Info<DataType>::System_Info::Set_BraKet_Ops(shared_ptr<vector<string>> Op_names, string term_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout <<  "System_Info::System_Info::Build_BraKet(shared_ptr<vector<string>> BraKet_names, string expression_name ) " << endl;
  
  shared_ptr<vector<shared_ptr<TensOp<double>>>> BraKet_Ops = make_shared<vector< shared_ptr<TensOp<double>> > >();

  for ( string name : *Op_names ) 
    BraKet_Ops->push_back(T_map->at(name));

  BraKet_map->emplace(term_name, BraKet_Ops);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
string System_Info<DataType>::System_Info::Build_Expression( vector<Term_Info<DataType>>&  term_info_list  ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "System_Info::System_Info::Build_Expression" << endl;                                                                                     
  
  shared_ptr<vector<string>> BraKet_name_list = make_shared<vector<string>>(0);
  
  for ( Term_Info<DataType> BraKet_info : term_info_list ) {

    for (string op_name : BraKet_info.op_list ) {
      bool get_op = true;
      for ( string defined_op : Op_list ) {
        if ( op_name == defined_op ) {
          cout << " already defined " << op_name << " will not redefine" << endl;
          get_op =false; 
          break;
        }
      }
      if( get_op ) {
        Initialize_Tensor_Op_Info( op_name );
        Op_list.push_back(op_name);
      }
    }

    string BraKet_name = "";
    
    if ( BraKet_info.type == "ci_derivative" ) {
      BraKet_name += "c_{I}^{" + BraKet_info.Bra_name + "} < "+BraKet_info.Bra_name  + " | ";
    } else { 
      BraKet_name += "< "+ BraKet_info.Bra_name +" | ";
    }
    
    for ( string op : BraKet_info.op_list ) 
      BraKet_name += op;
    
    BraKet_name += " | " + BraKet_info.Ket_name + " > " ;
   
    Set_BraKet_Ops( make_shared<vector<string>>(BraKet_info.op_list), BraKet_name ) ;

    BraKet_name_list->push_back(BraKet_name); 
    cout << "new BraKet_name = " << BraKet_name << endl;

  }    
 
  string expression_name = accumulate(BraKet_name_list->begin(), BraKet_name_list->end(), string(""));
  cout << "new_expression_name = " << expression_name << endl;

  auto BraKet_List = make_shared<std::vector<std::shared_ptr<std::vector<std::shared_ptr<TensOp<DataType>>>>>>( BraKet_name_list->size() ); 

  for ( int ii = 0 ; ii != BraKet_name_list->size() ; ii++ ) 
    BraKet_List->at(ii) = BraKet_map->at(BraKet_name_list->at(ii)) ;

  shared_ptr<Expression<DataType>> new_expression = make_shared<Expression<DataType>>(BraKet_List, TargetStates);

  expression_map->emplace(expression_name, new_expression);

  return expression_name;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Info<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
