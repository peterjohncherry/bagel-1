#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/expression.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression<DataType>::Expression( shared_ptr<vector< shared_ptr<Term_Info<DataType>>>>                      Term_list,
                                  shared_ptr<map< string, shared_ptr<MultiTensOp::MultiTensOp<DataType>>>>  MT_map,      
                                  shared_ptr<map< string, shared_ptr<CtrTensorPart<DataType>> >>            CTP_map,     
                                  shared_ptr<map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>       CMTP_map,    
                                  shared_ptr<map< string, shared_ptr<vector<shared_ptr<CtrOp_base>> >>>     ACompute_map,
                                  shared_ptr<map< string, shared_ptr<GammaInfo> > >                         GammaMap ):    
                                  Term_list_(Term_list), MT_map_(MT_map), CTP_map_(CTP_map), CMTP_map_(CMTP_map),
                                  ACompute_map_(ACompute_map), GammaMap_(GammaMap){
                                  //TODO GammaMap should be generated outside and not fed in here. It constains _no_ Expression specific information.
                                  //   G_to_A_map does, and it is G_to_A_map which should be used for looping in Expression_Computer.
                                  //   Doing it this way will result in a larger GammaMap, but there will only be one of them, the only issue is how large?
                                  //   Potentially can generate all expressions, and then generate elements of gamma map as need be from G_to_A_map
                                  //   info.
                                  //   ACompute_map contains information only used by this expression, but I do not think it is expression specific.
                                  //   However, state information must be incorporated.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression<DataType>::Expression (new constructor) " << endl; 
  G_to_A_map_ = make_shared<map< string, shared_ptr< map<string, AContribInfo >>>>();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression<DataType>::Expression( shared_ptr< vector< pair< string, DataType >>> BraKet_list,
                                  shared_ptr< map < string, shared_ptr< vector < shared_ptr< TensOp::TensOp<DataType> >>>>> BraKet_map,
                                  shared_ptr< StatesInfo<DataType> > target_states ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map        = make_shared<map< string, shared_ptr<TensOp::TensOp<DataType>>>>();
  CTP_map      = make_shared<map< string, shared_ptr<CtrTensorPart<DataType>> >>();    
  CMTP_map     = make_shared<map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>(); 
  ACompute_map = make_shared<map< string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >>(); 

  GammaMap     = make_shared< map< string, shared_ptr<GammaInfo> > >(); 
  G_to_A_map   = make_shared< map< string, shared_ptr< map<string, AContribInfo  >>>>(); 

  target_states_ = target_states;
  
  for ( pair<string,DataType> bk_name_factor :  *BraKet_list )
    Build_BraKet( BraKet_map->at( bk_name_factor.first), bk_name_factor.second );
  
  for (auto braket : BraKet_Terms)
    Get_CMTP_Compute_Terms();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class DataType>
void Expression<DataType>::Initialize(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map  = make_shared<map< string, shared_ptr<TensOp::TensOp<DataType>>>>();
  CTP_map    = make_shared<map< string, shared_ptr<CtrTensorPart<DataType>> >>();    
  CMTP_map   = make_shared<map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>(); 
  ACompute_map = make_shared<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >>(); 
  G_to_A_map   = make_shared< map<string, shared_ptr< map<string, AContribInfo >>>>(); 
  GammaMap = make_shared< map<string, shared_ptr<GammaInfo> > >(); 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression<DataType>::Build_BraKet( shared_ptr<vector<shared_ptr<TensOp::TensOp<DataType>>>> Tens_vec, DataType factor  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //TODO Feed MultiTensOp in as arg, no need for SubOps/Build_TotalOp() 
  shared_ptr<BraKet<DataType>> New_BraKet = make_shared<BraKet<DataType>>(G_to_A_map, GammaMap, target_states_, factor );
  New_BraKet->Sub_Ops = Tens_vec;

  New_BraKet->Build_TotalOp();
  cout << New_BraKet->Total_Op->name();   print_vector(*( New_BraKet->Total_Op->aops()), " aops " ) ; cout << endl;
  New_BraKet->Build_Gamma_SpinFree(New_BraKet->Total_Op->aops(), New_BraKet->Total_Op->idxs()); 

  CTP_map->insert( New_BraKet->Total_Op->CTP_map->begin(),  New_BraKet->Total_Op->CTP_map->end());
  CMTP_map->insert(New_BraKet->Total_Op->CMTP_map->begin(), New_BraKet->Total_Op->CMTP_map->end());

  BraKet_Terms.push_back(New_BraKet);   
  cout << "pushed back into braket terms" << endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Builds a list of the tensor blocks which are needed for evaluation of this expression.
// This effectively occurs in other functions, but is useful to have seperately here for allocation purposes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression<DataType>::necessary_tensor_blocks(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression::necessary_tensor_blocks" << endl;  

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map->begin(); G2A_mapit != G_to_A_map->end(); G2A_mapit++) {
    
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      string CMTP_name  = A_map_it->first;
      auto CMTP_loc = CMTP_map->find(CMTP_name);
      if ( CMTP_loc == CMTP_map->end())
        throw std::logic_error( CMTP_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;

      //awkward, but this will avoid issues where members of CTP_vec are formed from multiple tensors
      shared_ptr<CtrMultiTensorPart<DataType>> CMTP = CMTP_loc->second;     
      shared_ptr<vector<string>> ranges = CMTP->full_id_ranges;
      shared_ptr<vector<string>> idxs = CMTP->full_idxs;
      int t_num = 0;
      int op_size = 1;
      char op_name = idxs->at(0)[0];
      for ( int ii = 1 ; ii != idxs->size() ; ii++ ) {
        if ( op_name == idxs->at(ii)[0] ){ 
          op_size+=1; 
        } else { 
          vector<string> ctp_idxs(op_size);
          vector<string> ctp_ranges(op_size);
          int pos = 0;
          for ( int jj = ii - op_size; jj != ii ; jj++, pos++ ) {
            ctp_idxs[pos] = idxs->at(pos);
            ctp_ranges[pos] = ranges->at(pos);
          }
          op_name = idxs->at(ii)[0]; 
          op_size = 1; 
        }
      }
 
//      shared_ptr<vector<shared_ptr<CtrTensorPart<DataType>>>> CTP_vec = CMTP->CTP_vec;
      

    }
    cout << "X" << endl;
  }
  cout << "leaving Get_CMTP_compute_Terms" << endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma into the map
// Note this 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression<DataType>::Get_CMTP_Compute_Terms(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression::Get_CMTP_Compute_Terms" << endl;  

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map->begin(); G2A_mapit != G_to_A_map->end(); G2A_mapit++) {
    
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      
      string CMTP_name  = A_map_it->first;

      shared_ptr<vector<shared_ptr<CtrOp_base>>>  ACompute_list; 
      if ( CMTP_map->find(CMTP_name) == CMTP_map->end())
        throw std::logic_error( CMTP_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;

      auto ACompute_list_loc = ACompute_map->find(CMTP_name);
      if ( ACompute_list_loc != ACompute_map->end() ){
        cout << "Expression::Get_CMTP_Compute_Terms::already built compute list for " << CMTP_name << " during generation of earlier compute list" << endl;
        cout << CMTP_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map->at(CMTP_name)->size() << "  --- Still in if " << endl;
        continue;
      } else {  
        ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >>(0);
        CMTP_map->at(CMTP_name)->FullContract(CTP_map, ACompute_list, ACompute_map);
        ACompute_map->emplace(CMTP_name, ACompute_list);
        CMTP_map->at(CMTP_name)->got_compute_list = true; 
      }
  
      cout << CMTP_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map->at(CMTP_name)->size() << endl;

    }
    cout << "X" << endl;
  }
  cout << "leaving Get_CMTP_compute_Terms" << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
