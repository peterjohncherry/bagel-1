#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/expression.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
Expression<DType>::Expression( std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DType>>>>>> BraKet_list,
                           std::shared_ptr<StatesInfo<DType>> TargetStates_in ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map                 = make_shared<map< string, shared_ptr<TensOp::TensOp<DType>>>>();
  CTP_map               = make_shared<map< string, shared_ptr<CtrTensorPart<DType>> >>();    
  CMTP_map              = make_shared<map< string, shared_ptr<CtrMultiTensorPart<DType>> >>(); 
  ACompute_map          = make_shared<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >>(); 
  CMTP_Eqn_Compute_List = make_shared<map< vector<string>, shared_ptr<vector<pair<shared_ptr<vector<string>>, pair<int,int> >>> >>();

  GammaMap              = make_shared< map<string, shared_ptr<GammaInfo> > >(); 
  G_to_A_map            = make_shared< map<string, shared_ptr< map<string, AContribInfo  >>>>(); 

  TargetStates = TargetStates_in;
 
  for (auto BraKet_Tensors : *BraKet_list) 
    Build_BraKet( BraKet_Tensors );
  
  for (auto braket : BraKet_Terms)
    Get_CMTP_Compute_Terms();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class DType>
void Expression<DType>::Initialize(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T_map  = make_shared<map< string, shared_ptr<TensOp::TensOp<DType>>>>();
  CTP_map    = make_shared<map< string, shared_ptr<CtrTensorPart<DType>> >>();    
  CMTP_map   = make_shared<map< string, shared_ptr<CtrMultiTensorPart<DType>> >>(); 
  ACompute_map = make_shared<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >>(); 
  CMTP_Eqn_Compute_List = make_shared<map< vector<string>, shared_ptr<vector<pair<shared_ptr<vector<string>>, pair<int,int> >>> >>();
  G_to_A_map   = make_shared< map<string, shared_ptr< map<string, AContribInfo >>>>(); 
  GammaMap = make_shared< map<string, shared_ptr<GammaInfo> > >(); 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void Expression<DType>::Build_BraKet(shared_ptr<vector<shared_ptr<TensOp::TensOp<DType>>>> Tens_vec ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<BraKet<DType>> New_BraKet = make_shared<BraKet<DType>>(G_to_A_map, GammaMap, TargetStates );
  New_BraKet->Sub_Ops = Tens_vec;

  New_BraKet->Build_TotalOp();
  cout << New_BraKet->Total_Op->name();   print_vector(*( New_BraKet->Total_Op->aops()), " aops " ) ; cout <<  endl;
  New_BraKet->Build_Gamma_SpinFree(New_BraKet->Total_Op->aops(), New_BraKet->Total_Op->idxs()); 

  CTP_map->insert( New_BraKet->Total_Op->CTP_map->begin(),  New_BraKet->Total_Op->CTP_map->end());

  CMTP_map->insert(New_BraKet->Total_Op->CMTP_map->begin(), New_BraKet->Total_Op->CMTP_map->end());

  BraKet_Terms.push_back(New_BraKet);   
  cout << "pushed back into braket terms" << endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma into the map
// Note this 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void Expression<DType>::Get_CMTP_Compute_Terms(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression::Get_CMTP_Compute_Terms" << endl;  

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto  G2A_mapit = G_to_A_map->begin(); G2A_mapit != G_to_A_map->end(); G2A_mapit++) {
    
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
        CMTP_map->at(CMTP_name)->got_compute_list =true; 
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
