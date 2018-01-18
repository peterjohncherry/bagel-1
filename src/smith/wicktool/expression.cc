#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/expression.h>

using namespace std;

//TODO GammaMap should be generated outside and not fed in here. It constains _no_ Expression specific information.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression<DataType>::Expression( vector< BraKet<DataType>>&  Term_list,
                                  shared_ptr<StatesInfo<DataType>> target_states,
                                  shared_ptr<map< string, shared_ptr<MultiTensOp::MultiTensOp<DataType>>>>  MT_map,      
                                  shared_ptr<map< string, shared_ptr<CtrTensorPart<DataType>> >>            CTP_map,     
                                  shared_ptr<map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>       CMTP_map,    
                                  shared_ptr<map< string, shared_ptr<vector<shared_ptr<CtrOp_base>> >>>     ACompute_map,
                                  shared_ptr<map< string, shared_ptr<GammaInfo> > >                         GammaMap ):    
                                  Term_list_(Term_list), target_states_(target_states), MT_map_(MT_map), CTP_map_(CTP_map), CMTP_map_(CMTP_map),
                                  ACompute_map_(ACompute_map), GammaMap_(GammaMap) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression<DataType>::Expression (new constructor) " << endl; 

  G_to_A_map_ = make_shared<map< string, shared_ptr< map<string, AContribInfo >>>>();

  // Will loop through terms and then generate compute list for each. I've split it up like this
  // for  merging BraKet with BraKet.
  GammaMap_ = GammaMap;
  for ( BraKet<DataType> braket : Term_list_ ) 
    braket.generate_gamma_Atensor_contractions( MT_map_, G_to_A_map_, GammaMap_, target_states_ ); //TODO replace with init gammas func
  
  Get_CMTP_Compute_Terms();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma (as determined by BraKet) into the map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression<DataType>::Get_CMTP_Compute_Terms(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression::Get_CMTP_Compute_Terms" << endl;  

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map_->begin(); G2A_mapit != G_to_A_map_->end(); G2A_mapit++) {
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      string CMTP_name  = A_map_it->first;     
      shared_ptr<vector<shared_ptr<CtrOp_base>>>  ACompute_list; 
      if ( CMTP_map_->find(CMTP_name) == CMTP_map_->end())
        throw std::logic_error( CMTP_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;

      auto ACompute_list_loc = ACompute_map_->find(CMTP_name);
      if ( ACompute_list_loc != ACompute_map_->end() ){
        cout << "Expression::Get_CMTP_Compute_Terms::already built compute list for " << CMTP_name << " during generation of earlier compute list" << endl;
        cout << CMTP_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map_->at(CMTP_name)->size() << "  --- Still in if " << endl;
        continue;
      } else {  
        ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >>(0);
        CMTP_map_->at(CMTP_name)->FullContract(CTP_map_, ACompute_list, ACompute_map_);
        ACompute_map_->emplace(CMTP_name, ACompute_list);
        CMTP_map_->at(CMTP_name)->got_compute_list = true; 
      }
      cout << CMTP_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map_->at(CMTP_name)->size() << endl;
    }
  }
  cout << "leaving Get_CMTP_compute_Terms" << endl;
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
  for (auto G2A_mapit = G_to_A_map_->begin(); G2A_mapit != G_to_A_map_->end(); G2A_mapit++) {
    
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      string CMTP_name  = A_map_it->first;
      auto CMTP_loc = CMTP_map_->find(CMTP_name);
      if ( CMTP_loc == CMTP_map_->end())
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
