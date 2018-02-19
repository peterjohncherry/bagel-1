#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h>

using namespace std;

//TODO gamma_info_map should be generated outside and not fed in here. It constains _no_ Expression specific information.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
Expression<DataType>::Expression( shared_ptr<vector< BraKet<DataType>>> braket_list,
                                  shared_ptr<StatesInfo<DataType>> states_info,
                                  shared_ptr<map< string, shared_ptr<MultiTensOp::MultiTensOp<DataType>>>>  MT_map,
                                  shared_ptr<map< string, shared_ptr<CtrTensorPart_Base> >>            CTP_map,
                                  shared_ptr<map< string, shared_ptr<vector<shared_ptr<CtrOp_base>> >>>     ACompute_map,
                                  shared_ptr<map< string, shared_ptr<GammaInfo> > >                         gamma_info_map ):
                                  braket_list_(braket_list), states_info_(states_info), MT_map_(MT_map), CTP_map_(CTP_map),
                                  ACompute_map_(ACompute_map), gamma_info_map_(gamma_info_map) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression<DataType>::Expression (new constructor) " << endl;

  //Note that this G_to_A_map_ is expression specific
  G_to_A_map_ = make_shared<map< string, shared_ptr< map<string, AContribInfo >>>>();

  name_ = "";
  for ( BraKet<DataType>& bk : *braket_list_ ) {
    if (bk.factor() != 0.0 )
      name_ += "(" + to_string(bk.factor()) + ")" + bk.name() + " + ";
  }
  name_.pop_back();
  name_.pop_back();

  // Will loop through terms and then generate mathematical task map. It's split into
  // two functions as this will gives more control over merging together of different BraKets G_to_A_maps.
  for ( BraKet<DataType>& braket : *braket_list_ )
    braket.generate_gamma_Atensor_contractions( MT_map_, G_to_A_map_, gamma_info_map_, states_info_ );
  get_gamma_Atensor_contraction_list();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma (as determined by BraKet) into the map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression<DataType>::get_gamma_Atensor_contraction_list(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression::get_gamma_Atensor_contraction_list" << endl;

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map_->begin(); G2A_mapit != G_to_A_map_->end(); G2A_mapit++) {
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      string cmtp_name  = A_map_it->first;
      shared_ptr<vector<shared_ptr<CtrOp_base>>>  ACompute_list;
      if ( CTP_map_->find(cmtp_name) == CTP_map_->end())
        throw std::logic_error( cmtp_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;

      auto ACompute_list_loc = ACompute_map_->find(cmtp_name);
      if ( ACompute_list_loc != ACompute_map_->end() ){
        cout << "Expression::get_gamma_Atensor_contraction_list::already built compute list for " << cmtp_name << " during generation of earlier compute list" << endl;
        cout << cmtp_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map_->at(cmtp_name)->size() << "  --- Still in if " << endl;
        continue;
      } else {
        ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >>(0);
        CTP_map_->at(cmtp_name)->FullContract(CTP_map_, ACompute_list, ACompute_map_);
        ACompute_map_->emplace(cmtp_name, ACompute_list);
        CTP_map_->at(cmtp_name)->got_compute_list( true );
      }
      cout << cmtp_name << " has a compute list of length : "; cout.flush() ; cout << ACompute_map_->at(cmtp_name)->size() << endl;
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
      string cmtp_name  = A_map_it->first;
      auto CMTP_loc = CTP_map_->find(cmtp_name);
      if ( CMTP_loc == CTP_map_->end())
        throw std::logic_error( cmtp_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;

      //awkward, but this will avoid issues where members of CTP_vec are formed from multiple tensors
      shared_ptr<CtrTensorPart_Base> CMTP = CMTP_loc->second;
      shared_ptr<vector<string>> ranges = CMTP->full_id_ranges();
      shared_ptr<vector<string>> idxs = CMTP->full_idxs();
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

//      shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> CTP_vec = CMTP->CTP_vec;
    }
  }
  cout << "leaving Expression::necessary_tensor_blocks" << endl;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
