#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/expression_full.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma (as determined by BraKet) into the map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Full<DataType>::generate_algebraic_task_list(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void Expression_Full<DataType>::generate_algebraic_task_list()" << endl;
    
 
  // Will loop through terms and then generate mathematical task map. It's split into
  // two functions as this will gives more control over merging together of different BraKets' G_to_A_maps.
  required_blocks_ = make_shared<set<shared_ptr<Range_Block_Info>>>();
  for ( shared_ptr<BraKet_Base>& braket : *braket_list_ )
    braket->generate_gamma_Atensor_contractions( MT_map_, G_to_A_map_, gamma_info_map_, states_info_, required_blocks_, CTP_map_ );
 
  this->get_gamma_Atensor_contraction_list();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Full<DataType>::get_gamma_Atensor_contraction_list(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Full::get_gamma_Atensor_contraction_list" << endl;

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map_->begin(); G2A_mapit != G_to_A_map_->end(); G2A_mapit++) {
    auto A_map = G2A_mapit->second;
    for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
      string cmtp_name  = A_map_it->first;
      shared_ptr<vector<shared_ptr<CtrOp_base>>>  ACompute_list;
      if ( CTP_map_->find(cmtp_name) == CTP_map_->end()){
        throw std::logic_error( cmtp_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;
      }
      cout << "cmtp_name = " << cmtp_name << endl;
      auto ACompute_list_loc = ACompute_map_->find(cmtp_name);
      if ( ACompute_list_loc != ACompute_map_->end() ){
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Full<double>;
template class Expression_Full<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
