#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/expression_orb_exc_deriv.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma (as determined by BraKet) into the map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Expression_Orb_Exc_Deriv<DataType>::generate_algebraic_task_list(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void Expression_Orb_Exc_Deriv<DataType>::generate_algebraic_task_list()" << endl;
 
  required_blocks_ = make_shared<set<string>>();

  cout << " braket_list_->front().target_op_ = " <<  braket_list_->front()->target_op() << endl;

  auto  exc_block_G_to_A_map = make_shared<map< string, shared_ptr< map<string, shared_ptr<AContribInfo_Base>> >>>();
  for ( shared_ptr<BraKet_Base>& braket : *braket_list_ )
    braket->generate_gamma_Atensor_contractions( MT_map_, exc_block_G_to_A_map, gamma_info_map_, states_info_,  required_blocks_, CTP_map_ );

  this->get_gamma_Atensor_contraction_list( exc_block_G_to_A_map );
//    target_to_G_to_A_map_->emplace( exc_block_name, exc_block_G_to_A_map);
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
Expression_Orb_Exc_Deriv<DataType>::get_gamma_Atensor_contraction_list( shared_ptr<map< string, shared_ptr< map<string, shared_ptr<AContribInfo_Base>> >>> exc_block_G_to_A_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Orb_Exc_Deriv::get_gamma_Atensor_contraction_list" << endl;

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = exc_block_G_to_A_map->begin(); G2A_mapit != exc_block_G_to_A_map->end(); G2A_mapit++) {
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Orb_Exc_Deriv<double>;
template class Expression_Orb_Exc_Deriv<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
