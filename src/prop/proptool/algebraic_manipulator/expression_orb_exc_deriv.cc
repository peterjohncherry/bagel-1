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
 
  required_blocks_ = make_shared<set<shared_ptr<Range_Block_Info>>>();

  cout << " braket_list_->front().target_op_ = " <<  braket_list_->front()->target_op() << endl;

  //auto  exc_block_G_to_A_map = make_shared<map< string, shared_ptr< map<string, shared_ptr<AContribInfo_Base>> >>>();
  exc_block_G_to_A_map_ = make_shared<map<string, shared_ptr<map< string, shared_ptr< map<string, shared_ptr<AContribInfo_Base>> >>>>>();
  for ( shared_ptr<BraKet_Base>& braket : *braket_list_ )
    braket->generate_gamma_Atensor_contractions( MT_map_, exc_block_G_to_A_map_, gamma_info_map_, states_info_,  required_blocks_, CTP_map_ );

  //TODO This is repeatedly generating the task list for different A; this is not necessary; define seperate vector or map for all required A.
  for ( auto& elem : *exc_block_G_to_A_map_ )
    this->get_gamma_Atensor_contraction_list( elem.second );
  
  cout << "EOED::gatl1 " <<endl;
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
Expression_Orb_Exc_Deriv<DataType>::get_gamma_Atensor_contraction_list( shared_ptr< map< string, shared_ptr< map<string, shared_ptr<AContribInfo_Base>> >>>& G_to_A_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Orb_Exc_Deriv::get_gamma_Atensor_contraction_list" << endl;

  //loop through G_to_A_map ; get all A-tensors associated with a given gamma
  for (auto G2A_mapit = G_to_A_map->begin(); G2A_mapit != G_to_A_map->end(); G2A_mapit++) {

    auto post_reorder_map = G2A_mapit->second;
    for (auto prm_it = post_reorder_map->begin(); prm_it != post_reorder_map->end(); prm_it++) {
      auto gamma_pos_map = prm_it->second->gamma_pos_map(); 

      for (auto gpm_it = gamma_pos_map->begin(); gpm_it != gamma_pos_map->end(); gpm_it++) {
        auto A_map = gpm_it->second;

        for (auto A_map_it = A_map->begin(); A_map_it != A_map->end(); A_map_it++){
          string cmtp_name  = A_map_it->first;

          if ( CTP_map_->find(cmtp_name) == CTP_map_->end())
            throw std::logic_error( cmtp_name + " is not yet in the map!! Generation of Gamma contributions probably has problems!! " ) ;
       
          auto ACompute_list_loc = ACompute_map_->find(cmtp_name);
          if ( ACompute_list_loc != ACompute_map_->end() ){
            continue;
          } else {
            shared_ptr<vector<shared_ptr<CtrOp_base>>>  ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >>(0);
            CTP_map_->at(cmtp_name)->build_contraction_sequence(CTP_map_, ACompute_list, ACompute_map_);
            ACompute_map_->emplace(cmtp_name, ACompute_list);
            CTP_map_->at(cmtp_name)->got_compute_list( true );
          }
   
        }
      }
    }
  }
  return;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Expression_Orb_Exc_Deriv<double>;
template class Expression_Orb_Exc_Deriv<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
