#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket_orb_exc_deriv.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet_OrbExcDeriv<DataType>::generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string, std::shared_ptr<TensOp_Base>>> MT_map,                
                                                                        std::shared_ptr<std::map<std::string,
                                                                                        std::shared_ptr<std::map<std::string,
                                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> >> block_G_to_A_map,
                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                                                        std::shared_ptr<StatesInfo_Base> target_states,
                                                                        std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks,
                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map ) {  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_BRAKET_ORBEXCDERIV
cout << "BraKet_OrbExcDeriv<DataType>::generate_gamma_Atensor_contractions : " << name_ << " "; cout.flush(); cout << "target_op_ = " <<  target_op_ << endl;
cout << "multiop_info_->name_ = "; cout.flush(); cout <<  multiop_info_->name_ << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Total_Op_ = MT_map->at( multiop_info_->op_name_ );

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( *(multiop_info_->op_order_),  *(multiop_info_->transformations_) );

  Total_Op_->generate_ranges( multiop_info_ );

  shared_ptr<GammaGenerator_OrbExcDeriv<DataType>> GGen = make_shared<GammaGenerator_OrbExcDeriv<DataType>>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, block_G_to_A_map, factor_, target_op_ );

  auto all_ranges = Total_Op_->all_ranges_state_specific_->at( multiop_info_->name_ );

  if ( multiop_info_->op_info_vec()->size() > 1 ) { 
 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
        GGen->add_gamma( range_map_it->second, trans_aops );
        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
              vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();

              for (  auto& block :  *(range_map_it->second->unique_block_->range_blocks()) ){
                MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
                required_blocks->emplace( block );
              }

            }
          }
        }
      }
    }
  } else { 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
        GGen->add_gamma( range_map_it->second, trans_aops );
        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
                MT_map->at( Total_Op_->name() )->add_required_block( range_map_it->second->unique_block_);
                required_blocks->emplace( range_map_it->second->unique_block_ );
              }
            }
          }
        }
      }
    }

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

  return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet_OrbExcDeriv<double>;
template class BraKet_OrbExcDeriv<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
