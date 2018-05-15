#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket_full.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet_Full<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<TensOp_Base>>> MT_map,
                                                                 shared_ptr<map<string, shared_ptr<map<string,shared_ptr<AContribInfo_Base> >>>> G_to_A_map,
                                                                 shared_ptr<map<string, shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                                                 shared_ptr<StatesInfo_Base> target_states,
                                                                 shared_ptr<set<shared_ptr<Range_Block_Info>>> required_blocks,
                                                                 shared_ptr<map<string, shared_ptr<CtrTensorPart_Base>>> ctp_map  ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet_Full::generate_gamma_Atensor_contractions : " << name_ << endl;

  Total_Op_ = MT_map->at( multiop_info_->op_name_ );

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( *(multiop_info_->op_order_),  *(multiop_info_->transformations_) );

  Total_Op_->generate_ranges( multiop_info_ );

  shared_ptr<GammaGeneratorRedux<DataType>> GGen = make_shared<GammaGeneratorRedux<DataType>>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, G_to_A_map, factor_ );

  auto all_ranges = Total_Op_->all_ranges_state_specific_->at( multiop_info_->name_ );

  if ( multiop_info_->op_info_vec()->size() > 1 ) { 
 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
  
      // TODO if is only here for non-relativistic case ; constraints should be specified in input
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
  
        GGen->add_gamma( range_map_it->second, trans_aops );

        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){

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
              for (  auto& block :  *(range_map_it->second->unique_block_->range_blocks()) ){
                MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
                required_blocks->emplace( block );
              }
            }
          }
        }
      }
    }
  }

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

  print_gamma_Atensor_contractions( G_to_A_map, false );

  return; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet_Full<double>;
template class BraKet_Full<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
