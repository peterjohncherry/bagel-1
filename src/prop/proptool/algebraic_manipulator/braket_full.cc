#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket_full.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace WickUtils;

#define __DEBUG_PROPTOOL_BRAKET_FULL
#define __DEBUG_PROPTOOL_BRAKET_FULL_VERBOSE

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet_Full<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<TensOp_Base>>> MT_map,
                                                                 shared_ptr<map<string, shared_ptr<map<string,shared_ptr<AContribInfo_Base> >>>> G_to_A_map,
                                                                 shared_ptr<map<string, shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                                                 shared_ptr<StatesInfo_Base> target_states,
                                                                 shared_ptr<set<shared_ptr<Range_Block_Info>>> required_blocks,
                                                                 shared_ptr<map<string, shared_ptr<CtrTensorPart_Base>>> ctp_map  ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_BRAKET_FULL
cout << "BraKet_Full::generate_gamma_Atensor_contractions : " << name_ << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Total_Op_ = MT_map->at( multiop_info_->op_name_ );

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( *(multiop_info_->op_order_),  *(multiop_info_->transformations_) );
  Total_Op_->generate_ranges( multiop_info_ );

  factor_.first = 1.0; 
  factor_.second = 0.0; 
  shared_ptr<GammaGeneratorRedux<DataType>> GGen = make_shared<GammaGeneratorRedux<DataType>>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, G_to_A_map, factor_ );

  auto all_ranges = Total_Op_->all_ranges_state_specific_->at( multiop_info_->name_ );

  if ( multiop_info_->op_info_vec()->size() > 1 ) { 
 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){

      #ifdef __DEBUG_PROPTOOL_BRAKET_FULL_VERBOSE
        print_braket_block( multiop_info_, range_map_it->second ); 
      #endif
  
      // TODO if is only here for non-relativistic case ; constraints should be specified in input
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
  
        GGen->add_gamma( range_map_it->second, trans_aops );

        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          cout << "done ANO" << endl;
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            cout << "done NO" << endl;
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
              cout << "done AO" << endl;
             
              for (  auto& block :  *(range_map_it->second->unique_block_->range_blocks()) ){
                MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
                required_blocks->emplace( block );
                cout << "require_block =  " <<  block->name() << endl;
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
          cout << "done ANO" << endl;
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            cout << "done NO" << endl;
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
              cout << "done AO" << endl;

              auto& block = range_map_it->second->unique_block_; 
              MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
              required_blocks->emplace( block );
              cout << "require_block =  " <<  block->name() << endl;
              
            }
          }
        }
      }
    }
  }
  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );
  
#ifdef __DEBUG_PROPTOOL_BRAKET_FULL_VERBOSE
  print_gamma_Atensor_contractions( G_to_A_map, false );
#endif
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet_Full<double>;
template class BraKet_Full<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
