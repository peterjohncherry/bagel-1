#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BraKet_Base::BraKet_Base( std::shared_ptr<Op_Info> multiop_info, 
                          std::pair<double, double> factor, int bra_num, int ket_num,  std::string type) :
                          multiop_info_(multiop_info), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
                          type_(type) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet_Base::BraKet_Base" << endl;

  cout << "BraKet_Base->type = " << type << endl;

  if (type_[0] == 'c' ) {
    name_ = "c_{I}"; 
  } else {
    name_ = ""; 
  }
  name_ += "<" + std::to_string(bra_num)+ "| ";   name_ += multiop_info->name_;   name_ += " |"+ std::to_string(ket_num) + ">";
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BraKet_Base::print_gamma_Atensor_contractions(shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo_Base> >>>> G_to_A_map,
		                                        bool has_orb_exc ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout <<  "BraKet_Base::print_gamma_Atensor_contractions()" << endl; 
 
  cout << "no proj" << endl;
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){
  
    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;
  
    if ( map_it->second->size() == 0 ) {
  
      cout << endl << "      - No Contributions! -" <<  endl << endl;
  
    } else {
  
      for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
        cout <<  A_map_it->first << "  "; cout.flush();
        for ( int qq = 0; qq !=A_map_it->second->id_orders().size(); qq++) {
          cout << "{ " ; print_vector( A_map_it->second->id_order(qq), "" ); 
          cout << " (" << A_map_it->second->factor(qq).first <<  "," <<  A_map_it->second->factor(qq).first<<  ") }" ; cout <<endl;
        }
        cout << endl;
      }
    }
  }

}
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

  cout << "hello" << endl; 
  print_vector( *(multiop_info_->op_order_) ," op_info->op_order" ); cout << endl;

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
 
              cout << "We need these blocks : " ; cout.flush(); cout << " Total_Op_->sub_tensops().size() = " ; cout.flush(); cout << Total_Op_->sub_tensops().size() << endl;
              vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();

              for (  auto& block :  *(range_map_it->second->range_blocks_canonical()) ){
                cout << block->name() << " is a required block "; cout.flush(); print_vector( *(block->orig_rngs_), "rb_it->orig_rngs_" ); cout << endl;
                MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
                required_blocks->emplace( block );
              }

              cout << endl;
            }
          }
        }
      }
    }
  } else { 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
      // TODO if is only here for non-relativistic case ; constraints should be specified in input
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
        GGen->add_gamma( range_map_it->second, trans_aops );
        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
  

              for (  auto& block :  *(range_map_it->second->range_blocks_canonical()) ){
                cout << block->name() << " is a required block "; cout.flush(); print_vector( *(block->orig_rngs_), "rb_it->orig_rngs_" ); cout << endl;
                MT_map->at( block->op_info_->op_name_ )->add_required_block( block );
                required_blocks->emplace( block );
              }


              }
              cout << endl;
            }
          }
        }
      }
    }

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

  print_gamma_Atensor_contractions( G_to_A_map, false );

  return; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet_OrbExcDeriv<DataType>::generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string, std::shared_ptr<TensOp_Base>>> MT_map,                
                                                                        std::shared_ptr<std::map<std::string,
                                                                                        std::shared_ptr<std::map<std::string,
                                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> >> block_G_to_A_map,
                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                                                        std::shared_ptr<StatesInfo_Base> target_states,
                                                                        std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks,
                                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map ) {  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet_OrbExcDeriv<DataType>::generate_gamma_Atensor_contractions : " << name_ << endl;

  Total_Op_ = MT_map->at( multiop_info_->op_name_ );

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( *(multiop_info_->op_order_),  *(multiop_info_->transformations_) );

  Total_Op_->generate_ranges( multiop_info_ );

  shared_ptr<GammaGenerator_OrbExcDeriv<DataType>> GGen = make_shared<GammaGenerator_OrbExcDeriv<DataType>>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, block_G_to_A_map, factor_, target_op_ );

  auto all_ranges = Total_Op_->all_ranges_state_specific_->at( multiop_info_->name_ );

  if ( multiop_info_->op_info_vec()->size() > 1 ) { 
 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
      // TODO if is only here for non-relativistic case ; constraints should be specified in input
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
        GGen->add_gamma( range_map_it->second, trans_aops );
        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
              cout << "We need these blocks : " ; cout.flush(); cout << " Total_Op_->sub_tensops().size() = " ; cout.flush(); cout << Total_Op_->sub_tensops().size() << endl;
              vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();
              for (  auto& rb_it :  *(range_map_it->second->range_blocks()) ){
                cout << rb_it->unique_block_->name() << " is a required block " << endl;
                MT_map->at( rb_it->unique_block_->op_info_->op_name_ )->add_required_block( rb_it->unique_block_ );
                required_blocks->emplace( rb_it->unique_block_ );
              }
              cout << endl;
            }
          }
        }
      }
    }
  } else { 
    for ( auto range_map_it = all_ranges->begin(); range_map_it != all_ranges->end(); range_map_it++ ){
      // TODO if is only here for non-relativistic case ; constraints should be specified in input
      if ( !(range_map_it->second->ci_sector_transition_ ) ) {
        GGen->add_gamma( range_map_it->second, trans_aops );
        if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
          if ( GGen->generic_reorderer( "normal order", false, false ) ) {
            if ( GGen->generic_reorderer( "alternating order", false, true ) ){
  
                cout << range_map_it->second->unique_block_->name() << " is a required block " << endl;
                MT_map->at( Total_Op_->name() )->add_required_block( range_map_it->second->unique_block_);
                required_blocks->emplace( range_map_it->second->unique_block_ );
              }
              cout << endl;
            }
          }
        }
      }
    }

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );
  
  //TODO This looks like it's going to overlap a load of stuff, also,  not sure if this is state specific
  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

  cout << " ggac 8  " << endl;
//  print_gamma_Atensor_contractions( G_to_A_map, false );

  return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet_Full<double>;
template class BraKet_Full<std::complex<double>>;
template class BraKet_OrbExcDeriv<double>;
template class BraKet_OrbExcDeriv<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
