#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/debugging_utils.h>

using namespace std;
using namespace WickUtils;


#define __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::add_gamma( const shared_ptr<Range_Block_Info> block_info, shared_ptr<vector<bool>> trans_aops ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
cout << "void GammaGeneratorRedux<DataType>::add_gamma " << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  block_info_ = block_info;
  block_aops_ = trans_aops;
  block_aops_rngs_ = block_info->orig_rngs_ch();

  op_info_ = block_info_->op_info();
  idxs_trans_ = block_info_->idxs_trans();
  std_rngs_ = *(block_info_->unique_block_->orig_rngs_);
  standard_order_ = *(block_info_->idxs_trans());

  {
  vector<int> ids_pos_raw( std_rngs_.size() );
  iota( ids_pos_raw.begin(), ids_pos_raw.end(), 0 );
  gamma_vec_.push_back( make_unique<GammaIntermediate_Redux_Raw<DataType>>( ids_pos_raw, vector<pair<int,int>>(0), block_info_->factors() ) );
  }

#ifdef __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
  WickUtils::print_vector( gamma_vec_.front()->ids_pos_ , "gamma_vec_.front()->ids_pos_"); cout << endl;
  WickUtils::print_vector( standard_order_, "standard_order_"); cout << endl;
  WickUtils::print_vector( *std_ids_, "standard_ids_"); cout << endl;
  Debugging_Utils::print_vector_at_pos( *std_ids_, standard_order_, "std_ids standard_order_"); cout << endl;
  Debugging_Utils::print_vector_at_pos( *block_aops_rngs_, standard_order_, "block_aops_rngs , standard_order_"); cout << endl;
  block_idxs_ = vector<string>( std_ids_->size() );
  vector<int>::iterator it_it = block_info_->idxs_trans()->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ ) 
    *bi_it = (*std_ids_)[ *it_it ];
  print_gamma_intermediate(gamma_vec_.front(), "initial gamma intermediate" );  
#endif

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
cout << "void GammaGeneratorRedux<DataType>::add_Acontrib_to_map" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  unique_ptr<GammaIntermediate_Base_Raw>& gint = gamma_vec_[kk];
 
  string Gname_alt = get_gamma_name(chrvec_to_strvec(*block_aops_rngs_), *block_aops_, gint->ids_pos_, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo_Base>>>() );

  Gamma_map_->emplace( Gname_alt, make_shared<GammaInfo<DataType>>( target_states_->civec_info(bra_name), target_states_->civec_info(ket_name),
                                                                    *block_aops_, chrvec_to_strvec(*block_aops_rngs_), gint->ids_pos_, Gamma_map_ ) );

  vector<pair<int,int>> idxs_deltas_pos = gint->deltas_pos_; 
  transform_to_canonical_ids_pos( idxs_deltas_pos ); 

  string Aname_alt = get_ctp_name( op_info_->op_state_name_canonical(), *std_ids_,  std_rngs_, idxs_deltas_pos );

  if ( total_op_->CTP_map()->find(Aname_alt) == total_op_->CTP_map()->end() ) 
    total_op_->enter_cmtps_into_map(idxs_deltas_pos, std_rngs_, op_info_ );

  vector<int> standardized_ids_pos = gint->ids_pos_;  
  transform_to_canonical_ids_pos( standardized_ids_pos ); 
 
  vector<int> Aid_order_new = get_Aid_order( standardized_ids_pos );
  pair<double,double> new_fac = gint->factors_; 
  auto AInfo_loc =  G_to_A_map->at( Gname_alt )->find(Aname_alt);
  if ( AInfo_loc == G_to_A_map->at( Gname_alt )->end() ) {
    auto AInfo = make_shared<AContribInfo_Full<DataType>>( Aname_alt, Aid_order_new, new_fac );
    G_to_A_map->at( Gname_alt )->emplace(Aname_alt, AInfo) ;

  } else {
    shared_ptr<AContribInfo_Base> AInfo = AInfo_loc->second;
    for ( int qq = 0 ; qq != AInfo->id_orders().size(); qq++ ) {
      if( Aid_order_new == AInfo->id_order(qq) ){
        AInfo->combine_factors( qq, new_fac );
        AInfo->remaining_uses_ += 1;
        AInfo->total_uses_ += 1;
        break;

      } else if ( qq == AInfo->id_orders().size()-1) {
        AInfo->add_id_order(Aid_order_new);
        AInfo->add_factor(new_fac);
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::block_trans_test( unique_ptr<GammaIntermediate_Base_Raw>& gint ){
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_GENERATOR_REDUX
cout << "GammaGeneratorRedux<DataType>::block_trans_test " << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////

  block_idxs_ = vector<string>( std_ids_->size() );
  vector<int>::iterator it_it = block_info_->idxs_trans()->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ ) 
    *bi_it = (*std_ids_)[ *it_it ];

  print_vector( block_idxs_, " block_idxs " ); cout << endl;
  print_vector( *std_ids_,  " std_ids_   " ); cout << endl;

  print_pair_vector( gint->deltas_pos_ , "gint->deltas_pos_" ); cout.flush();
  cout << "  : [ ";
  for ( auto ctr : gint->deltas_pos_ ) {
    cout << "(" << block_idxs_[ctr.first] <<","<<  block_idxs_[ctr.second] << ")"; cout.flush();
  }
  cout << "]" << endl;

  vector<pair<int,int>> idxs_deltas_pos = gint->deltas_pos_;
  transform_to_canonical_ids_pos( idxs_deltas_pos );

  print_pair_vector( idxs_deltas_pos , "idxs_deltas_pos_" );
  cout << "  : [ ";
  for ( auto ctr : idxs_deltas_pos ) {
    cout << "(" << (*std_ids_)[ctr.first] <<","<<  (*std_ids_)[ctr.second] << ")"; cout.flush();
  }
  cout << "]" << endl;
 
  cout << "  : [ ";
  for ( auto ctr : idxs_deltas_pos ) {
    cout << "(" << (std_rngs_)[ctr.first] <<","<<  (std_rngs_)[ctr.second] << ")"; cout.flush();
  }
  cout << "]" << endl;

  return;
}
//////////////////////////////////////////////////////////////////////////
template class GammaGeneratorRedux<double>;
template class GammaGeneratorRedux<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
