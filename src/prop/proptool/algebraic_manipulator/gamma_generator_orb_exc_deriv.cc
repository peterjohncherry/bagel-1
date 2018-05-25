#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>

using namespace std;
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::add_gamma( const shared_ptr<Range_Block_Info> block_info, shared_ptr<vector<bool>> trans_aops  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
cout << "void GammaGenerator_OrbExcDeriv::add_gamma " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////

  idxs_trans_ = block_info->idxs_trans();
  block_aops_ = trans_aops;
  block_aops_rngs_ = block_info->orig_rngs_ch(); // May be transformed by Bra or Ket.

  op_info_ = block_info->op_info();

  block_rngs_ = make_shared<vector<string>>( block_info->orig_rngs()->size() ) ; // After the below definition, this remains unchanged.
  vector<string>::iterator br_it = block_rngs_->begin();
  for ( vector<int>::iterator it_it = idxs_trans_->begin() ; it_it != idxs_trans_->end(); it_it++ , br_it++ )
    *br_it = (*(block_info->unique_block_->orig_rngs_))[*it_it];

  std_rngs_ = *(block_info->unique_block_->orig_rngs_);
  standard_order_  = *(block_info->idxs_trans());

  // TODO change; this will be problematic when we have multiple operators with the same name but different states.
  target_block_start_ = 0 ;
  target_block_end_ = 0 ;
  for ( vector<int>::iterator oo_it = block_info->op_info_->op_order_->begin() ; oo_it != block_info->op_info_->op_order_->end(); oo_it++ ) {
    shared_ptr<Range_Block_Info>& rb = (*(block_info->range_blocks()))[*oo_it];
    target_block_end_ += rb->num_idxs_; 
    if ( rb->op_info_->op_state_name_[0] == target_op_[0] ) {
      target_block_name_ = rb->name(); 
      target_block_size_ = rb->num_idxs_; 
      break; 
    }
    target_block_start_ += rb->num_idxs_; 
  }

  // TMP.
  std_rngs_target_op_free_ = vector<string>(block_aops_->size() - target_block_size_ ); 
  std_idxs_target_op_free_ = vector<string>(block_aops_->size() - target_block_size_ ); 
  { 
  std_name_target_op_free_ = "";
  shared_ptr<vector<shared_ptr<Range_Block_Info>>> range_blocks = block_info->range_blocks() ; 
  vector<string>::const_iterator si_it = std_ids_->begin();
  vector<string>::iterator sr_it = std_rngs_.begin();
  vector<string>::iterator sitof_it = std_idxs_target_op_free_.begin();
  vector<string>::iterator srtof_it = std_rngs_target_op_free_.begin();

  for( vector<int>::iterator oo_it = op_info_->op_order_->begin(); oo_it != op_info_->op_order_->end(); oo_it++ ) {
    if ( (*range_blocks)[*oo_it]->full_op_name()[0] != target_op_[0] ){ 
      for ( int rr = 0  ; rr != (*range_blocks)[*oo_it]->num_idxs_; sr_it++, si_it++, srtof_it++, sitof_it++, rr++ ) {
        *sitof_it = *si_it;
        *srtof_it = *sr_it;
      }  
      std_name_target_op_free_ +=  (*range_blocks)[*oo_it]->op_state_name_;
    } else { 
      si_it += target_block_size_; 
      sr_it += target_block_size_;
    }
  }
  }

  shared_ptr<vector<int>> ids_pos = make_shared<vector<int>>( std_rngs_.size() );
  iota( ids_pos->begin(), ids_pos->end(), 0 );

  shared_ptr<vector<pair<int,int>>> A_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> target_target_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> target_A_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  
  //TODO  Change to specialized class
  gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediate_Base >>>(1);
  gamma_vec_->front() = make_shared<GammaIntermediate_OrbExcDeriv<DataType>>( ids_pos, A_deltas_pos, target_A_deltas_pos, target_target_deltas_pos, block_info->factors());

  final_gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediate_Base >>>(0);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::swap( int ii, int jj, int kk ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
cout << "GammaGenerator_OrbExcDeriv<DataType>::swap ii = " << ii << " jj = " << jj << " kk = " << kk << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<GammaIntermediate_Base>& gint =  gamma_vec_->at(kk);

  // only need one buffer really, but this looks clearer
  int j_pos = (*((gint)->ids_pos_))[jj];
  int i_pos = (*((gint)->ids_pos_))[ii];

  (*((gint)->ids_pos_))[ii] = j_pos;
  (*((gint)->ids_pos_))[jj] = i_pos;

  if ( ( (*block_aops_rngs_)[ j_pos ] == (*block_aops_rngs_)[ i_pos ]) &&( (*block_aops_)[ i_pos ] != (*block_aops_)[ j_pos ]) ){

      i_pos = standard_order_[i_pos ]; 
      j_pos = standard_order_[j_pos ]; 

      shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( gint->ids_pos_->size()-2);
      vector<int>::iterator nip_it = new_ids_pos->begin();
      vector<int>::iterator gip_it = gint->ids_pos_->begin();
      for( int qq = 0 ; qq != gint->ids_pos_->size() ; qq++, gip_it++) {
        if( (qq != ii) && (qq != jj) ){
         *nip_it = *gip_it;
          ++nip_it;
        }
      }

      shared_ptr<pint_vec> new_target_target_deltas_pos = make_shared<pint_vec>( *(gint->target_target_deltas_pos()) );
      shared_ptr<pint_vec> new_target_A_deltas_pos = make_shared<pint_vec>( *(gint->target_A_deltas_pos()) );
      shared_ptr<pint_vec> new_A_A_deltas_pos = make_shared<pint_vec>( *(gint->deltas_pos_) );

      //Make it so T is first index, if all or no T, make it so creation op is first index
      if ( (j_pos >= target_block_start_) && ( j_pos < target_block_end_ )  ) {

        if (  (i_pos >= target_block_start_) &&  ( i_pos < target_block_end_ ) ) {
          pair<int,int> new_delta = (*block_aops_)[ j_pos ] ? make_pair( j_pos, i_pos ) : make_pair( i_pos , j_pos);
	  new_target_target_deltas_pos->push_back(new_delta);
	  
        } else {
          pair<int,int> new_delta = make_pair( j_pos, i_pos );
	  new_target_A_deltas_pos->push_back(new_delta);
        }

      } else if ( i_pos >= target_block_start_ && i_pos < target_block_end_ ) {
        pair<int,int> new_delta = make_pair( i_pos, j_pos );
        new_target_A_deltas_pos->push_back(new_delta);

      } else {  
        pair<int,int> new_delta = (*block_aops_)[ j_pos ] ? make_pair( j_pos, i_pos ) : make_pair( i_pos , j_pos);
        new_A_A_deltas_pos->push_back(new_delta);
      } 
  
      shared_ptr<GammaIntermediate_OrbExcDeriv<DataType>> new_gamma = 
           make_shared<GammaIntermediate_OrbExcDeriv<DataType>>( new_ids_pos, new_A_A_deltas_pos, new_target_A_deltas_pos, new_target_target_deltas_pos, (gint)->factors_ );
      
      gamma_vec_->push_back(new_gamma);

  }
//  gint->factors_.first*=-1.0;
//  gint->factors_.second*=-1.0;

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
cout << "GammaGenerator_OrbExcDeriv<DataType>::add_Acontrib_to_map" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<GammaIntermediate_Base> gint = (*gamma_vec_)[kk];
  vector<int> gamma_ids_pos = *(gint->ids_pos_);

  string Gname_alt = get_gamma_name( chrvec_to_strvec(*block_aops_rngs_), *block_aops_, gamma_ids_pos, bra_name, ket_name );
  if ( Gamma_map_->find(Gname_alt) == Gamma_map_->end() ) 
    Gamma_map_->emplace( Gname_alt, make_shared<GammaInfo<DataType>>( target_states_->civec_info(bra_name), target_states_->civec_info(ket_name),
                                                                     block_aops_, make_shared<const vector<string>>(chrvec_to_strvec(*block_aops_rngs_)),
                                                                     gint->ids_pos_, Gamma_map_) );

  //This is the name of the CTP we need to contract (perhaps partially) with gamma
  string Aname_alt = get_ctp_name( std_name_target_op_free_, std_idxs_target_op_free_, std_rngs_target_op_free_, *(gint->deltas_pos_) );
  if ( total_op_->CTP_map()->find(Aname_alt) == total_op_->CTP_map()->end() )
    total_op_->enter_cmtps_into_map( *(gint->deltas_pos_), std_rngs_, op_info_ );

  transform_to_canonical_ids_pos( gamma_ids_pos );
 
  vector<int> gamma_contraction_pos(0);
  vector<int> A_contraction_pos(0);
  shared_ptr<vector<int>> T_pos = make_shared<vector<int>>(0);
  shared_ptr<vector<string>> post_gamma_contraction_rngs = make_shared<vector<string>>(0);
  // Get A indexes which must be contracted with gamma indexes
  // Get T indexes which are not contracted with A, and are basically gamma indexes
  int pos_in_gamma = 0;
  for( vector<int>::iterator gip_it = gamma_ids_pos.begin(); gip_it != gamma_ids_pos.end() ; gip_it++, pos_in_gamma++ ){
    if ( ( *gip_it < target_block_start_) ||  ( *gip_it >= target_block_end_ )  ){ 

      if ( *gip_it < target_block_start_)  { 
        A_contraction_pos.push_back(*gip_it);
      } else { 
        A_contraction_pos.push_back(*gip_it - target_block_size_ );
      }
      gamma_contraction_pos.push_back(pos_in_gamma);

    } else {
      post_gamma_contraction_rngs->push_back(std_rngs_[*gip_it]);
      T_pos->push_back(*gip_it - target_block_start_ );
    }
  }
 
  // Get A indexes which are not contracted with gamma; these are T indexes.
  vector<int> A_T_pos(0);
  if ( gint->target_A_deltas_pos()->size() != 0 ) {
    for ( pair<int,int>& ctr : *(gint->target_A_deltas_pos()) ) { 

      if ( ctr.second < target_block_start_ ) {
        A_T_pos.push_back(ctr.second);
      } else {
        A_T_pos.push_back(ctr.second + target_block_end_);
      } 
      post_gamma_contraction_rngs->push_back(std_rngs_[ctr.first]);
      T_pos->push_back(ctr.first - target_block_start_ );

    }
  }

  // Might have
  // Tw Tx Ty Tz  = < | ax ay ai al | J > Ai A(j/w) A(k/z) Al   
  //
  // Sort A to contract i and l with gamma
  // Tw Tx Ty Tz  = < | ax ay | J > A'(j/w) A'(k/z)       : [ ijkl -> iljk ]
  //
  // Sort remaining indexes so they match up with T order
  // Tw Tx Ty Tz  = < | ax ay | J > A'(j/w) A'(k/z)       : [ xyjk ->jxyk ] or [ xywz -> wxyz ]
  //
  
  // The latter case corresponds to this reordering (not we get x and y in the first loop, and j/w and k/z in the second
  shared_ptr<vector<int>> post_contraction_reordering = get_ascending_order (*T_pos );

  shared_ptr<vector<int>> pre_contraction_reordering;
  vector<int> A_ids_pos( A_contraction_pos.size() + A_T_pos.size() );
  copy ( A_contraction_pos.begin(), A_contraction_pos.end(), A_ids_pos.begin() );
  copy ( A_T_pos.begin(), A_T_pos.end(), A_ids_pos.begin() + A_contraction_pos.size() );
 
  pre_contraction_reordering  =  make_shared<vector<int>> ( get_Aid_order( A_ids_pos ) ) ;

  shared_ptr<map<string, shared_ptr<map<string, shared_ptr<AContribInfo_Base>>>>> G_to_A_map;
  auto map_loc = block_G_to_A_map_->find(target_block_name_);
  if ( map_loc == block_G_to_A_map_->end() ) {
    G_to_A_map = make_shared<map<string, shared_ptr<map<string, shared_ptr<AContribInfo_Base>>>>>() ;
    block_G_to_A_map_->emplace( target_block_name_,  G_to_A_map );
  } else {  
    G_to_A_map = map_loc->second;
  } 

#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
  print_target_block_info( gamma_ids_pos, A_ids_pos, *T_pos, A_T_pos, A_contraction_pos, gamma_contraction_pos,       
                           *pre_contraction_reordering, *post_contraction_reordering, *post_gamma_contraction_rngs);
#endif

  //Should be bk_factor, but something is going wrong...
  pair<double,double> new_fac = make_pair(1.0, 0.0);
  pair_fac_mult( gint->factors_, new_fac );
   
  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt,  make_shared<map<string, shared_ptr<AContribInfo_Base>>>() );

  string final_reordering_name = get_final_reordering_name( Gname_alt, *post_contraction_reordering );

  shared_ptr<vector<string>> pre_contraction_ranges = make_shared<vector<string>>( A_ids_pos.size());  
  { 
    vector<string>::iterator pcr_it = pre_contraction_ranges->begin();
    for ( vector<int>::iterator aip_it =  A_ids_pos.begin(); aip_it != A_ids_pos.end(); aip_it++, pcr_it++) 
      *pcr_it = std_rngs_[ *aip_it ];
  }

  auto a_info_loc =  G_to_A_map->at( Gname_alt )->find(final_reordering_name);
  if ( a_info_loc == G_to_A_map->at( Gname_alt )->end() ) {

     auto a_info = make_shared<AContribInfo_OrbExcDeriv<DataType>>( final_reordering_name, target_block_name_, T_pos,  post_gamma_contraction_rngs );
     a_info->add_reordering( Aname_alt, gamma_contraction_pos, *pre_contraction_reordering, pre_contraction_ranges, new_fac );
     G_to_A_map->at( Gname_alt )->emplace( final_reordering_name, a_info );

  } else {
//    shared_ptr<AContribInfo_OrbExcDeriv<DataType>> AInfo = std::dynamic_pointer_cast<AContribInfo_OrbExcDeriv<DataType>>( AInfo_loc->second );
//    AInfo->add_reordering( *post_contraction_reordering, *pre_contraction_reordering, new_fac ); 
  }
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void
GammaGenerator_OrbExcDeriv<DataType>::print_target_block_info( const vector<int>& gamma_ids_pos, const vector<int>&  A_ids_pos,                  
                                                               const vector<int>& T_pos, const vector<int>& A_T_pos,                    
                                                               const vector<int>& A_contraction_pos,
                                                               const vector<int>& gamma_contraction_pos,      
                                                               const vector<int>& pre_contraction_reordering, 
                                                               const vector<int>& post_contraction_reordering,
                                                               const vector<string>& post_gamma_contraction_rngs) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
cout << "GammaGenerator_OrbExcDeriv::print_target_block_info" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "------------" << target_block_name_ << "----------"<< endl;
  print_vector( gamma_ids_pos,                "Gamma_ids_pos              " ); cout << endl;
  print_vector( A_ids_pos,                    "A_ids_pos                  " ); cout << endl;
  print_vector( T_pos,                       "T_pos                      " ); cout << endl;
  print_vector( A_T_pos,                      "A_T_pos                    " ); cout << endl;
  print_vector( A_contraction_pos,            "A_contraction_pos          " ); cout << endl;
  print_vector( gamma_contraction_pos,        "gamma_contraction_pos      " ); cout << endl;
  print_vector( pre_contraction_reordering,  "pre_contraction_reordering " ); cout << endl;
  print_vector( post_contraction_reordering, "post_contraction_reordering" ); cout << endl;
  print_vector( post_gamma_contraction_rngs, "post_gamma_contraction_rngs" ); cout << endl;
  cout << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
string
GammaGenerator_OrbExcDeriv<DataType>::get_final_reordering_name ( string Gname_alt,
                                                                  const vector<int>& post_contraction_reordering ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMA_GENERATOR_ORB_EXC_DERIV
cout << "GammaGenerator_OrbExcDeriv::get_final_reordering_name" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////

  string final_reordering_name = Gname_alt; 
  final_reordering_name += " pcr[";
  for ( auto elem  : post_contraction_reordering ){
    final_reordering_name += ( '0' + elem );
    final_reordering_name += ' ';
  }
  final_reordering_name += ']';

  return final_reordering_name;

} 

////////////////////////////////////////////////////////////////////////////
template class GammaGenerator_OrbExcDeriv<double>;
template class GammaGenerator_OrbExcDeriv<std::complex<double>>;
////////////////////////////////////////////////////////////////////////////
