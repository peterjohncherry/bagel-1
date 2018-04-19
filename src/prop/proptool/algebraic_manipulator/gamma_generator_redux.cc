#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaGeneratorRedux::GammaGeneratorRedux( shared_ptr<StatesInfo<double>> target_states, int Bra_num, int Ket_num,
                                          shared_ptr<TensOp_Base> total_op,
                                          shared_ptr<map<string, shared_ptr<GammaInfo>>>& Gamma_map_in,
                                          shared_ptr<map<string, shared_ptr<map<string, shared_ptr<AContribInfo> >>>>& G_to_A_map_in,
                                          double bk_factor                                                           ):
                                          target_states_(target_states),
                                          Bra_names_(target_states_->civec_names( Bra_num )),
                                          Ket_names_(target_states_->civec_names( Ket_num )),
                                          total_op_(total_op), std_ids_(total_op->idxs()), std_aops_(total_op->aops()),
                                          G_to_A_map(G_to_A_map_in), Gamma_map(Gamma_map_in), 
                                          bk_factor_(bk_factor), orig_aops_half_size_( std_aops_->size()/2 ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::GammaGeneratorRedux" << endl;

  print_vector(*std_aops_, "std_aops"); cout << endl; // This should be constant for all range blocks, but this is not the same as the MT aops
  print_vector(*std_ids_, "std_ids"); cout << endl; // ids are not (necessarily) be constant for all range blocks
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::add_gamma( const shared_ptr<Range_Block_Info> block_info, shared_ptr<vector<bool>> trans_aops ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGeneratorRedux::add_gamma " << endl;

  block_idxs_ = vector<string>( std_ids_->size() );  // TODO get this with reordering

  {
  vector<int>::iterator it_it = block_info->idxs_trans()->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ ) {
    *bi_it = (*std_ids_)[ *it_it ];
  }
  }

  block_aops_ = trans_aops;
  block_aops_rngs_ = block_info->orig_rngs_ch();

  block_rngs_ = block_info->orig_rngs();
  idxs_trans_ = block_info->idxs_trans();
  shared_ptr<vector<int>>  idxs_trans_inverse_ = block_info->idxs_trans_inverse();

  std_rngs_ = *(block_info->unique_block_); // This still needs to be transformed... 

  standard_order_ = *(block_info->idxs_trans()); // Note that we should only have block_trans, not range trans; rng_trans does not
                                                 // necessarily transform unique block into the original block (e.g.,  aabb -> bbaa requires no reordering )

  cout << endl;
  cout << "--------------- gamma def -------------------" << endl;
  print_vector( std_rngs_ ,        " unique_block_      "); cout <<endl;
  print_vector( standard_order_ ,  " range_reordering   "); cout << endl;
  print_vector(*block_rngs_ ,      " orig_rngs          "); cout <<endl;
  cout << endl;

  { // TEST for transformations
  vector<string> unique_block_dupe( std_rngs_.size());
  vector<string>::const_iterator br_it = block_rngs_->begin();
  for ( vector<int>::iterator it_it = idxs_trans_->begin() ;  it_it != idxs_trans_->end() ; it_it++, br_it++ ) { 
    unique_block_dupe[*it_it]  = *br_it;
  }

  if ( unique_block_dupe != std_rngs_ ) { 
    print_vector(unique_block_dupe, "unique_block_dupe" ) ;
    print_vector(std_rngs_, " != std_rngs_" ) ;
    throw logic_error( " reordering is broken " ) ;
  }
  }

  int ii = 0 ;
  block_to_std_order_ = vector<int>(standard_order_.size());
  for ( vector<int>::iterator so_it = standard_order_.begin() ; so_it != standard_order_.end() ; ++so_it, ++ii ) 
    block_to_std_order_[*so_it] = (ii);

  shared_ptr<vector<int>> ids_pos = make_shared<vector<int>>( std_rngs_.size() );
  iota( ids_pos->begin(), ids_pos->end(), 0 );

  shared_ptr<vector<pair<int,int>>> deltas_pos = make_shared<vector<pair<int,int>>>(0);
  int my_sign = 1; // TODO should be double from range_block

  gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>( 1, make_shared<GammaIntermediateRedux>( ids_pos, deltas_pos, my_sign ) );
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);
  cout << " leaving add_gamma" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::generic_reorderer( string reordering_name, bool first_reordering, bool final_reordering ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux::generic_reorderer" << endl; 
  
  int kk = 0;
  bool does_it_contribute = false;
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);
  
   for ( string bra_name : *Bra_names_ )
     for ( string ket_name : *Ket_names_ )
       does_it_contribute = generic_reorderer_different_sector( reordering_name, bra_name, ket_name, final_reordering );

   return does_it_contribute;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::generic_reorderer_different_sector( string reordering_name, string bra_name,
                                                              string ket_name, bool final_reordering   ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::generic_reorderer_different_sector" << endl;
//
//  if ( final_reordering ) { 
//    cout << "final_reordering !!!" << endl; 
//  } else { 
//    cout << "NOT final_reordering !!!" << endl; 
//  }
 
  shared_ptr<map<char,int>> bra_hole_map = target_states_->hole_range_map(bra_name);;
  shared_ptr<map<char,int>> bra_elec_map = target_states_->elec_range_map(bra_name);;
  shared_ptr<map<char,int>> ket_hole_map = target_states_->hole_range_map(ket_name);;
  shared_ptr<map<char,int>> ket_elec_map = target_states_->elec_range_map(ket_name);;

  if ( reordering_name == "normal order" ) {
    int kk = 0;
    while ( kk != gamma_vec->size()) {
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) 
        normal_order(kk);
      kk++;
    }
    kk = 0;
    while ( kk != gamma_vec->size()){
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) && ( gamma_vec->at(kk)->ids_pos->size() != 0 ) )
        final_gamma_vec->push_back( gamma_vec->at(kk ));
      kk++;
    }

  } else if ( reordering_name == "anti-normal order" ) {
    int kk= 0 ;
    while ( kk != gamma_vec->size()){
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) 
        anti_normal_order(kk);
      kk++;
    }
    
    kk = 0;
    while ( kk != gamma_vec->size()){
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) && ( gamma_vec->at(kk)->ids_pos->size() != 0 ) ) 
        final_gamma_vec->push_back( gamma_vec->at(kk ));
      kk++ ;
    } 
  } else if ( reordering_name == "alternating order" ) {
    int kk = 0;
    while ( kk != gamma_vec->size()) {
      alternating_order(kk);
      kk++;
    } 
    for ( shared_ptr<GammaIntermediateRedux>& gint : *gamma_vec ){
      if ( proj_onto_map( gint, *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ){ 
        final_gamma_vec->push_back( gint );
      }
    }
  }
  gamma_vec = final_gamma_vec;
  bool does_it_contribute = (gamma_vec->size() > 0 );
 

  int kk = 0;
  if ( does_it_contribute ) 
    for ( auto gint : *gamma_vec )
      print_gamma_intermediate( gint );


  if ( final_reordering && does_it_contribute ) { 
    while ( kk != gamma_vec->size()){
      add_Acontrib_to_map( kk, bra_name, ket_name );
      kk++;
    } 
  }

  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::proj_onto_map( shared_ptr<GammaIntermediateRedux> gint, 
                                         map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                         map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux::proj_onto_map" << endl;

  shared_ptr<vector<int>> idxs_pos =  gint->ids_pos;
 
  for ( vector<int>::reverse_iterator ip_it = gint->ids_pos->rbegin(); ip_it !=  gint->ids_pos->rend(); ip_it++ ) {
    char rng = (*block_aops_rngs_)[*ip_it];  
    if( ! ( (*block_aops_)[*ip_it] ) ){
      auto ket_elec_map_loc = ket_elec_map.find( rng );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        return false;
      } else if ( (ket_elec_map_loc->second -= 1 ) == -1  ) {
        return false;
      }
      auto ket_hole_map_loc = ket_hole_map.find( rng );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        ket_hole_map.emplace( rng, 1 );
      } else {
        ket_hole_map_loc->second += 1;
      }

    } else {
      auto ket_hole_map_loc = ket_hole_map.find( rng );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        return false;
      } else if ( (ket_hole_map_loc->second -= 1 ) == -1  ) {
        return false;
      }
      auto ket_elec_map_loc = ket_elec_map.find( rng );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        ket_elec_map.emplace( rng, 1 );
      } else {
        ket_elec_map_loc->second += 1;
      }
    }
  }

  for ( auto elem : ket_elec_map )  
    if ( elem.second != bra_elec_map.at(elem.first) )
       return false;     

  for ( auto elem : ket_hole_map )  
    if ( elem.second != bra_hole_map.at(elem.first) )
       return false;     

  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::anti_normal_order( int kk ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::anti_normal_order" << endl;

  shared_ptr<vector<int>> ids_pos  = gamma_vec->at(kk)->ids_pos;
  int num_kill = 0;
  for ( int pos : *ids_pos )
    if (!block_aops_->at(pos)) 
      num_kill++; 

  num_kill--; 

  for (int ii = ids_pos->size()-1 ; ii != -1; ii--){

    if ( ii > num_kill ) {
      if (block_aops_->at(ids_pos->at(ii)))
        continue;

      while(!block_aops_->at( ids_pos->at(ii) )){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if (block_aops_->at( ids_pos->at(jj) )){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }

    } else if (ii <= num_kill) {
      if (!block_aops_->at(ids_pos->at(ii)))
        continue;

      while(block_aops_->at(ids_pos->at(ii))){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if(!block_aops_->at(ids_pos->at(jj)) ){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::normal_order( int kk ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux::normal_order" << endl;

  shared_ptr<vector<int>> ids_pos  = gamma_vec->at(kk)->ids_pos;
  int num_plus = 0;
  for ( int pos : *ids_pos )
    if (block_aops_->at(pos)) 
      num_plus++; 

  num_plus--; 
 
  for (int ii = ids_pos->size()-1 ; ii != -1; ii--){

    if ( ii > num_plus ) {
      if (!block_aops_->at(ids_pos->at(ii)))
        continue;

      while(block_aops_->at( ids_pos->at(ii) )){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if (!block_aops_->at( ids_pos->at(jj) )){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }

    } else if (ii <= num_plus) {
      if (block_aops_->at(ids_pos->at(ii)))
        continue;

      while(!block_aops_->at(ids_pos->at(ii))){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if(block_aops_->at(ids_pos->at(jj)) ){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::alternating_order( int kk ) {  // e.g. +-+-+-+-
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux::alternating_order" << endl;

  shared_ptr<vector<int>> ids_pos = gamma_vec->at(kk)->ids_pos; 
  vector<int> new_ids_pos(ids_pos->size()); 
  set_standardized_alt_order_unranged( kk, new_ids_pos );

  for (int ii = ids_pos->size()-1 ; ii != -1; ii-- ){
    if (ids_pos->at(ii) == new_ids_pos[ii])
      continue;

    while( ids_pos->at(ii) != new_ids_pos[ii] ){
      for ( int jj = (ii-1); jj != -1 ; jj--) {
        if (ids_pos->at(jj) == new_ids_pos[ii]){
          swap( jj, jj+1, kk, gamma_vec);
          break;
        }
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGeneratorRedux::add_Acontrib_to_map" << endl;

  shared_ptr<GammaIntermediateRedux> gamma_int = gamma_vec->at(kk);

  double my_sign = bk_factor_*gamma_int->my_sign;
  shared_ptr<vector<pair<int,int>>> deltas_pos     = gamma_int->deltas_pos;
  shared_ptr<vector<int>> ids_pos        = gamma_int->ids_pos;

  vector<int> standardized_ids_pos( ids_pos->size() ); 
  {
  vector<int>::iterator si_it = standardized_ids_pos.begin();
  for ( vector<int>::iterator ip_it = ids_pos->begin(); ip_it != ids_pos->end(); ip_it++, si_it++ ) 
    *si_it = standard_order_[*ip_it];
  } 

  print_vector( *ids_pos,  "ids_pos" ) ; cout << endl;
  print_vector( standardized_ids_pos,  "standardized_ids_pos" ) ; cout << endl;
  print_vector( standard_order_,  "standard_order" ) ; cout << endl;

  //TODO should standardize deltas pos when building gamma intermediate so as to avoid repeated transformation
  vector<pair<int,int>> idxs_deltas_pos(deltas_pos->size());
  vector<pair<int,int>>::iterator  idp_it = idxs_deltas_pos.begin();
  for ( vector<pair<int,int>>::iterator dp_it = deltas_pos->begin(); dp_it != deltas_pos->end(); dp_it++, idp_it++ ) { 
    *idp_it = make_pair( (*idxs_trans_)[dp_it->first], (*idxs_trans_)[dp_it->second]) ;  
  } 
  print_pair_vector( *deltas_pos, "deltas_pos"); cout <<endl;
  print_pair_vector( idxs_deltas_pos, "idxs_deltas_pos"); cout <<endl;

  string Aname_alt = get_ctp_name( *std_ids_, std_rngs_, idxs_deltas_pos );
  cout << "Aname_alt = " << Aname_alt << endl;
   
  if ( total_op_->CTP_map()->find(Aname_alt) == total_op_->CTP_map()->end() ) {
    pair<double,double> ctp_factor = make_pair(1.0, 1.0); // should be my_sign/my_factor from gint
    total_op_->enter_cmtps_into_map(idxs_deltas_pos, ctp_factor, std_rngs_ );
  }

  string Gname_alt = get_gamma_name( chrvec_to_strvec(*block_aops_rngs_), *block_aops_, *ids_pos, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo>>>() );

  //TODO do this reordering w.r.t. standardized orders DQ : Is this ok? 
  vector<int> Aid_order_new = get_Aid_order( standardized_ids_pos );
  pair<double,double> new_fac = make_pair(my_sign, my_sign);
  auto AInfo_loc =  G_to_A_map->at( Gname_alt )->find(Aname_alt);
  if ( AInfo_loc == G_to_A_map->at( Gname_alt )->end() ) {
    auto AInfo = make_shared<AContribInfo_Full>( Aname_alt, Aid_order_new, new_fac );
    G_to_A_map->at( Gname_alt )->emplace(Aname_alt, AInfo) ;

  } else {
    shared_ptr<AContribInfo> AInfo = AInfo_loc->second;
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

  Gamma_map->emplace( Gname_alt, make_shared<GammaInfo>( target_states_->civec_info(bra_name), target_states_->civec_info(ket_name),
                                                         block_aops_, make_shared<const vector<string>>(chrvec_to_strvec(*block_aops_rngs_)), ids_pos, Gamma_map) );
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Swaps indexes round, flips sign, and if ranges are the same puts new density matrix in the list.
// CAREFUL : always keep creation operator as the left index in the contraction .
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::swap( int ii, int jj, int kk, shared_ptr<vector<shared_ptr<GammaIntermediateRedux>>> gamma_vec  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux::swap" << endl;

  shared_ptr<vector<int>> ids_pos              = (*gamma_vec)[kk]->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = (*gamma_vec)[kk]->deltas_pos;

  int idx_buff = (*ids_pos)[ii];  
  (*ids_pos)[ii] = (*ids_pos)[jj];
  (*ids_pos)[jj] = idx_buff;      

  if ( ( (*block_aops_rngs_)[ (*ids_pos)[jj] ] == (*block_aops_rngs_)[ (*ids_pos)[ii] ]) &&
       (*block_aops_)[ (*ids_pos)[ii] ] != (*block_aops_)[ (*ids_pos)[jj] ] ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>(*deltas_pos);
    
    //this choice (based on creation/annihilation and the reordering may be irrelevant in new scheme; reordering has to be redone at end anyhow
    pair<int,int> new_delta = (*block_aops_)[ (*ids_pos)[jj] ]  ? make_pair( (*ids_pos)[jj], (*ids_pos)[ii] ): make_pair( (*ids_pos)[ii], (*ids_pos)[jj] );
    new_deltas_tmp->push_back( new_delta );
//  shared_ptr<pint_vec> new_deltas = WickUtils::standardize_delta_ordering_generic( *new_deltas_tmp, block_idxs_ );
                                                                                           
    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>();
    for( int qq = 0 ; qq !=ids_pos->size() ; qq++)
      if ( (qq !=ii) && (qq!=jj))
        new_ids_pos->push_back( (*ids_pos)[qq] );

    shared_ptr<GammaIntermediateRedux> new_gamma = make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas_tmp, (*gamma_vec)[kk]->my_sign );
    gamma_vec->push_back(new_gamma);

  }
  (*gamma_vec)[kk]->my_sign *= -1;

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec>
GammaGeneratorRedux::standardize_delta_ordering_generic( shared_ptr<pint_vec> deltas_pos  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGeneratorRedux::standardize_delta_ordering_generic" << endl;
//TODO must order by indexes, not just by initial position
//     If one of the indexes is X, cannot "just" not contract
//     must also account for reordering ; T_{ijkl} = ... + <I| ijklmnop | J> A_{mnop} delta_{lm}
//     T_{ijkl} = .... + < I | ijklmnop | J> A_{lmop} 
//     so must flip order of A based on T (X) position

  shared_ptr<vector<pair<int,int>>> new_deltas_pos;
 
  if (deltas_pos->size() > 1 ) {
    new_deltas_pos =  make_shared <vector<pair<int,int>>>( deltas_pos->size());
    vector<int> posvec(deltas_pos->size(),0);
    for (int ii = 0 ; ii != deltas_pos->size() ; ii++)
      for (int jj = 0 ; jj != deltas_pos->size() ; jj++)
        if (block_idxs_[deltas_pos->at(ii).first] > block_idxs_[deltas_pos->at(jj).first] )
          posvec[ii]++;
 
    for (int ii = 0 ; ii != deltas_pos->size(); ii++)
      new_deltas_pos->at(posvec[ii]) = deltas_pos->at(ii);
 
  } else {
    new_deltas_pos = deltas_pos;
 
  }
  return new_deltas_pos;
}
//////////////////////////////////////////////////////////////////////////////
vector<int> GammaGeneratorRedux::get_standard_idx_order(const vector<string>&idxs) {
//////////////////////////////////////////////////////////////////////////////

  vector<int> pos(idxs.size());
  iota(pos.begin(), pos.end(), 0);

  auto idx_order_tmp = standard_order_;
  sort(pos.begin(), pos.end(), [&idxs, &idx_order_tmp](int i1, int i2){
                                   return (bool)( idxs[i1] < idxs[i2] );
                                 }
                            );
  return pos;
}
///////////////////////////////////////////////////////////////////////////////////////
// Getting reordering vec to go from Atens uncids to gamma uncids
// Running sort twice to get inverse; seems weird and there's probably a better way...
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGeneratorRedux::get_Aid_order ( const vector<int>& id_pos ) {
///////////////////////////////////////////////////////////////////////////////////////
//   cout << "GammaGeneratorRedux::get_Aid_order " << endl;

  vector<int> new_id_pos(id_pos.size());
  vector<int> tmp_order = get_position_order(id_pos);
  vector<int> new_order = get_position_order(tmp_order) ;

  return new_order;
}
///////////////////////////////////////////////////////////////////////////////////////
//Returns ordering vector for ids_pos, e.g., if ids_pos = {6,4,5} , then pos = {1,2,0}
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGeneratorRedux::get_position_order(const vector<int> &ids_pos) {
///////////////////////////////////////////////////////////////////////////////////////

  vector<int> pos(ids_pos.size());
  iota(pos.begin(), pos.end(), 0);
  sort(pos.begin(), pos.end(), [&ids_pos](int i1, int i2){return (bool)( ids_pos[i1] < ids_pos[i2] ); });

  return pos;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is preferable for ci-derivative terms; the derivative gamma matrices are by far the largest object, so 
// we want the ordering which minimizes the number of distinct gamma matrices, even if this means increasing the number of
// MO-tensor index transpositions we have to do. Consequently, we order into a canonical range sequence, by
// range, rather than by operator index. 
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGeneratorRedux::get_standardized_alt_order ( const vector<string>& rngs ,const vector<bool>& aops ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::get_standardized_alt_order " << endl;
  // TODO this should use standardized ordering 

   
  vector<int> standard_order = get_standard_idx_order(rngs) ;
  vector<int> standard_order_plus(standard_order.size()/2);
  vector<int> standard_order_kill(standard_order.size()/2);
  vector<int>::iterator standard_order_plus_it = standard_order_plus.begin();
  vector<int>::iterator standard_order_kill_it = standard_order_kill.begin();
  for ( int pos : standard_order) {
    if ( aops[pos] ) {
      *standard_order_plus_it++ = pos;
    } else {
      *standard_order_kill_it++ = pos;
    }
  }

  vector<int> standard_alt_order(standard_order.size());
  vector<int>::iterator standard_alt_order_it = standard_alt_order.begin();
  standard_order_plus_it = standard_order_plus.begin();
  standard_order_kill_it = standard_order_kill.begin();
  while ( standard_alt_order_it != standard_alt_order.end()){
    *standard_alt_order_it++ = *standard_order_plus_it++;
    *standard_alt_order_it++ = *standard_order_kill_it++;
  }

  return standard_alt_order;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This final ordering is preferable for orbital excitation derivative terms; it ensures that the indexes of the
// perturbation tensor are always in the same order, which makes combining terms easier. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::set_standardized_alt_order_unranged ( int kk , vector<int>& standard_alt_order) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::set_standardized_alt_order_unranged" << endl;
  // TODO this should use standardized ordering 

  shared_ptr<vector<int>>           ids_pos = gamma_vec->at(kk)->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 

  vector<string> unc_idxs( ids_pos->size() );
  vector<bool> unc_aops( block_aops_->size() - 2*deltas_pos->size() );

  vector<string>::iterator ui_it = unc_idxs.begin();
  vector<bool>::iterator ua_it = unc_aops.begin();
  for ( vector<int>::iterator ip_it= ids_pos->begin(); ip_it != ids_pos->end(); ip_it++, ui_it++, ua_it++ ) {
    *ui_it = block_idxs_[*ip_it];
    *ua_it = (*block_aops_)[*ip_it];
  }

  vector<int> good_order( ids_pos->size() );
  vector<int> ids_pos_standardized( ids_pos->size() );
  vector<int>::iterator ips_it = ids_pos_standardized.begin();
  for ( vector<int>::iterator ip_it  =ids_pos->begin(); ip_it != ids_pos->end(); ip_it++, ips_it++ )
    *ips_it = standard_order_[*ip_it];

  iota( good_order.begin(), good_order.end(), 0);
  sort( good_order.begin(), good_order.end(), [&ids_pos_standardized ] ( int i1, int i2) { return (bool)( ids_pos_standardized[i1] < ids_pos_standardized[i2]); });  

  vector<int> ids_reordered_pos(ids_pos->size());
  vector<int>::iterator irp_it = ids_reordered_pos.begin();
  for ( auto elem : good_order ) 
    *irp_it++ =  ids_pos->at(elem); 
  
  vector<string> ids_reordered(ids_pos->size());
  vector<string>::iterator ir_it = ids_reordered.begin();
  for ( auto elem : ids_reordered_pos ) 
    *ir_it++ =  block_idxs_[elem]; 
  
  vector<int> standard_order_plus(ids_reordered_pos.size()/2);
  vector<int> standard_order_kill(ids_reordered_pos.size()/2);
  vector<int>::iterator standard_order_plus_it = standard_order_plus.begin();
  vector<int>::iterator standard_order_kill_it = standard_order_kill.begin();
  for ( int pos : ids_reordered_pos) {
    if ( block_aops_->at(pos) ) {
      *standard_order_plus_it++ = pos;
    } else {
      *standard_order_kill_it++ = pos;
    }
  }

  vector<int>::iterator standard_alt_order_it = standard_alt_order.begin();
  standard_order_plus_it = standard_order_plus.begin();
  standard_order_kill_it = standard_order_kill.begin();
  while ( standard_alt_order_it != standard_alt_order.end()){
    *standard_alt_order_it++ = *standard_order_plus_it++;
    *standard_alt_order_it++ = *standard_order_kill_it++;
  }

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void
GammaGeneratorRedux::print_gamma_intermediate( shared_ptr<GammaIntermediateRedux> gint ) { 
/////////////////////////////////////////////////////////////////////////////////////////////
     print_vector( *(gint->ids_pos), "gint_ids_pos" ); cout  << endl;
     cout << "gint_aops = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << block_aops_->at(pos) << " " ; cout.flush();  }  cout << "] " << endl;
     cout << "gint_rngs = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << (*block_aops_rngs_)[pos] << " " ; cout.flush();  }   cout << "]" << endl;
     cout << "gint_ids  = [ "; cout.flush();   for ( auto pos : *(gint->ids_pos) ) { cout << block_idxs_[pos] << " " ; cout.flush();  }   cout << "] " <<  endl;
    return;
} 
