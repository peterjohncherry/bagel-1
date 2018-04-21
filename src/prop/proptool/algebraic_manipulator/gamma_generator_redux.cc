#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
GammaGeneratorRedux<DataType>::GammaGeneratorRedux( shared_ptr<StatesInfo<DataType>> target_states, int Bra_num, int Ket_num,
                                                    shared_ptr<TensOp_Base> total_op,
                                                    shared_ptr<map<string, shared_ptr<GammaInfo<DataType>>>>& Gamma_map_in,
                                                    shared_ptr<map<string, shared_ptr<map<string, shared_ptr<AContribInfo> >>>>& G_to_A_map_in,
                                                    double bk_factor                                                           ):
                                                    target_states_(target_states),
                                                    Bra_names_(target_states_->civec_names( Bra_num )),
                                                    Ket_names_(target_states_->civec_names( Ket_num )),
                                                    total_op_(total_op), std_ids_(total_op->idxs()), std_aops_(total_op->aops()),
                                                    G_to_A_map(G_to_A_map_in), Gamma_map(Gamma_map_in), 
                                                    bk_factor_(bk_factor), orig_aops_half_size_( std_aops_->size()/2 ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::GammaGeneratorRedux" << endl;

  print_vector(*std_aops_, "std_aops"); cout << endl; // This should be constant for all range blocks, but this is not the same as the MT aops
  print_vector(*std_ids_, "std_ids"); cout << endl; // ids are not (necessarily) be constant for all range blocks
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::add_gamma( const shared_ptr<Range_Block_Info> block_info, shared_ptr<vector<bool>> trans_aops ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGeneratorRedux<DataType>::add_gamma " << endl;

  block_idxs_ = vector<string>( std_ids_->size() );
  {
  vector<int>::iterator it_it = block_info->idxs_trans()->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ ) 
    *bi_it = (*std_ids_)[ *it_it ];
  }

  block_aops_ = trans_aops;
  block_aops_rngs_ = block_info->orig_rngs_ch();

  block_rngs_ = block_info->orig_rngs();
  idxs_trans_ = block_info->idxs_trans();
  shared_ptr<vector<int>>  idxs_trans_inverse_ = block_info->idxs_trans_inverse();

  std_rngs_ = *(block_info->unique_block_);
  standard_order_ = *(block_info->idxs_trans());

  cout << endl;
  cout << "--------------- gamma def -------------------" << endl;
  print_vector( std_rngs_ ,        " unique_block_      "); cout <<endl;
  print_vector( standard_order_ ,  " range_reordering   "); cout << endl;
  print_vector(*block_rngs_ ,      " orig_rngs          "); cout <<endl;
  cout << endl;

  int ii = 0 ;
  block_to_std_order_ = vector<int>(standard_order_.size());
  for ( vector<int>::iterator so_it = standard_order_.begin() ; so_it != standard_order_.end() ; ++so_it, ++ii ) 
    block_to_std_order_[*so_it] = (ii);

  shared_ptr<vector<int>> ids_pos = make_shared<vector<int>>( std_rngs_.size() );
  iota( ids_pos->begin(), ids_pos->end(), 0 );

  shared_ptr<vector<pair<int,int>>> deltas_pos = make_shared<vector<pair<int,int>>>(0);
  
  pair<double,double>  factors = block_info->factors();
  cout << "factors = (" << factors.first << "," << factors.second << ")" << endl;

  gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>( 1, make_shared<GammaIntermediateRedux>( ids_pos, deltas_pos, factors ) );
  final_gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
bool GammaGeneratorRedux<DataType>::generic_reorderer( string reordering_name, bool first_reordering, bool final_reordering ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux<DataType>::generic_reorderer" << endl; 
  
  int kk = 0;
  bool does_it_contribute = false;
  final_gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);
  for ( string bra_name : *Bra_names_ ) {
    for ( string ket_name : *Ket_names_ ) {
      bra_name_ = bra_name;
      ket_name_ = ket_name;
      bra_hole_map_ = target_states_->hole_range_map(bra_name);
      bra_elec_map_ = target_states_->elec_range_map(bra_name);
      ket_hole_map_ = target_states_->hole_range_map(ket_name);
      ket_elec_map_ = target_states_->elec_range_map(ket_name);
      does_it_contribute = generic_reorderer_different_sector( reordering_name, final_reordering );
    }
  }
  return does_it_contribute;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
bool GammaGeneratorRedux<DataType>::generic_reorderer_different_sector( string reordering_name, bool final_reordering ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::generic_reorderer_different_sector" << endl;

  if ( reordering_name == "normal order" ) {
    normal_order();
    for ( vector<shared_ptr<GammaIntermediateRedux>>::iterator gv_it = gamma_vec_->begin(); gv_it != gamma_vec_->end(); gv_it++ ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) ) 
        final_gamma_vec_->push_back( *gv_it );
    }

  } else if ( reordering_name == "anti-normal order" ) {
    anti_normal_order();

    vector<shared_ptr<GammaIntermediateRedux>>::iterator gv_it = gamma_vec_->begin();
    while ( gv_it != gamma_vec_->end() ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) )
        final_gamma_vec_->push_back( *gv_it );
      gv_it++;
    }

  } else if ( reordering_name == "alternating order" ) {
    alternating_order();
    for ( vector<shared_ptr<GammaIntermediateRedux>>::iterator gv_it = gamma_vec_->begin(); gv_it != gamma_vec_->end(); gv_it++ ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) )
        final_gamma_vec_->push_back( *gv_it );
    }
  }

  gamma_vec_ = final_gamma_vec_;
  bool does_it_contribute = ( gamma_vec_->size() > 0 );

  int kk = 0;
  if ( does_it_contribute ) 
    for ( auto& gint : *gamma_vec_ )
      print_gamma_intermediate( gint );


  if ( final_reordering && does_it_contribute ) { 
    while ( kk != gamma_vec_->size()){
      add_Acontrib_to_map( kk, bra_name_, ket_name_ );
      kk++;
    } 
  }

  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
bool GammaGeneratorRedux<DataType>::proj_onto_map( shared_ptr<GammaIntermediateRedux> gint, 
                                         map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                         map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::proj_onto_map" << endl;

  shared_ptr<vector<int>> idxs_pos =  gint->ids_pos_;

  for ( vector<int>::reverse_iterator ip_it = idxs_pos->rbegin(); ip_it !=  idxs_pos->rend(); ip_it++ ) {
    char rng = (*block_aops_rngs_)[*ip_it];  

    if( !((*block_aops_)[*ip_it]) ){
      auto ket_elec_map_loc = ket_elec_map.find( rng );
      if ( ket_elec_map_loc == ket_elec_map.end() ) { cout << "failed  proj_onto_map 1" << endl;
        return false;
      } else if ( (ket_elec_map_loc->second -= 1 ) == -1  ) {cout << "failed proj_onto_map 2" << endl;
        return false;
      }
      auto ket_hole_map_loc = ket_hole_map.find( rng );

      if ( ket_hole_map_loc == ket_hole_map.end() ) {

        ket_hole_map.emplace( rng, 1 );
      } else {
        ket_hole_map_loc->second += 1;
      }
      cout << "added [" << rng << "] hole " << endl;

    } else {
      auto ket_hole_map_loc = ket_hole_map.find( rng );
      if ( ket_hole_map_loc == ket_hole_map.end() ) { cout << "failed_proj_onto_map 3" << endl;
        return false;
      } else if ( (ket_hole_map_loc->second -= 1 ) == -1  ) { cout << "failed proj_onto_map 4" << endl;
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

  // TODO  may have problems if bra and ket involve different ranges, but not relevant for now.
  for ( auto& elem : *ket_elec_map_ )  
    if ( elem.second != bra_elec_map.at(elem.first) )
      return false;     
    

  for ( auto& elem : *ket_hole_map_ )                     
    if ( elem.second != bra_hole_map.at(elem.first) )
      return false;     
    
  cout << " passed proj_onto_map " << endl;
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::normal_order" << endl;
 
  int kk = 0;
  while ( kk != gamma_vec_->size() ) {
 
    if ( proj_onto_map( gamma_vec_->at(kk), *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) ){ 

      shared_ptr<vector<int>> ids_pos = gamma_vec_->at(kk)->ids_pos_;
      int num_plus = 0;
      for ( int pos : *ids_pos )
        if ( (*block_aops_)[ pos ]) 
          num_plus++; 
      
      num_plus--; 
     
      for (int ii = ids_pos->size()-1 ; ii != -1; ii--){
      
        if ( ii > num_plus ) {
          if ( !(*block_aops_)[ (*ids_pos)[ii] ] )
            continue;
      
          while( (*block_aops_)[(*ids_pos)[ii]] ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if ( !(*block_aops_)[(*ids_pos)[jj]] ){
                swap( jj, jj+1, kk);
                break;
              }
            }
          }
      
        } else if (ii <= num_plus) {
          if ( (*block_aops_)[(*ids_pos)[ii]])
            continue;
      
          while(!(*block_aops_)[(*ids_pos)[ii]]){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if((*block_aops_)[(*ids_pos)[jj]] ){
                swap( jj, jj+1, kk);
                break;
              }
            }
          }
        }
      }
    }
    kk++; 
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::anti_normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::anti_normal_order" << endl;
 
  int kk = 0;
  cout <<" gamma_vec_->size() = " << gamma_vec_->size() << " kk = " << kk << endl; 
  while ( kk != gamma_vec_->size() ) {
    if ( proj_onto_map( gamma_vec_->at(kk), *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) ){ 
     
      shared_ptr<vector<int>> ids_pos  = gamma_vec_->at(kk)->ids_pos_;
      int num_kill = 0;
      for ( int pos : *ids_pos )
        if (!(*block_aops_)[pos]) 
          num_kill++; 
      num_kill--; 
      
      for (int ii = ids_pos->size()-1 ; ii != -1; ii--){
        if ( ii > num_kill ) {
          if ( (*block_aops_)[ (*ids_pos)[ii]] ) 
            continue;
      
          while( !(*block_aops_)[(*ids_pos)[ii]] ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if ( (*block_aops_)[(*ids_pos)[jj]] ){
                cout << "2kk  = " << kk<< endl;
                swap( jj, jj+1, kk );
                break;
              }
            }
          }
      
      cout << "3kk = " << kk << endl;
        } else if (ii <= num_kill) {
          if (!(*block_aops_)[(*ids_pos)[ii]] )
            continue;
      
          while((*block_aops_)[(*ids_pos)[ii]] ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if(!(*block_aops_)[(*ids_pos)[jj]] ) {
                cout << "4kk = " << kk << endl;
                print_vector( *(gamma_vec_->at(kk)->ids_pos_) , " gv_it->ids_pos"); cout << endl;
                swap( jj, jj+1, kk );
                break;
              }
            }
          }
        }
      }
    }
    cout << "end kk = " << kk++ << endl;
  }
  cout << "leaving anti-normal order " << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::alternating_order() {  // e.g. +-+-+-+-
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGeneratorRedux<DataType>::alternating_order" << endl;
 
//  for ( vector<shared_ptr<GammaIntermediateRedux>>::iterator gv_it = gamma_vec_->begin(); gv_it != gamma_vec_->end(); gv_it++ ) {
  int kk = 0;
  while ( kk != gamma_vec_->size() ) {
    shared_ptr<GammaIntermediateRedux> gint = gamma_vec_->at(kk);
    if ( proj_onto_map( gint, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) ){ 
      shared_ptr<vector<int>> ids_pos = gint->ids_pos_; 

      vector<int> new_ids_pos( ids_pos->size() ); 
      set_standardized_alt_order_unranged( gint, new_ids_pos );
    
      for (int ii = ids_pos->size()-1 ; ii != -1; ii-- ){
        if ( (*ids_pos)[ii] == new_ids_pos[ii])
          continue;
    
        while( (*ids_pos)[ii] != new_ids_pos[ii] ){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
            if ( (*ids_pos)[jj] == new_ids_pos[ii] ){
              swap( jj, jj+1, kk );
              break;
            }
          }
        }
      }
    }
    kk++;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGeneratorRedux<DataType>::add_Acontrib_to_map" << endl;

  shared_ptr<GammaIntermediateRedux> gint = gamma_vec_->at(kk);

  shared_ptr<vector<pair<int,int>>> deltas_pos     = gint->deltas_pos_;
  shared_ptr<vector<int>> ids_pos        = gint->ids_pos_;

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
    pair<double,double> ctp_factor = gint->factors_; // should be my_sign/my_factor from gint
    total_op_->enter_cmtps_into_map(idxs_deltas_pos, ctp_factor, std_rngs_ );
  }

  string Gname_alt = get_gamma_name( chrvec_to_strvec(*block_aops_rngs_), *block_aops_, *ids_pos, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo>>>() );

  //TODO do this reordering w.r.t. standardized orders DQ : Is this ok? 
  vector<int> Aid_order_new = get_Aid_order( standardized_ids_pos );
  pair<double,double> new_fac = make_pair( bk_factor_*gint->factors_.first , bk_factor_*gint->factors_.second );  

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

  Gamma_map->emplace( Gname_alt, make_shared<GammaInfo<DataType>>( target_states_->civec_info(bra_name), target_states_->civec_info(ket_name),
                                                         block_aops_, make_shared<const vector<string>>(chrvec_to_strvec(*block_aops_rngs_)), ids_pos, Gamma_map) );
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::swap( int ii, int jj, shared_ptr<GammaIntermediateRedux>& gint  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::swap ii = " << ii << " jj = " << jj << endl;

  int idx_buff = (*(gint->ids_pos_))[ii];  
 (*(gint->ids_pos_))[ii] = (*(gint->ids_pos_))[jj];
 (*(gint->ids_pos_))[jj] = idx_buff;      

  if ( ( (*block_aops_rngs_)[ (*(gint->ids_pos_))[jj] ] == (*block_aops_rngs_)[ (*(gint->ids_pos_))[ii] ]) &&
       (*block_aops_)[ (*(gint->ids_pos_))[ii] ] != (*block_aops_)[(*(gint->ids_pos_))[jj] ] ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>( *(gint->deltas_pos_) );
    pair<int,int> new_delta = (*block_aops_)[ ( *(gint->ids_pos_))[jj] ]  ?
                               make_pair( (*(gint->ids_pos_))[jj], ( *(gint->ids_pos_))[ii] ): make_pair((*(gint->ids_pos_))[ii],(*(gint->ids_pos_))[jj]);

    new_deltas_tmp->push_back( new_delta );

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( gint->ids_pos_->size()-2);
    vector<int>::iterator nip_it = new_ids_pos->begin();
    for( int qq = 0 ; qq != gint->ids_pos_->size() ; qq++) {
      if ( (qq != ii) && (qq != jj)){
       *nip_it = (*(gint->ids_pos_))[qq];
        ++nip_it;
      }
    }

    shared_ptr<GammaIntermediateRedux> new_gamma = make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas_tmp, gint->factors_ );

    print_vector( *(gint->ids_pos_) , " after_new_gamma gint->ids_pos "); cout << endl;
    gamma_vec_->push_back(new_gamma);
    cout << " gint->factors_ = ("<<   gint->factors_.first << "," <<  gint->factors_.second << ")" << endl;
    print_pair_vector( *(gint->deltas_pos_) , " after push_back gint->ids_pos "); cout << endl;
    print_vector( *(gint->ids_pos_) , " after push_back gint->ids_pos "); cout << endl;

  }
  pair<double,double> anti_herm_fac = make_pair( -1.0 , 1.0 );
  WickUtils::pair_fac_mult( anti_herm_fac, gint->factors_  );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::swap( int ii, int jj, vector<shared_ptr<GammaIntermediateRedux>>::iterator gv_it  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::swap ii = " << ii << " jj = " << jj << endl;

  int idx_buff = (*((*gv_it)->ids_pos_))[ii];  
 (*((*gv_it)->ids_pos_))[ii] = (*((*gv_it)->ids_pos_))[jj];
 (*((*gv_it)->ids_pos_))[jj] = idx_buff;      

  if ( ( (*block_aops_rngs_)[ (*((*gv_it)->ids_pos_))[jj] ] == (*block_aops_rngs_)[ (*((*gv_it)->ids_pos_))[ii] ]) &&
       (*block_aops_)[ (*((*gv_it)->ids_pos_))[ii] ] != (*block_aops_)[(*((*gv_it)->ids_pos_))[jj] ] ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>( *((*gv_it)->deltas_pos_) );
    pair<int,int> new_delta = (*block_aops_)[ ( *((*gv_it)->ids_pos_))[jj] ]  ?
                               make_pair( (*((*gv_it)->ids_pos_))[jj], ( *((*gv_it)->ids_pos_))[ii] ): make_pair((*((*gv_it)->ids_pos_))[ii],(*((*gv_it)->ids_pos_))[jj]);

    new_deltas_tmp->push_back( new_delta );

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( (*gv_it)->ids_pos_->size()-2);
    vector<int>::iterator nip_it = new_ids_pos->begin();
    for( int qq = 0 ; qq != (*gv_it)->ids_pos_->size() ; qq++) {
      if ( (qq != ii) && (qq != jj)){
       *nip_it = (*((*gv_it)->ids_pos_))[qq];
        ++nip_it;
      }
    }

    shared_ptr<GammaIntermediateRedux> new_gamma = make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas_tmp, (*gv_it)->factors_ );

    print_vector( *((*gv_it)->ids_pos_) , " after_new_gamma (*gv_it)->ids_pos "); cout << endl;
    gamma_vec_->push_back(new_gamma);
    print_vector( *((*gv_it)->ids_pos_) , " after push_back (*gv_it)->ids_pos "); cout << endl;

  }
  pair<double,double> anti_herm_fac = make_pair( -1.0 , 1.0 );
  WickUtils::pair_fac_mult( anti_herm_fac, (*gv_it)->factors_  );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGeneratorRedux<DataType>::swap( int ii, int jj, int kk ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::swap ii = " << ii << " jj = " << jj << " kk = " << kk << endl;

  shared_ptr<GammaIntermediateRedux> gint =  gamma_vec_->at(kk); 

  int idx_buff = (*((gint)->ids_pos_))[ii];  
 (*((gint)->ids_pos_))[ii] = (*((gint)->ids_pos_))[jj];
 (*((gint)->ids_pos_))[jj] = idx_buff;      

  if ( ( (*block_aops_rngs_)[ (*((gint)->ids_pos_))[jj] ] == (*block_aops_rngs_)[ (*((gint)->ids_pos_))[ii] ]) &&
       (*block_aops_)[ (*((gint)->ids_pos_))[ii] ] != (*block_aops_)[(*((gint)->ids_pos_))[jj] ] ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>( *((gint)->deltas_pos_) );
    pair<int,int> new_delta = (*block_aops_)[ ( *((gint)->ids_pos_))[jj] ]  ?
                               make_pair( (*((gint)->ids_pos_))[jj], ( *((gint)->ids_pos_))[ii] ): make_pair((*((gint)->ids_pos_))[ii],(*((gint)->ids_pos_))[jj]);

    new_deltas_tmp->push_back( new_delta );

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( (gint)->ids_pos_->size()-2);
    vector<int>::iterator nip_it = new_ids_pos->begin();
    for( int qq = 0 ; qq != (gint)->ids_pos_->size() ; qq++) {
      if ( (qq != ii) && (qq != jj)){
       *nip_it = (*((gint)->ids_pos_))[qq];
        ++nip_it;
      }
    }

    shared_ptr<GammaIntermediateRedux> new_gamma = make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas_tmp, (gint)->factors_ );

    print_vector( *((gint)->ids_pos_) , " after_new_gamma (gint)->ids_pos "); cout << endl;
    gamma_vec_->push_back(new_gamma);
    print_vector( *(gamma_vec_->at(kk)->ids_pos_) , " after push_back (gamma_vec_->at(" +to_string(kk) +")->ids_pos) "); cout << endl;
    print_vector( *((gint)->ids_pos_) , " after push_back (gint)->ids_pos "); cout << endl;

  }
  pair<double,double> anti_herm_fac = make_pair( -1.0 , 1.0 );
  WickUtils::pair_fac_mult( anti_herm_fac, (gint)->factors_  );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
shared_ptr<pint_vec>
GammaGeneratorRedux<DataType>::standardize_delta_ordering_generic( shared_ptr<pint_vec> deltas_pos  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGeneratorRedux<DataType>::standardize_delta_ordering_generic" << endl;
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
template<typename DataType> 
vector<int> GammaGeneratorRedux<DataType>::get_standard_idx_order(const vector<string>&idxs) {
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
template<typename DataType> 
vector<int> GammaGeneratorRedux<DataType>::get_Aid_order ( const vector<int>& id_pos ) {
///////////////////////////////////////////////////////////////////////////////////////
//   cout << "GammaGeneratorRedux<DataType>::get_Aid_order " << endl;

  vector<int> new_id_pos(id_pos.size());
  vector<int> tmp_order = get_position_order(id_pos);
  vector<int> new_order = get_position_order(tmp_order) ;

  return new_order;
}
///////////////////////////////////////////////////////////////////////////////////////
//Returns ordering vector for ids_pos, e.g., if ids_pos = {6,4,5} , then pos = {1,2,0}
///////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
vector<int> GammaGeneratorRedux<DataType>::get_position_order(const vector<int> &ids_pos) {
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
template<typename DataType> 
vector<int> GammaGeneratorRedux<DataType>::get_standardized_alt_order ( const vector<string>& rngs ,const vector<bool>& aops ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux<DataType>::get_standardized_alt_order " << endl;
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
template<typename DataType> 
void GammaGeneratorRedux<DataType>::set_standardized_alt_order_unranged ( shared_ptr<GammaIntermediateRedux>& gint,
                                                                vector<int>& standard_alt_order ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::set_standardized_alt_order_unranged" << endl;
  // TODO this should use standardized ordering 

  shared_ptr<vector<int>>           ids_pos = gint->ids_pos_;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gint->deltas_pos_; 

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
template<typename DataType> 
void
GammaGeneratorRedux<DataType>::print_gamma_intermediate( shared_ptr<GammaIntermediateRedux> gint ) { 
/////////////////////////////////////////////////////////////////////////////////////////////
     print_vector( *(gint->ids_pos_), "gint_ids_pos" ); cout  << endl;
     cout << "gint_aops = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos_) ) { cout << block_aops_->at(pos) << " " ; cout.flush();  }  cout << "] " << endl;
     cout << "gint_rngs = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos_) ) { cout << (*block_aops_rngs_)[pos] << " " ; cout.flush();  }   cout << "]" << endl;
     cout << "gint_ids  = [ "; cout.flush();   for ( auto pos : *(gint->ids_pos_) ) { cout << block_idxs_[pos] << " " ; cout.flush();  }   cout << "] " <<  endl;
    return;
} 

 // { // TEST for transformations
 // vector<string> unique_block_dupe( std_rngs_.size());
 // vector<string>::const_iterator br_it = block_rngs_->begin();
 // for ( vector<int>::iterator it_it = idxs_trans_->begin() ;  it_it != idxs_trans_->end() ; it_it++, br_it++ ) { 
 //   unique_block_dupe[*it_it]  = *br_it;
 // }

 //  if ( unique_block_dupe != std_rngs_ ) { 
 //    print_vector(unique_block_dupe, "unique_block_dupe" ) ;
 //    print_vector(std_rngs_, " != std_rngs_" ) ;
 //    throw logic_error( " reordering is broken " ) ;
 //  }
 // }
//////////////////////////////////////////////////////////////////////////
template class GammaGeneratorRedux<double>;
template class GammaGeneratorRedux<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
