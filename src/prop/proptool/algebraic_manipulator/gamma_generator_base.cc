#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>
#include <src/prop/proptool/proputils.h>

#define __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE
#define __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP
using namespace std;
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator_Base::generic_reorderer( string reordering_name, bool first_reordering, bool final_reordering ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE
cout << "GammaGenerator_Base::generic_reorderer" << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int kk = 0;
  bool does_it_contribute = false;
  final_gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediate_Base>>>(0);
  for ( string bra_name : *Bra_names_ ) {
    for ( string ket_name : *Ket_names_ ) {
      bra_name_ = bra_name;
      ket_name_ = ket_name;
      bra_hole_map_ = target_states_->hole_range_map(bra_name_);
      bra_elec_map_ = target_states_->elec_range_map(bra_name_);
      ket_hole_map_ = target_states_->hole_range_map(ket_name_);
      ket_elec_map_ = target_states_->elec_range_map(ket_name_);
      if( generic_reorderer_different_sector( reordering_name, final_reordering ) ) //TODO find a better way; have contribs specified for individual BraKet combinations?
        does_it_contribute = true;
    }
  }
  return does_it_contribute;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator_Base::generic_reorderer_different_sector( string reordering_name, bool final_reordering ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::generic_reorderer_different_sector" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( reordering_name == "normal order" ) {
    normal_order();
    typename vector<shared_ptr<GammaIntermediate_Base >>::iterator gv_it = gamma_vec_->begin();
    while ( gv_it != gamma_vec_->end() ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) ) 
        final_gamma_vec_->push_back( *gv_it );
      gv_it++;
    }

  } else if ( reordering_name == "anti-normal order" ) {
    anti_normal_order();

    typename vector<shared_ptr<GammaIntermediate_Base>>::iterator gv_it = gamma_vec_->begin();
    while ( gv_it != gamma_vec_->end() ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) )
        final_gamma_vec_->push_back( *gv_it );
      ++gv_it;
    }

  } else if ( reordering_name == "alternating order" ) {
    alternating_order();
    typename vector<shared_ptr<GammaIntermediate_Base>>::iterator gv_it = gamma_vec_->begin();
    while ( gv_it != gamma_vec_->end() ) {
      if ( proj_onto_map( *gv_it, *bra_hole_map_, *bra_elec_map_, *ket_hole_map_, *ket_elec_map_ ) )
        final_gamma_vec_->push_back( *gv_it );
      ++gv_it;
    }
  }
  cout << " ====  gammas left after transformation to " << reordering_name << " ===== " << endl;
  for ( const auto& gint : *final_gamma_vec_ ) { print_gamma_intermediate( gint , "" ); }

  gamma_vec_ = final_gamma_vec_;
  bool does_it_contribute = ( gamma_vec_->size() > 0 );

  int kk = 0;
  if ( final_reordering && does_it_contribute ) { 
    while ( kk != gamma_vec_->size()){
      add_Acontrib_to_map( kk, bra_name_, ket_name_ );
      ++kk;
    } 
  }
  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator_Base::proj_onto_map( shared_ptr<GammaIntermediate_Base> gint, 
                                         map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                         map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::proj_onto_map" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<int>> idxs_pos =  gint->ids_pos_;

  for ( vector<int>::reverse_iterator ip_it = idxs_pos->rbegin(); ip_it !=  idxs_pos->rend(); ip_it++ ) {
    char rng = (*block_aops_rngs_)[*ip_it];  

    if( !((*block_aops_)[*ip_it]) ){
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

  // TODO  may have problems if bra and ket involve different ranges, but not relevant for now.
  for ( auto& elem : *ket_elec_map_ )  
    if ( elem.second != bra_elec_map.at(elem.first) )
      return false;     
    

  for ( auto& elem : *ket_hole_map_ )                     
    if ( elem.second != bra_hole_map.at(elem.first) )
      return false;     
    
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::normal_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
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
void GammaGenerator_Base::anti_normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::anti_normal_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
  int kk = 0;
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
                swap( jj, jj+1, kk );
                break;
              }
            }
          }
      
        } else if (ii <= num_kill) {
          if (!(*block_aops_)[(*ids_pos)[ii]] )
            continue;
      
          while((*block_aops_)[(*ids_pos)[ii]] ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if(!(*block_aops_)[(*ids_pos)[jj]] ) {
                swap( jj, jj+1, kk );
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
void GammaGenerator_Base::alternating_order() {  // e.g. +-+-+-+-
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::alternating_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
  int kk = 0;
  while ( kk != gamma_vec_->size() ) {
    shared_ptr<GammaIntermediate_Base> gint = gamma_vec_->at(kk);
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
    ++kk;
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::swap( int ii, int jj, int kk ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP
cout << "GammaGenerator_Base::swap ii = " << ii << " jj = " << jj << " kk = " << kk << endl;
print_gamma_intermediate (gamma_vec_->at(kk) , "pre_swap gamma" );
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<GammaIntermediate_Base> gint =  gamma_vec_->at(kk);
 
  int i_pos = (*(gint->ids_pos_))[ii];
  int j_pos = (*(gint->ids_pos_))[jj];
 
  (*(gint->ids_pos_))[ii] = j_pos;
  (*(gint->ids_pos_))[jj] = i_pos;

  if ( ( (*block_aops_rngs_)[j_pos] == (*block_aops_rngs_)[i_pos]) && (*block_aops_)[i_pos] != (*block_aops_)[ j_pos ] ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>( *(gint->deltas_pos_) );
    pair<int,int> new_delta = (*block_aops_)[ j_pos ]  ?  make_pair( j_pos, i_pos ): make_pair( i_pos, j_pos);

    new_deltas_tmp->push_back( new_delta );

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( gint->ids_pos_->size()-2);
    vector<int>::iterator nip_it = new_ids_pos->begin();
    vector<int>::iterator gip_it = gint->ids_pos_->begin();
    for( int qq = 0 ; qq != gint->ids_pos_->size(); ++qq, ++gip_it ) {
      if ( (qq != ii) && (qq != jj)){
       *nip_it = *gip_it;
        ++nip_it;
      }
    }
 
    auto new_fac =  make_pair( gint->factors_.first * -1.0,  0.0 );
    shared_ptr<GammaIntermediate_Base> new_gamma = make_shared<GammaIntermediate_Base>( new_ids_pos, new_deltas_tmp, new_fac );
    gamma_vec_->push_back(new_gamma);

#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP
    print_gamma_intermediate (new_gamma , "new gamma" );
  }
  print_gamma_intermediate (gamma_vec_->at(kk) , "post_swap gamma" );
#else
  }
#endif
 
  return;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec>
GammaGenerator_Base::standardize_delta_ordering_generic( shared_ptr<pint_vec> deltas_pos  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::standardize_delta_ordering_generic" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////////
// Getting reordering vec to go from Atens uncids to gamma uncids
// Running sort twice to get inverse; seems weird and there's probably a better way...
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator_Base::get_Aid_order ( const vector<int>& id_pos ) {
///////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::get_Aid_order " << endl;
#endif ////////////////////////////////////////////////////////////////////////////////

  vector<int> new_id_pos(id_pos.size());
  vector<int> tmp_order = get_position_order(id_pos);
  vector<int> new_order = get_position_order(tmp_order) ;

  return new_order;
}
///////////////////////////////////////////////////////////////////////////////////////
//Returns ordering vector for ids_pos, e.g., if ids_pos = {6,4,5} , then pos = {1,2,0}
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator_Base::get_position_order(const vector<int> &ids_pos) {
///////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::get_position_order" << endl;
#endif ////////////////////////////////////////////////////////////////////////////////
  vector<int> pos(ids_pos.size());
  iota(pos.begin(), pos.end(), 0);
  sort(pos.begin(), pos.end(), [&ids_pos](int i1, int i2){return (bool)( ids_pos[i1] < ids_pos[i2] ); });

  return pos;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This final ordering is preferable for orbital excitation derivative terms; it ensures that the indexes of the
// perturbation tensor are always in the same order, which makes combining terms easier. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::set_standardized_alt_order_unranged ( shared_ptr<GammaIntermediate_Base>& gint,
                                                                          vector<int>& standard_alt_order ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_VERBOSE_GAMMAGENERATOR_BASE_PROPTOOL
cout << "GammaGenerator_Base::set_standardized_alt_order_unranged" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> ids_reordered_pos =  *(gint->ids_pos_);
  { 
  vector<int>& tmp = standard_order_; // TODO cannot remember how to avoid this without asking lambda to capture whole class...
  sort( ids_reordered_pos.begin(), ids_reordered_pos.end(), [tmp] ( int i1, int i2) { return (bool)( tmp[i1] < tmp[i2]); });  
  } 
 
  vector<int> standard_order_plus(ids_reordered_pos.size()/2);
  vector<int> standard_order_kill(ids_reordered_pos.size()/2);
  vector<int>::iterator sop_it = standard_order_plus.begin();
  vector<int>::iterator sok_it = standard_order_kill.begin();
  for ( int pos : ids_reordered_pos) {
    if ( block_aops_->at(pos) ) {
      *sop_it = pos;
      ++sop_it;
    } else {
      *sok_it = pos;
      ++sok_it;
    }
  }

  sop_it = standard_order_plus.begin();
  sok_it = standard_order_kill.begin();
  vector<int>::iterator sao_it = standard_alt_order.begin();
  while ( sao_it != standard_alt_order.end()){
    *sao_it = *sop_it;
    ++sop_it;
    ++sao_it;
    *sao_it = *sok_it;
    ++sok_it;
    ++sao_it;
  }

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::transform_to_canonical_ids_pos( vector<int>& ids_pos ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::transform_to_canonical_ids_pos" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  vector<int> buff = ids_pos;
  vector<int>::iterator b_it = buff.begin();
  for ( vector<int>::iterator ip_it = ids_pos.begin(); ip_it != ids_pos.end() ; ip_it++, b_it++ )
    *ip_it = standard_order_[*b_it];

  return;
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::transform_to_canonical_ids_pos( vector<pair<int,int>>& deltas_pos ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::transform_to_canonical_ids_pos" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<pair<int,int>> buff = deltas_pos;
  vector<pair<int,int>>::iterator  b_it = deltas_pos.begin();
  for ( vector<pair<int,int>>::iterator dp_it = deltas_pos.begin(); dp_it != deltas_pos.end(); dp_it++, b_it++ ) 
    *b_it = make_pair( standard_order_[dp_it->first], standard_order_[dp_it->second]);

  return;
} 
/////////////////////////////////////////////////////////////////////////////////////////////
void
GammaGenerator_Base::print_gamma_intermediate( const shared_ptr<GammaIntermediate_Base>& gint, string gamma_name ) { 
/////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::print_gamma_intermediate" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////

  if ( gamma_name != "" ) 
    cout << "-------- " << gamma_name << " -------" << endl;

  cout << "gint_ids_pos = [ ";
  cout.flush(); for(const auto& pos : *(gint->ids_pos_)){ cout << pos << "  " ; cout.flush();  }  cout << "] " << endl;
  cout << "gint_aops    = [ ";
  cout.flush(); for(const auto& pos : *(gint->ids_pos_)){ cout << block_aops_->at(pos) << "  " ; cout.flush();  }  cout << "] " << endl;
  cout << "gint_rngs    = [ ";
  cout.flush(); for(const auto& pos : *(gint->ids_pos_)){ cout << (*block_aops_rngs_)[pos] << "  " ; cout.flush();  }   cout << "]" << endl;
  cout << "gint_ids     = [ ";
  cout.flush(); for(const auto& pos : *(gint->ids_pos_)){ cout << block_idxs_[pos] << " " ; cout.flush();  }   cout << "] " <<  endl;
  cout << "gint_deltas_pos   = [ ";
  cout.flush(); for(const auto& dta : *(gint->deltas_pos_)){ cout << "("<< dta.first << "," << dta.second << ")"; cout.flush();  }cout << "] " <<  endl;
  cout << "gint_deltas_idxs  = [ ";
  cout.flush(); for(const auto& dta : *(gint->deltas_pos_)){ cout << "("<< block_idxs_[dta.first] << "," << block_idxs_[dta.second] << ")"; cout.flush();  }cout << "] " <<  endl;
  cout << "factor  = ( " ; cout.flush(); cout << gint->factors_.first << ", " << gint->factors_.second << " )" << endl << endl;
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::transformation_tester( shared_ptr<GammaIntermediate_Base>& gint  ){  // e.g. ++++----
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE
cout << "GammaGenerator_Base::transformation_tester" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  block_idxs_ = vector<string>( std_ids_->size() );
  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ )
    *bi_it = (*std_ids_)[ *it_it ];
  }

  shared_ptr<vector<pair<int,int>>> deltas_pos = gint->deltas_pos_;

  cout << endl << endl << endl;
  cout << "----------- gamma orb deriv Contribs -------------------" << endl;
  print_vector( *(gint->ids_pos_) ,                         " ids_pos     _            "); cout <<endl;
  print_pair_vector( *( gint->deltas_pos_),                 " A_A_deltas_pos           "); cout << endl;
  print_pair_vector( *( gint->target_target_deltas_pos()),  " target_target_deltas_pos "); cout << endl;
  print_pair_vector( *( gint->target_A_deltas_pos()),       " target_A_deltas_pos      "); cout << endl;
  cout << endl;

  if ( gint->A_A_deltas_pos()->size() != 0 ) { 

    vector<pair<string,string>> A_A_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->deltas_pos_) )
      A_A_deltas_rngs.push_back(make_pair((*block_rngs_)[ctr.first],(*block_rngs_)[ctr.second])) ;
    
    vector<pair<string,string>> A_A_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->deltas_pos_) ) 
      A_A_deltas_idxs.push_back(make_pair(block_idxs_[ctr.first] , block_idxs_[ctr.second] )) ;


    vector<pair<int,int>> A_A_deltas_pos( gint->deltas_pos_->size() ); 
    vector<pair<int,int>>::iterator aadp_it =  A_A_deltas_pos.begin(); 
    for ( auto gaadp_it = gint->deltas_pos_->begin(); gaadp_it != gint->deltas_pos_->end(); aadp_it++, gaadp_it++ ){
       (*aadp_it).first  =  (*idxs_trans_)[ (*gaadp_it).first  ];
       (*aadp_it).second =  (*idxs_trans_)[ (*gaadp_it).second ];
    }

    vector<pair<string,string>> A_A_deltas_rngs_trans;
    for ( pair<int,int>& ctr :  A_A_deltas_pos )
      A_A_deltas_rngs_trans.push_back(make_pair(std_rngs_[ctr.first] ,std_rngs_[ctr.second])) ;
    
    vector<pair<string,string>> A_A_deltas_idxs_trans;
    for ( pair<int,int>& ctr :  A_A_deltas_pos ) 
      A_A_deltas_idxs_trans.push_back(make_pair((*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] )) ;

    print_pair_vector( A_A_deltas_idxs_trans , "A_A_deltas_idxs_trans " ); 
    print_pair_vector( A_A_deltas_idxs , "   A_A_deltas_idxs " ); 

    assert( A_A_deltas_idxs_trans == A_A_deltas_idxs ); 

  } 

  if ( gint->target_target_deltas_pos()->size() != 0 ) { 

    vector<pair<string,string>> T_T_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) ){
       ctr.first  =  (*idxs_trans_)[ctr.first];
       ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_T_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );
  
    print_pair_vector( T_T_deltas_trans_idxs , "T_T_deltas_trans_idxs " ); 
    print_pair_vector( T_T_deltas_idxs , "   T_T_deltas_idxs " ); 

    assert( T_T_deltas_trans_idxs == T_T_deltas_idxs ) ;
  }

  if ( gint->target_A_deltas_pos()->size() != 0 ) { 
 
    vector<pair<string,string>> T_A_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) ){
      ctr.first  =  (*idxs_trans_)[ctr.first];
      ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_A_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    print_pair_vector( T_A_deltas_trans_idxs , "T_A_deltas_trans_idxs " ); 
    print_pair_vector( T_A_deltas_idxs , "   T_A_deltas_idxs " ); 

    assert( T_A_deltas_trans_idxs == T_A_deltas_idxs )  ;
 
  } 
  cout << endl << endl;
  return;
} 
