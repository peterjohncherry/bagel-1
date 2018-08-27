#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/debugging_utils.h>

//#define __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE
//#define __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_UNQ
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
  for ( string bra_name : *Bra_names_ ) {
   for ( string ket_name : *Ket_names_ ) {
      bra_name_ = bra_name;
      ket_name_ = ket_name;
      bra_hole_map_ = target_states_->hole_range_map(bra_name_);
      bra_elec_map_ = target_states_->elec_range_map(bra_name_);
      ket_hole_map_ = target_states_->hole_range_map(ket_name_);
      ket_elec_map_ = target_states_->elec_range_map(ket_name_);
      bool unq_does_it_contribute =  generic_reorderer_different_sector( reordering_name, final_reordering );
      if( unq_does_it_contribute ) 
        does_it_contribute = true;
    }
  }
  return does_it_contribute;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator_Base::generic_reorderer_different_sector( string reordering_name, bool final_reordering ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_UNQ
cout << "GammaGenerator_Base::generic_reorderer_different_sector : " << reordering_name ; cout.flush();
if (final_reordering)  { cout << " : final_reordering" << endl; } else { cout << " : intermediate_reordering" << endl; }
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( reordering_name == "normal order" ) {
    normal_order();

  } else if ( reordering_name == "anti-normal order" ) {
    anti_normal_order();

  } else if ( reordering_name == "alternating order" ) {
    alternating_order();
  }

  bool does_it_contribute = ( gamma_vec_.size() > 0 );
  if ( final_reordering && does_it_contribute ) { 
    int kk = 0;
    while ( kk != gamma_vec_.size()){
      add_Acontrib_to_map( kk, bra_name_, ket_name_ );
      ++kk;
    } 
  }

#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_UNQ
  cout << " ==================  UNIQUE_PTR VERSION WITH CONSTANT CHECKING =================== " << endl;
  cout << " ====  gammas left after transformation to " << reordering_name << " ===== " << endl;
  for ( const auto& gint : gamma_vec_ ) { print_gamma_intermediate( gint , "" ); }
#endif
  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::proj_onto_map( unique_ptr<GammaIntermediate_Base_Raw>& gint ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE_UNQ
cout << "GammaGenerator_Base::proj_onto_map" << endl;
WickUtils::print_vector( gint->ids_pos_ , "idxs_pos       " ); cout << endl;
Debugging_Utils::print_vector_at_pos(*block_aops_rngs_,  gint->ids_pos_, "aops_rngs"); cout << endl;
Debugging_Utils::print_vector_at_pos(*block_aops_,  gint->ids_pos_,      "aops_    "); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////

  assert( gint->survives_ == true ); // Should not be running this test if the gamma is already known to die.

  const vector<int>& idxs_pos =  gint->ids_pos_;
  map<char,int> ket_elec_map_mod = *ket_elec_map_;
  map<char,int> ket_hole_map_mod = *ket_hole_map_;
   for ( vector<int>::const_reverse_iterator ip_it = idxs_pos.rbegin(); ip_it !=  idxs_pos.rend(); ++ip_it ) {
    char rng = (*block_aops_rngs_)[*ip_it];  

    if( !((*block_aops_)[*ip_it]) ){
      auto ket_elec_map_mod_loc = ket_elec_map_mod.find( rng );
      if ( ket_elec_map_mod_loc == ket_elec_map_mod.end() ) {
        gint->survives_ = false;
        return;
      } else if ( (ket_elec_map_mod_loc->second -= 1 ) == -1  ) {
        gint->survives_ = false;
        return;
      }
      auto ket_hole_map_mod_loc = ket_hole_map_mod.find( rng );

      if ( ket_hole_map_mod_loc == ket_hole_map_mod.end() ) {
        ket_hole_map_mod.emplace( rng, 1 );
      } else {
        ket_hole_map_mod_loc->second += 1;
      }

    } else {
      auto ket_hole_map_mod_loc = ket_hole_map_mod.find( rng );
      if ( ket_hole_map_mod_loc == ket_hole_map_mod.end() ) {
        gint->survives_ = false;
        return;
      } else if ( (ket_hole_map_mod_loc->second -= 1 ) == -1  ) {
        gint->survives_ = false;
        return;
      }
      
      auto ket_elec_map_mod_loc = ket_elec_map_mod.find( rng );
      if ( ket_elec_map_mod_loc == ket_elec_map_mod.end() ) {
        ket_elec_map_mod.emplace( rng, 1 );
      } else {
        ket_elec_map_mod_loc->second += 1;
      }
    }
  }
  
  // TODO Will break if bra and ket but built from different orbitals ( not just differing occupations, but different sets ).
  // Compares modified ket_map to original bra map
  for ( auto& elem : ket_elec_map_mod )  
    if ( elem.second != bra_elec_map_->at(elem.first) ){
      gint->survives_ = false;
      return;     
    }

  for ( auto& elem : ket_hole_map_mod )                     
    if ( elem.second != bra_hole_map_->at(elem.first) ){
      gint->survives_ = false;
      return;     
    }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator_Base::proj_onto_map( const unique_ptr<GammaIntermediate_Base_Raw>& gint, 
                                             map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                             map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE_UNQ
cout << "GammaGenerator_Base::proj_onto_map" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////

  const vector<int>& idxs_pos =  gint->ids_pos_;
  for ( vector<int>::const_reverse_iterator ip_it = idxs_pos.rbegin(); ip_it !=  idxs_pos.rend(); ++ip_it ) {
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
  
  // TODO Will break if bra and ket but built from different orbitals ( not just differing occupations, but different sets ).
  // Compares modified ket_map to original bra map
  for ( auto& elem : ket_elec_map )  
    if ( elem.second != bra_elec_map_->at(elem.first) ){
      return false;     
    }

  for ( auto& elem : ket_hole_map )                     
    if ( elem.second != bra_hole_map_->at(elem.first) ){
      return false;     
    }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE_UNQ
cout << "GammaGenerator_Base::normal_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
  int kk = 0;
  int number_of_surviving_gammas = 0;
  vector<unique_ptr<GammaIntermediate_Base_Raw>> surviving_gammas;
  while ( kk != gamma_vec_.size() ) {
 
    proj_onto_map( gamma_vec_[kk] ); 

    if ( gamma_vec_[kk]->survives_ ) { 
      vector<int>& ids_pos = gamma_vec_[kk]->ids_pos_;
      int num_plus = 0;
      for ( int pos : ids_pos )
        if ( (*block_aops_)[ pos ]) 
          num_plus++; 
      
      num_plus--; 
     
      for (int ii = ids_pos.size()-1 ; ii != -1; ii--){
      
        if ( ii > num_plus ) {
          while( (*block_aops_)[ids_pos[ii]] &&  gamma_vec_[kk]->survives_ ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if ( !(*block_aops_)[ids_pos[jj]] ){
                swap( jj, jj+1, kk);
                proj_onto_map( gamma_vec_[kk] );
                break;
              }
            }
          }
          if ( !(gamma_vec_[kk]->survives_) )break;
      
        } else if (ii <= num_plus) {
          if ( (*block_aops_)[ids_pos[ii]])
            continue;
      
          while(!(*block_aops_)[ids_pos[ii]] && gamma_vec_[kk]->survives_ ){
            for ( int jj = (ii-1); jj != -1 ; jj--) {
              if((*block_aops_)[ids_pos[jj]] ){
                swap( jj, jj+1, kk);
                proj_onto_map( gamma_vec_[kk] );
                break;
              }
            }
          }
        }
        if ( !(gamma_vec_[kk]->survives_) )break;
      }
    }
    if ( gamma_vec_[kk]->survives_ ){
      ++number_of_surviving_gammas;
       //print_gamma_intermediate (gamma_vec_[kk] , "SURVIVING GAMMA !! NORMAL ORDER !!" );
    }
    ++kk;
  }
  if ( number_of_surviving_gammas > 0 ){ 
    vector<unique_ptr<GammaIntermediate_Base_Raw>> surviving_gammas(number_of_surviving_gammas);
    vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator sg_it = surviving_gammas.begin();  
    for ( vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator gv_it = gamma_vec_.begin(); gv_it != gamma_vec_.end(); gv_it++){ 
      if ( (*gv_it)->survives_ ) {
        *sg_it = make_unique<GammaIntermediate_Base_Raw>(*(gv_it->get())); // TODO can this be move, or would that screw up the iterators in gvu?
        ++sg_it;
      }
    }
    gamma_vec_ = move(surviving_gammas);
  } else { 
    gamma_vec_.clear();
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::anti_normal_order() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE_UNQ
cout << "GammaGenerator_Base::anti_normal_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
  int kk = 0;
  int number_of_surviving_gammas = 0;
  while ( kk != gamma_vec_.size() ) {
    proj_onto_map( gamma_vec_[kk] ); 

    if ( gamma_vec_[kk]->survives_ ) { 
      vector<int>& ids_pos  = gamma_vec_[kk]->ids_pos_;
      int num_kill = 0;
      for ( int pos : ids_pos )
        if (!(*block_aops_)[pos]) 
          num_kill++; 
      num_kill--; 
      
      for (int ii = ids_pos.size()-1 ; ii != -1; ii--){
        if ( ii > num_kill ) {
          while( ((!(*block_aops_)[ids_pos[ii]]) && (gamma_vec_[kk]->survives_)) ) {
            for ( int jj = ii-1; jj != -1 ; jj--) {
              if ( (*block_aops_)[ids_pos[jj]]  ){
                swap( jj, jj+1, kk );
                proj_onto_map( gamma_vec_[kk] );
                break;
              }
            }
          }

        } else if (ii <= num_kill) {
          while((*block_aops_)[ids_pos[ii]]  && gamma_vec_[kk]->survives_){
            for ( int jj = (ii-1); (jj != -1 && gamma_vec_[kk]->survives_); jj--) {
              if( !(*block_aops_)[ids_pos[jj]] ) {
                swap( jj, jj+1, kk );
                proj_onto_map( gamma_vec_[kk] ); 
                break;
              }
            }
          }
        }
        if (!gamma_vec_[kk]->survives_ ) { /* print_gamma_intermediate( gamma_vec_[kk], "Gamma dies in anti-normal"); */ break;}
      }
    } 
    if ( gamma_vec_[kk]->survives_ ){
      ++number_of_surviving_gammas;
    }
    ++kk; 
  }

  cout << "number of surviving gamma = " << number_of_surviving_gammas << endl;
  if ( number_of_surviving_gammas > 0 ){ 
    vector<unique_ptr<GammaIntermediate_Base_Raw>> surviving_gammas(number_of_surviving_gammas);
    vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator sg_it = surviving_gammas.begin();  
    for ( vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator gv_it = gamma_vec_.begin(); gv_it != gamma_vec_.end(); gv_it++){ 
      if ( (*gv_it)->survives_ ) {
        *sg_it = make_unique<GammaIntermediate_Base_Raw>(*(gv_it->get())); // TODO can this be move, or would that screw up the iterators in gvu?
        ++sg_it;
      }
    }
    gamma_vec_ = move(surviving_gammas);
  } else { 
    gamma_vec_.clear();
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::alternating_order() {  // e.g. +-+-+-+-
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE_UNQ
cout << "GammaGenerator_Base::alternating_order" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////
 
  int kk = 0;
  int number_of_surviving_gammas = 0;
  while ( kk != gamma_vec_.size() ) {
    proj_onto_map( gamma_vec_[kk] ); 
    if ( gamma_vec_[kk]->survives_ ){ 
      vector<int> ids_pos = gamma_vec_[kk]->ids_pos_; 
      vector<int> new_ids_pos = get_standardized_alt_order_unranged( ids_pos );

      for ( int ii = ids_pos.size()-1; ((ii!=-1) && (gamma_vec_[kk]->survives_)); --ii ){
        if ( ids_pos[ii] == new_ids_pos[ii])
          continue;

        while( ids_pos[ii] != new_ids_pos[ii] && gamma_vec_[kk]->survives_ ){
          for ( int jj = ii-1; jj != -1 ; jj--) {
            if ( ids_pos[jj] == new_ids_pos[ii] ){
              swap( jj, jj+1, kk );
              proj_onto_map( gamma_vec_[kk] );
              break;
            }
          }
        }
      }
    }
    if (gamma_vec_[kk]->survives_ )
      ++number_of_surviving_gammas;
    
    ++kk;
  }
  if ( number_of_surviving_gammas > 0 ){
    vector<unique_ptr<GammaIntermediate_Base_Raw>> surviving_gammas(number_of_surviving_gammas);
    vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator sg_it = surviving_gammas.begin();  
    for ( vector<unique_ptr<GammaIntermediate_Base_Raw>>::iterator gv_it = gamma_vec_.begin(); gv_it != gamma_vec_.end(); gv_it++){ 
      if ( (*gv_it)->survives_ ) {
        *sg_it = make_unique<GammaIntermediate_Base_Raw>(*(gv_it->get())); // TODO can this be move, or would that screw up the iterators in gvu?
        ++sg_it;
      }
    }
    gamma_vec_ = move(surviving_gammas);
  } else {
    gamma_vec_.clear();
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr ( why ???? ).
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::swap( int ii, int jj, int kk ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP_UNQ
cout << "GammaGenerator_Base::swap ii = " << ii << " jj = " << jj << " kk = " << kk << endl;
print_gamma_intermediate(gamma_vec_[kk] , "pre_swap gamma" );
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  unique_ptr<GammaIntermediate_Base_Raw>& gint =  gamma_vec_[kk];

  vector<int>& ids_pos_kk  =  gamma_vec_[kk]->ids_pos_;
  // using an extra buffer variable, but makes routine in if clearer 
  int i_pos = ids_pos_kk[ii];
  int j_pos = ids_pos_kk[jj];
 
  ids_pos_kk[ii] = j_pos;
  ids_pos_kk[jj] = i_pos;

  if ( ( (*block_aops_rngs_)[j_pos] == (*block_aops_rngs_)[i_pos]) && (*block_aops_)[i_pos] != (*block_aops_)[ j_pos ] ){

    vector<pair<int, int>>& deltas_pos_kk  =  gamma_vec_[kk]->deltas_pos_;
    vector<pair<int,int>> new_deltas_tmp( deltas_pos_kk.size()+1);
    copy ( deltas_pos_kk.begin(), deltas_pos_kk.end(), new_deltas_tmp.begin());
    new_deltas_tmp.back() = (*block_aops_)[ j_pos ]  ?  make_pair( j_pos, i_pos ): make_pair( i_pos, j_pos);

    vector<int> new_ids_pos(ids_pos_kk.size()-2);
    {
    vector<int>::iterator gip_it =  ( ii < jj ) ? ids_pos_kk.begin()+ii : ids_pos_kk.begin()+jj ; 
    copy( gip_it+2 , ids_pos_kk.end(), copy( ids_pos_kk.begin(), gip_it, new_ids_pos.begin() ));
    } 

    auto new_fac =  make_pair( gamma_vec_[kk]->factors_.first * 1.0,  0.0 );
    gamma_vec_.push_back(make_unique<GammaIntermediate_Base_Raw>( new_ids_pos, new_deltas_tmp, new_fac ));
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP_UNQ
    print_gamma_intermediate( gamma_vec_.back() , "new gamma" ); cout << endl;
#endif

  } 
  gamma_vec_[kk]->factors_ =  make_pair( gamma_vec_[kk]->factors_.first * -1.0,  0.0 );

#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_SWAP_UNQ
  print_gamma_intermediate( gamma_vec_[kk] , "post_swap gamma" );
#endif
  return;
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
vector<int> GammaGenerator_Base::get_position_order(const vector<int>& ids_pos) {
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
vector<int> GammaGenerator_Base::get_standardized_alt_order_unranged( vector<int> ids_reordered_pos ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_VERBOSE_GAMMAGENERATOR_BASE_PROPTOOL
cout << "GammaGenerator_Base::get_standardized_alt_order_unranged" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  vector<int> standard_alt_order( ids_reordered_pos.size() ); 
  vector<int>::iterator sao_it = standard_alt_order.begin();
  while ( sao_it != standard_alt_order.end()){
    *sao_it = *sop_it;
    ++sop_it;
    ++sao_it;
    *sao_it = *sok_it;
    ++sok_it;
    ++sao_it;
  }

  return standard_alt_order;
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
  vector<pair<int,int>>::iterator  b_it = buff.begin();
  for ( vector<pair<int,int>>::iterator dp_it = deltas_pos.begin(); dp_it != deltas_pos.end(); dp_it++, b_it++ ) 
    *dp_it = make_pair( standard_order_[b_it->first], standard_order_[b_it->second]);

  return;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GammaGenerator_Base::print_gamma_intermediate( const unique_ptr<GammaIntermediate_Base_Raw>& gint, string gamma_name ) { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMAGENERATOR_BASE_VERBOSE
cout << "GammaGenerator_Base::print_gamma_intermediate (unique_ptr_version)" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( gamma_name != "" ) 
    cout << "-------- " << gamma_name << " -------" << endl;

  WickUtils::print_vector(block_idxs_,                                     "block_idxs     "); cout << endl;
  WickUtils::print_vector(gint->ids_pos_,                                  "gint_ids_pos   "); cout << endl;
  WickUtils::print_pair_vector( gint->deltas_pos_,                         "gint_deltas_pos"); cout << endl;
  Debugging_Utils::print_vector_at_pos( block_idxs_,  gint->ids_pos_,      "gint_idxs      "); cout << endl;

  cout << "gint_deltas_idxs  = [ ";  cout.flush();
  for(const auto& dta : gint->deltas_pos_){
    cout << "("<< block_idxs_[dta.first]  << ":" << (*block_aops_rngs_)[dta.first ] << ","; cout.flush();
    cout <<       block_idxs_[dta.second] << ":" << (*block_aops_rngs_)[dta.second] << ") "; cout.flush();
  }
  cout << "] " <<  endl;

  Debugging_Utils::print_vector_at_pos(*block_aops_,       gint->ids_pos_,      "gint_aops      "); cout << endl;
  Debugging_Utils::print_vector_at_pos(*block_aops_rngs_,  gint->ids_pos_, "gint_aops_rngs "); cout << endl;

  cout << "factor  = ( " ; cout.flush(); cout << gint->factors_.first << ", " << gint->factors_.second << " )" << endl << endl;
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator_Base::transformation_tester( unique_ptr<GammaIntermediate_Base_Raw>& gint  ){  // e.g. ++++----
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

  cout << endl << endl << endl;
  cout << "----------- gamma orb deriv Contribs -------------------" << endl;
  print_vector( gint->ids_pos_ ,                         " ids_pos     _           "); cout <<endl;
  print_pair_vector( gint->deltas_pos_,                 " A_A_deltas_pos           "); cout << endl;
  print_pair_vector( gint->target_target_deltas_pos(),  " target_target_deltas_pos "); cout << endl;
  print_pair_vector( gint->target_A_deltas_pos(),       " target_A_deltas_pos      "); cout << endl;
  cout << endl;

  if ( gint->A_A_deltas_pos().size() != 0 ) { 

    vector<pair<string,string>> A_A_deltas_rngs;
    for ( pair<int,int>& ctr :  gint->deltas_pos_ )
      A_A_deltas_rngs.push_back(make_pair((*block_rngs_)[ctr.first],(*block_rngs_)[ctr.second])) ;
    
    vector<pair<string,string>> A_A_deltas_idxs;
    for ( pair<int,int>& ctr :  gint->deltas_pos_ ) 
      A_A_deltas_idxs.push_back(make_pair(block_idxs_[ctr.first] , block_idxs_[ctr.second] )) ;


    vector<pair<int,int>> A_A_deltas_pos( gint->deltas_pos_.size() ); 
    vector<pair<int,int>>::iterator aadp_it =  A_A_deltas_pos.begin(); 
    for ( auto gaadp_it = gint->deltas_pos_.begin(); gaadp_it != gint->deltas_pos_.end(); aadp_it++, gaadp_it++ ){
       (*aadp_it).first  =  (*idxs_trans_)[ gaadp_it->first  ];
       (*aadp_it).second =  (*idxs_trans_)[ gaadp_it->second ];
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

  if ( gint->target_target_deltas_pos().size() != 0 ) { 

    vector<pair<string,string>> T_T_deltas_rngs;
    for ( pair<int,int>& ctr :  gint->target_target_deltas_pos() )
       T_T_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_idxs;
    for ( pair<int,int>& ctr :  gint->target_target_deltas_pos() )
       T_T_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr : gint->target_target_deltas_pos() ){
       ctr.first  =  (*idxs_trans_)[ctr.first];
       ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_T_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  gint->target_target_deltas_pos() )
       T_T_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  gint->target_target_deltas_pos() )
       T_T_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );
  
    print_pair_vector( T_T_deltas_trans_idxs , "T_T_deltas_trans_idxs " ); 
    print_pair_vector( T_T_deltas_idxs , "   T_T_deltas_idxs " ); 

    assert( T_T_deltas_trans_idxs == T_T_deltas_idxs ) ;
  }

  if ( gint->target_A_deltas_pos().size() != 0 ) { 
 
    vector<pair<string,string>> T_A_deltas_rngs;
    for ( pair<int,int>& ctr :  gint->target_A_deltas_pos() )
       T_A_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_idxs;
    for ( pair<int,int>& ctr : gint->target_A_deltas_pos() )
       T_A_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr :  gint->target_A_deltas_pos() ){
      ctr.first  =  (*idxs_trans_)[ctr.first];
      ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_A_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  gint->target_A_deltas_pos() )
       T_A_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  gint->target_A_deltas_pos() )
       T_A_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    print_pair_vector( T_A_deltas_trans_idxs , "T_A_deltas_trans_idxs " ); 
    print_pair_vector( T_A_deltas_idxs , "   T_A_deltas_idxs " ); 

    assert( T_A_deltas_trans_idxs == T_A_deltas_idxs )  ;
 
  } 
  cout << endl << endl;
  return;
}
