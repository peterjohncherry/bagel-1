#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaGeneratorRedux::GammaGeneratorRedux( shared_ptr<StatesInfo<double>> target_states, int Bra_num, int Ket_num,
                                          shared_ptr<const vector<string>> orig_ids, shared_ptr<const vector<bool>> orig_aops,
                                          shared_ptr<map<string, shared_ptr<GammaInfo>>>& Gamma_map_in,
                                          shared_ptr<map<string, shared_ptr<map<string, shared_ptr<AContribInfo> >>>>& G_to_A_map_in,
                                          double bk_factor                                                           ):
                                          target_states_(target_states),
                                          Bra_names_(target_states_->civec_names( Bra_num )),
                                          Ket_names_(target_states_->civec_names( Ket_num )),
                                          orig_ids_(orig_ids), orig_aops_(orig_aops),
                                          G_to_A_map(G_to_A_map_in), Gamma_map(Gamma_map_in), bk_factor_(bk_factor),
                                          orig_aops_half_size_( orig_aops_->size()/2 ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::GammaGeneratorRedux" << endl;

  print_vector(*orig_aops_, "orig_aops"); cout << endl; // This should be constant for all range blocks, but this is not the same as the MT aops
  print_vector(*orig_ids_, "orig_ids"); cout << endl; // ids are not (necessarily) be constant for all range blocks
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::add_gamma( shared_ptr<Range_BlockX_Info> block_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::add_gamma" << endl;

  orig_aops_ = block_info->trans_aops(); 

  standardized_full_id_ranges_ = *(block_info->orig_rngs()); 
  standardized_full_ids_ = *(block_info->orig_idxs()); 

  block_aops_ = block_info->trans_aops(); 
  block_rngs_ = block_info->trans_rngs(); 
  block_idxs_ = block_info->trans_idxs(); 
  
  standard_order_ = *(block_info->idxs_trans());

  int ii = 0 ;
  for ( vector<int>::iterator so_it = standard_order_.begin() ; so_it != standard_order_.end() ; ++so_it, ++ii ) 
    block_to_std_order_[*so_it] = (ii);

  shared_ptr<vector<int>> ids_pos = make_shared<vector<int>>( orig_ids_->size() );
  iota( ids_pos->begin(), ids_pos->end(), 0 );

  shared_ptr<vector<pair<int,int>>> deltas_pos = make_shared<vector<pair<int,int>>>(0);
  int my_sign = 1; // TODO should be double from range_block
  gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>( 1, make_shared<GammaIntermediateRedux>( ids_pos, deltas_pos, my_sign ) );
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::add_gamma( shared_ptr<Range_Block_Info> block_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::add_gamma" << endl;

  orig_aops_ = block_info->orig_aops(); 
  print_vector(  *(block_info->transformed_idxs()), " *(block_info->transformed_idxs())" ) ; cout << endl;

  //Reordering of the original ids; don't use the original order so we can merge h.c. terms, e.g.,  f1f0x1x0 and x0x1f0f1 . 
  standard_order_ = vector<int>(orig_ids_->size());
  {
    vector<string> tmp_orig_ids = *(block_info->transformed_idxs());
    vector<int> tmp_pos(tmp_orig_ids.size());

    iota( tmp_pos.begin(), tmp_pos.end(), 0);
    sort( tmp_pos.begin(), tmp_pos.end(), [ &tmp_orig_ids ] ( int i1, int i2) { return (bool)( tmp_orig_ids[i1] < tmp_orig_ids[i2]); });  
   
    iota( standard_order_.begin(), standard_order_.end(), 0);
    sort( standard_order_.begin(), standard_order_.end(), [&tmp_pos ] ( int i1, int i2) { return (bool)( tmp_pos[i1] < tmp_pos[i2]); });  

    block_to_std_order_ = tmp_pos; 
  }

  standardized_full_ids_ = vector<string>(orig_ids_->size());
  vector<int>::iterator so_it = standard_order_.begin() ;
  int ii = 0 ;
  for ( vector<string>::const_iterator oi_it = orig_ids_->begin() ; oi_it != orig_ids_->end();  ++oi_it, ++so_it, ++ii ) {
    standardized_full_ids_[*so_it] = (*oi_it);
    block_to_std_order_[*so_it] = (ii);
  }

  shared_ptr<vector<string>> id_ranges = make_shared<vector<string>>(*block_info->unique_block());
  standardized_full_id_ranges_ = vector<string>(*id_ranges);
  vector<int>::iterator otso_it = block_to_std_order_.begin() ;
  for ( vector<string>::iterator sfir_it = standardized_full_id_ranges_.begin() ; sfir_it != standardized_full_id_ranges_.end(); ++sfir_it, ++otso_it )
    *sfir_it = (*id_ranges)[*otso_it]; 

  shared_ptr<vector<pair<int,int>>> deltas_pos = make_shared<vector<pair<int,int>>>(0);
  int my_sign = 1; // TODO should be double from range_block
  gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>( 1, make_shared<GammaIntermediateRedux>( make_shared<vector<int>>( standard_order_ ), deltas_pos, my_sign ) );
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);

  shared_ptr<GammaIntermediateRedux> gint = gamma_vec->front();

  cout << "gint_aops = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << orig_aops_->at(pos) << " " ; cout.flush();  }  cout << "] " << endl;
  cout << "gint_rngs = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << standardized_full_id_ranges_[pos] << " " ; cout.flush();  }   cout << "]" << endl;
  cout << "gint_ids  = [ "; cout.flush();   for ( auto pos : *(gint->ids_pos) ) { cout << orig_ids_->at(pos) << " " ; cout.flush();  }   cout << "] " <<  endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::generic_reorderer( string reordering_name, bool first_reordering, bool final_reordering ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::generic_reorderer" << endl; 
  
  int kk = 0;
  bool does_it_contribute = false;
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediateRedux>>>(0);
  
   for ( string bra_name : *Bra_names_ )
     for ( string ket_name : *Ket_names_ )
       does_it_contribute = generic_reorderer_different_sector( reordering_name, bra_name, ket_name, final_reordering);

   return does_it_contribute;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::generic_reorderer_different_sector( string reordering_name, string bra_name,
                                                              string ket_name, bool final_reordering   ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::generic_reorderer_different_sector" << endl;

  if ( final_reordering ) { 
    cout << "final_reordering !!!" << endl; 
  } else { 
    cout << "NOT final_reordering !!!" << endl; 
  }
 
  shared_ptr<map<char,int>> bra_hole_map = target_states_->hole_range_map(bra_name);;
  shared_ptr<map<char,int>> bra_elec_map = target_states_->elec_range_map(bra_name);;
  shared_ptr<map<char,int>> ket_hole_map = target_states_->hole_range_map(ket_name);;
  shared_ptr<map<char,int>> ket_elec_map = target_states_->elec_range_map(ket_name);;

//  cout << "-------------bra_hole_map------------" << endl;  for ( auto elem : *bra_hole_map ) { cout << elem.first << " = " << elem.second << endl;} cout << endl;
//  cout << "-------------bra_elec_map------------" << endl;  for ( auto elem : *bra_elec_map ) { cout << elem.first << " = " << elem.second << endl;} cout << endl;
//  cout << "-------------ket_hole_map------------" << endl;  for ( auto elem : *ket_hole_map ) { cout << elem.first << " = " << elem.second << endl;} cout << endl;
//  cout << "-------------ket_elec_map------------" << endl;  for ( auto elem : *ket_elec_map ) { cout << elem.first << " = " << elem.second << endl;} cout << endl;
 
  if ( reordering_name == "normal order" ) {
    cout << "doing normal order" << endl;
    int kk = 0;
    cout << "kk = "; cout.flush();
    while ( kk != gamma_vec->size()) {
      cout << kk << " "; cout.flush();
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) )
        normal_order(kk);
      kk++;
      cout << endl;
    }
    cout << "finished normal ordering" << endl;
    kk = 0;
    while ( kk != gamma_vec->size()){
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) && ( gamma_vec->at(kk)->ids_pos->size() != 0 ) )
        final_gamma_vec->push_back( gamma_vec->at(kk ));
      kk++;
    }

  } else if ( reordering_name == "anti-normal order" ) {
    cout << "doing anti-normal order" << endl;
    int kk= 0 ;
    while ( kk != gamma_vec->size()){
      cout << "kk = " << kk << endl;
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) 
        anti_normal_order(kk);
      kk++;
    }
    
    kk = 0;
    cout << "gamma_vec->size() = " << gamma_vec->size() << endl;
    while ( kk != gamma_vec->size()){
      print_gamma_intermediate( gamma_vec->at(kk) );
      if ( proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) && ( gamma_vec->at(kk)->ids_pos->size() != 0 ) ) 
        final_gamma_vec->push_back( gamma_vec->at(kk ));
      kk++ ;
    } 
  } else if ( reordering_name == "alternating order" ) {
    cout << "doing alternating ordering" << endl;
    int kk = 0;
    cout << "kk = "; cout.flush();
    while ( kk != gamma_vec->size()) {
      cout << kk << " "; cout.flush();
      alternating_order(kk);
      kk++;
    } 
    cout << endl << " done alternating ordering" << endl; 
    for ( shared_ptr<GammaIntermediateRedux>& gint : *gamma_vec ){
      if ( proj_onto_map( gint, *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ){ 
        final_gamma_vec->push_back( gint );
      }
    }
  }
  cout << "out of loop" << endl;
  gamma_vec = final_gamma_vec;
  bool does_it_contribute = (gamma_vec->size() > 0 );
 

  int kk = 0;
   
  cout << "gamma_vec->size() = " << gamma_vec->size() << endl; 
  cout << "final_gamma_vec->size() = " << final_gamma_vec->size() << endl; 
  if ( final_reordering && does_it_contribute ) { 
    cout << "final_gamma_vec->size() = " << final_gamma_vec->size() <<  endl;
    while ( kk != gamma_vec->size()){
      add_Acontrib_to_map( kk, bra_name, ket_name );
      kk++;
    } 
  }

  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::braket_survival_check_normal_order( shared_ptr<Range_Block_Info> block_info, 
                                                         shared_ptr<CIVecInfo<double>> bra_info, shared_ptr<CIVecInfo<double>> ket_info ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << endl << endl << " GammaGeneratorRedux::braket_survival_check" << endl;

  return;
}  
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::proj_onto_map( shared_ptr<GammaIntermediateRedux> gint, 
                                    map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                    map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::proj_onto_map" << endl;

  shared_ptr<vector<int>> idxs_pos =  gint->ids_pos;

  shared_ptr<vector<string>> id_ranges = make_shared<vector<string>>(idxs_pos->size());
  shared_ptr<vector<string>> ids = make_shared<vector<string>>(idxs_pos->size());
  shared_ptr<vector<bool>>   aops      = make_shared<vector<bool>>(idxs_pos->size());

  { 
  vector<string>::iterator i_it = ids->begin();
  vector<string>::iterator ir_it = id_ranges->begin();
  vector<bool>::iterator a_it = aops->begin();
  for ( vector<int>::iterator ip_it = idxs_pos->begin(); ip_it != idxs_pos->end(); ip_it++, i_it++, ir_it++, a_it++ ){
    *ir_it = standardized_full_id_ranges_[*ip_it];
    *a_it  = orig_aops_->at(*ip_it);
    *i_it  = orig_ids_->at(*ip_it);
  }
  }
  print_gamma_intermediate(gint); 

  vector<string>::reverse_iterator id_it = ids->rbegin();
  vector<bool>::reverse_iterator aops_it = aops->rbegin();
  for ( auto rng = id_ranges->rbegin(); rng != id_ranges->rend(); rng++, aops_it++, id_it++ ) {
    cout << *id_it << " : " << *aops_it << " : " <<  *rng << endl;
    if( !(*aops_it) ){
      auto ket_elec_map_loc = ket_elec_map.find( (*rng)[0] );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        cout << "Z1" << endl;
        return false;
      } else if ( (ket_elec_map_loc->second -= 1 ) == -1  ) {
        cout << "Z2" << endl;
        return false;
      }
      auto ket_hole_map_loc = ket_hole_map.find( (*rng)[0] );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        ket_hole_map.emplace( (*rng)[0], 1 );
        cout << "Z3" << endl;
      } else {
        cout << "Z4" << endl;
        ket_hole_map_loc->second += 1;
      }
    } else {
      auto ket_hole_map_loc = ket_hole_map.find( (*rng)[0] );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        cout << "Z5" << endl;
        return false;
      } else if ( (ket_hole_map_loc->second -= 1 ) == -1  ) {
        cout << "Z6" << endl;
        return false;
      }
      auto ket_elec_map_loc = ket_elec_map.find( (*rng)[0] );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        cout << "Z7" << endl;
        ket_elec_map.emplace( (*rng)[0], 1 );
      } else {
        cout << "Z8" << endl;
        ket_elec_map_loc->second += 1;
      }
    }
  }
 cout << endl <<endl;
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::anti_normal_order( int kk ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::anti_normal_order" << endl;

  shared_ptr<vector<int>> ids_pos  = gamma_vec->at(kk)->ids_pos;
  int num_kill = 0;
  for ( int pos : *ids_pos )
    if (!orig_aops_->at(pos)) 
      num_kill++; 

  num_kill--; 
 
  for (int ii = ids_pos->size()-1 ; ii != -1; ii--){

    if ( ii > num_kill ) {
      if (orig_aops_->at(ids_pos->at(ii)))
        continue;

      while(!orig_aops_->at( ids_pos->at(ii) )){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if (orig_aops_->at( ids_pos->at(jj) )){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }

    } else if (ii <= num_kill) {
      if (!orig_aops_->at(ids_pos->at(ii)))
        continue;

      while(orig_aops_->at(ids_pos->at(ii))){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if(!orig_aops_->at(ids_pos->at(jj)) ){
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
    if (orig_aops_->at(pos)) 
      num_plus++; 

  num_plus--; 
 
  for (int ii = ids_pos->size()-1 ; ii != -1; ii--){

    if ( ii > num_plus ) {
      if (!orig_aops_->at(ids_pos->at(ii)))
        continue;

      while(orig_aops_->at( ids_pos->at(ii) )){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if (!orig_aops_->at( ids_pos->at(jj) )){
            swap( jj, jj+1, kk, gamma_vec);
            break;
          }
        }
      }

    } else if (ii <= num_plus) {
      if (orig_aops_->at(ids_pos->at(ii)))
        continue;

      while(!orig_aops_->at(ids_pos->at(ii))){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
          if(orig_aops_->at(ids_pos->at(jj)) ){
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
  cout << "Z5" << endl;
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

 
  //Must use standard id ordering here, to match up with entries in CTP_map, and avoid duplication
  //TODO should standardize deltas pos when building gamma intermediate so as to avoid repeated transformation
  vector<pair<int,int>> standardized_deltas_pos(deltas_pos->size());
  vector<pair<int,int>>::iterator sdp_it = standardized_deltas_pos.begin(); 
  for ( pair<int,int>& elem : *deltas_pos ){ 
    *sdp_it = make_pair ( block_to_std_order_[elem.first], block_to_std_order_[elem.second] );
    ++sdp_it;
  }
  
  string Aname_alt = get_Aname( standardized_full_ids_, standardized_full_id_ranges_, standardized_deltas_pos );


  string Gname_alt = get_gamma_name( standardized_full_id_ranges_, *orig_aops_, *ids_pos, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo>>>() );

  //TODO do this reordering w.r.t. standardized orders
  vector<int> Aid_order_new = get_Aid_order ( *ids_pos ) ;
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
                                                         orig_aops_, make_shared<vector<string>>(standardized_full_id_ranges_), ids_pos, Gamma_map) );
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::print_gamma_contributions( shared_ptr<vector<shared_ptr<GammaIntermediateRedux>>> final_gamma_vec, string name) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  print_gamma_contributions( final_gamma_vec, name, Bra_names_->at(0), Ket_names_->at(0) );
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::print_gamma_contributions( shared_ptr<vector<shared_ptr<GammaIntermediateRedux>>> final_gamma_vec,
                                                string name,  string bra_name, string ket_name ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "GammaGeneratorRedux::print_gamma_contributions" << endl;
  
  cout << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << "     LIST OF GAMMAS FOLLOWING " << name;  cout << "N = " << final_gamma_vec->size() << endl;
  cout << "-----------------------------------------------------" << endl;

  for ( shared_ptr<GammaIntermediateRedux> gint : *final_gamma_vec ) {
    string Gname_tmp = WickUtils::get_gamma_name( standardized_full_id_ranges_, *(orig_aops_),  *(gint->ids_pos), bra_name, ket_name) ;
    cout <<Gname_tmp <<  "   ("<< gint->my_sign <<","<< gint->my_sign << ")       " ;
    cout << get_Aname( *(orig_ids_), standardized_full_id_ranges_, *(gint->deltas_pos) ) << endl;
  }
  cout << "-----------------------------------------------------" << endl;

  return;
}

///////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::Contract_remaining_indexes( int kk ){
//////////////////////////////////////////////////////////////////////////////////////
cout << " GammaGeneratorRedux::Contract_remaining_indexes" << endl;
  
  shared_ptr<vector<int>>           ids_pos = gamma_vec->at(kk)->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos;
  int my_sign = gamma_vec->at(kk)->my_sign;


  print_vector( *ids_pos , " ids_pos" ); cout.flush();
  print_vector( standardized_full_id_ranges_ , " standardized_full_id_ranges" ); cout.flush();

  //vector of different index ranges, and vectors containing list of index positions
  //associated with each of these ranges
  vector<string> diff_rngs(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> make_ops_pos =  make_shared<vector<shared_ptr<vector<int>>>>(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> kill_ops_pos =  make_shared<vector<shared_ptr<vector<int>>>>(0);

  int start_pos = 0;
  while( start_pos!= ids_pos->size()  ){
    if ( Forbidden_Index(standardized_full_id_ranges_, ids_pos->at(start_pos)) )
      break;
    start_pos++;
  }

  // get the positions of the creation and annihilation operators associated with
  // each different range
  for ( int jj = start_pos;  jj != ids_pos->size() ; jj++){
    int ii = 0;
    string rng = standardized_full_id_ranges_[ids_pos->at(jj)];
    if ( Forbidden_Index(standardized_full_id_ranges_, ids_pos->at(jj) ) ) {
      do  {
        if( jj != start_pos  &&  rng == diff_rngs[ii] ){
          if ( orig_aops_->at(ids_pos->at(jj)) ){
            make_ops_pos->at(ii)->push_back(ids_pos->at(jj));
          } else {
            kill_ops_pos->at(ii)->push_back(ids_pos->at(jj));
          }
          break;
        }

        if (jj == start_pos || ii == diff_rngs.size()-1 ){
          diff_rngs.push_back(rng);
          shared_ptr<vector<int>> init_vec1_ = make_shared<vector<int>>(1,ids_pos->at(jj));
          shared_ptr<vector<int>> init_vec2_ = make_shared<vector<int>>(0);
          if ( orig_aops_->at(ids_pos->at(jj)) ){
            make_ops_pos->push_back(init_vec1_);
            kill_ops_pos->push_back(init_vec2_);
          } else {
            make_ops_pos->push_back(init_vec2_);
            kill_ops_pos->push_back(init_vec1_);
          }
          break;
        }
        ii++;
      } while (ii != diff_rngs.size());
    }
  }

cout << "Z5" << endl;
  // first index  = range.
  // second index = pair vector defining possible way of contracting indexes with that range.
  // third index  = pair defining contraction of indexes with the same range.
  vector<shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>> new_contractions(diff_rngs.size());

  cout << " new_contractions.size() = " << new_contractions.size() << endl;
  for (int ii =0 ; ii != new_contractions.size(); ii++){
    print_vector( *make_ops_pos->at(ii), " make_ops_pos->at(" + to_string(ii)+")" ); cout.flush(); 
    print_vector( *kill_ops_pos->at(ii), " make_ops_pos->at("+to_string(ii)+")" ); cout.flush(); 
    new_contractions[ii] = get_unique_pairs( make_ops_pos->at(ii), kill_ops_pos->at(ii), make_ops_pos->at(ii)->size() );
  }

  cout << "Z6" << endl;
  shared_ptr<vector<int>> forvec = make_shared<vector<int>>(diff_rngs.size(),0) ;
  shared_ptr<vector<int>> min = make_shared<vector<int>>(diff_rngs.size(),0) ;
  shared_ptr<vector<int>> max = make_shared<vector<int>>(diff_rngs.size()) ;

  for (int ii = 0; ii != max->size();  ii++)
    max->at(ii) = make_ops_pos->at(ii)->size()-1;

  cout << "Z7" << endl;
  do {

    shared_ptr<vector<pair<int,int>>> new_deltas_pos_tmp = make_shared<vector<pair<int,int>>>(*deltas_pos);
    for ( int qq = 0; qq != forvec->size(); qq++)
      new_deltas_pos_tmp->insert(new_deltas_pos_tmp->end(), new_contractions[qq]->at(forvec->at(qq))->begin(), new_contractions[qq]->at(forvec->at(qq))->end());

    //TODO check this is OK now that you've updated standardize delta ordering generic
    shared_ptr<vector<pair<int,int>>> new_deltas_pos = standardize_delta_ordering_generic( new_deltas_pos_tmp ) ;
    shared_ptr<vector<int>> new_ids_pos =  get_unc_ids_from_deltas_ids_comparison( ids_pos , new_deltas_pos );
    final_gamma_vec->push_back(make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas_pos, my_sign));

  } while ( fvec_cycle_skipper(forvec, max, min) ) ;

  cout << "Z8" << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace this with something more sophisticated which uses constraint functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::Forbidden_Index( shared_ptr<const vector<string>> full_id_ranges,  int position ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGeneratorRedux::Forbidden_Index" << endl;

  if ( full_id_ranges->at(position)[0] != 'a' && full_id_ranges->at(position)[0] != 'A'){
    return true;
  } else {
    return false;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace this with something more sophisticated which uses constraint functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::Forbidden_Index( const vector<string>& full_id_ranges,  int position ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGeneratorRedux::Forbidden_Index" << endl;

  if ( full_id_ranges[position][0] != 'a' && full_id_ranges[position][0] != 'A'){
    return true;
  } else {
    return false;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Swaps indexes round, flips sign, and if ranges are the same puts new density matrix in the list.
// CAREFUL : always keep creation operator as the left index in the contraction .
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGeneratorRedux::swap( int ii, int jj, int kk, shared_ptr<vector<shared_ptr<GammaIntermediateRedux>>> gamma_vec  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<int>>  ids_pos             = gamma_vec->at(kk)->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos;

  int idx_buff = ids_pos->at(ii);
  ids_pos->at(ii) = ids_pos->at(jj);
  ids_pos->at(jj) = idx_buff;

  if ( (standardized_full_id_ranges_[ids_pos->at(jj)] == standardized_full_id_ranges_[ids_pos->at(ii)]) &&
       (orig_aops_->at(ids_pos->at(ii)) !=  orig_aops_->at(ids_pos->at(jj))) ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>(*deltas_pos);

    pair<int,int> new_delta = orig_aops_->at(ids_pos->at(jj)) ? make_pair(ids_pos->at(jj), ids_pos->at(ii)) : make_pair(ids_pos->at(ii), ids_pos->at(jj));
    new_deltas_tmp->push_back(new_delta);
    shared_ptr<pint_vec> new_deltas = WickUtils::standardize_delta_ordering_generic( *new_deltas_tmp,  *orig_ids_ );
    int new_sign = gamma_vec->at(kk)->my_sign;

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>();
    for( int qq = 0 ; qq !=ids_pos->size() ; qq++)
      if ( (qq !=ii) && (qq!=jj))
        new_ids_pos->push_back(ids_pos->at(qq));

    shared_ptr<GammaIntermediateRedux> new_gamma = make_shared<GammaIntermediateRedux>( new_ids_pos, new_deltas, new_sign );
    gamma_vec->push_back(new_gamma);

  }
  gamma_vec->at(kk)->my_sign *= -1;

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
        if (orig_ids_->at(deltas_pos->at(ii).first) > orig_ids_->at(deltas_pos->at(jj).first) )
          posvec[ii]++;
 
    for (int ii = 0 ; ii != deltas_pos->size(); ii++)
      new_deltas_pos->at(posvec[ii]) = deltas_pos->at(ii);
 
  } else {
    new_deltas_pos = deltas_pos;
 
  }
  return new_deltas_pos;
}
////////////////////////////////////////////////////////////////////////////////////
//Returns false if gamma contains anything which isn't active alpha (a) or beta (A) 
////////////////////////////////////////////////////////////////////////////////////
bool GammaGeneratorRedux::all_active_ranges(shared_ptr<GammaIntermediateRedux> gint) {
////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGeneratorRedux::all_active_ranges" << endl;

 for ( vector<int>::iterator ip_it = gint->ids_pos->begin() ; ip_it != gint->ids_pos->end(); ip_it++ )
   if ( standardized_full_id_ranges_[*ip_it][0] != 'a' && standardized_full_id_ranges_[*ip_it][0] != 'A')
     return false;

 return true;

}
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGeneratorRedux::get_standard_order ( const vector<string>& rngs ) {
///////////////////////////////////////////////////////////////////////////////////////
cout <<"GammaGeneratorRedux::get_standard_order" << endl;

  vector<string> new_rngs(rngs.size());
  vector<string>::iterator new_rngs_it = new_rngs.begin();
  vector<int> new_order = get_standard_idx_order(rngs) ;

  for ( int pos : new_order ) { *new_rngs_it++ = rngs[pos] ; }
  return new_order;
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
vector<int> GammaGeneratorRedux::get_standardized_alt_order ( const vector<string>& rngs ,const vector<bool>& aops ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGeneratorRedux::get_standardized_alt_order " << endl;
   
  print_vector( rngs, "rngs" ) ; cout << endl;
  vector<int> standard_order = get_standard_idx_order(rngs) ;
  print_vector(standard_order , "standard_order" ) ; cout << endl;

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
void GammaGeneratorRedux::set_standardized_alt_order_unranged ( int kk , vector<int>& standard_alt_order) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGenerator::set_standardized_alt_order_unranged" << endl;

  shared_ptr<vector<int>>           ids_pos = gamma_vec->at(kk)->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 

  vector<string> unc_idxs( ids_pos->size() );
  vector<bool> unc_aops( orig_aops_->size() - 2*deltas_pos->size() );

  vector<string>::iterator ui_it = unc_idxs.begin();
  vector<bool>::iterator ua_it = unc_aops.begin();
  for ( vector<int>::iterator ip_it= ids_pos->begin(); ip_it != ids_pos->end(); ip_it++, ui_it++, ua_it++ ) {
    *ui_it = orig_ids_->at(*ip_it);
    *ua_it = orig_aops_->at(*ip_it);
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
    *ir_it++ =  orig_ids_->at(elem); 
  
  vector<int> standard_order_plus(ids_reordered_pos.size()/2);
  vector<int> standard_order_kill(ids_reordered_pos.size()/2);
  vector<int>::iterator standard_order_plus_it = standard_order_plus.begin();
  vector<int>::iterator standard_order_kill_it = standard_order_kill.begin();
  for ( int pos : ids_reordered_pos) {
    if ( orig_aops_->at(pos) ) {
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
     cout << "gint_aops = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << orig_aops_->at(pos) << " " ; cout.flush();  }  cout << "] " << endl;
     cout << "gint_rngs = [ "; cout.flush();  for ( auto pos : *(gint->ids_pos) ) { cout << standardized_full_id_ranges_[pos] << " " ; cout.flush();  }   cout << "]" << endl;
     cout << "gint_ids  = [ "; cout.flush();   for ( auto pos : *(gint->ids_pos) ) { cout << orig_ids_->at(pos) << " " ; cout.flush();  }   cout << "] " <<  endl;
    return;
}
