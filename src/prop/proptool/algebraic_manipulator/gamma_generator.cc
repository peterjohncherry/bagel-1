#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO extend for rel case (Bra and Ket can vary,  prev gammas and
// prev_Bra_info should be constructed from idxs_pos, full_idxs_ranges, full_aops and Ket.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaInfo::GammaInfo ( shared_ptr<CIVecInfo<double>> Bra_info, shared_ptr<CIVecInfo<double>> Ket_info,
                       shared_ptr<const vector<bool>> full_aops_vec, shared_ptr<const vector<string>> full_idx_ranges,
                       shared_ptr<vector<int>> idxs_pos  ,
                       shared_ptr<map<string, shared_ptr<GammaInfo>>> Gamma_map ) :
                       order_(idxs_pos->size()), Bra_info_(Bra_info), Ket_info_(Ket_info) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "GammaInfo:GammaInfo" << endl;
  id_ranges_ = make_shared<vector<string>>(idxs_pos->size());
  aops_      = make_shared<vector<bool>>(idxs_pos->size());
  for (int ii = 0 ; ii != idxs_pos->size(); ii++ ){
    id_ranges_->at(ii) = full_idx_ranges->at(idxs_pos->at(ii));
    aops_->at(ii)      = full_aops_vec->at(idxs_pos->at(ii));
  }

  map< string , vector<int> > shift_map ;
  vector<string> diff_ranges = { "act" } ; // obtain from diff ranges function
  if ( diff_ranges.size() > 1 ){
    for ( string rng : diff_ranges){
      vector<int> shift_vec = vector<int>(idxs_pos->size());
      vector<string>::reverse_iterator rngs_it = id_ranges_->rbegin();
      vector<bool>::reverse_iterator   aops_it = aops_->rbegin();
      int shift = 0;
      for ( vector<int>::reverse_iterator shift_it = shift_vec.rbegin(); shift_it != shift_vec.rend(); shift_it++){
        if ( *rngs_it++ == rng ){
          if ( ( *aops_it++ ) ) {
            shift++;
          } else {
            shift--;
          }
        }
        *shift_it = shift;
        rngs_it++;
      }
      shift_map.emplace(rng, shift_vec);
    }
  }

  sigma_id_ranges_ = make_shared<vector<string>>(id_ranges_->size()+1);
  sigma_id_ranges_->at(0) = Bra_info_->name();
  for ( int  ii = 1 ; ii != sigma_id_ranges_->size() ; ii++ )
    sigma_id_ranges_->at(ii) = id_ranges_->at(ii-1);

  name_ = get_gamma_name( full_idx_ranges, full_aops_vec, idxs_pos, Bra_info_->name(), Ket_info_->name() ); 
  sigma_name_ = "S_"+name_;

  if ( (idxs_pos->size() > 2 ) && ( idxs_pos->size() % 2 == 0 ) ) {

    prev_gammas_ = vector<string>(idxs_pos->size()/2);
    prev_sigmas_ = vector<string>(idxs_pos->size()/2);

    for (int ii = 2 ; ii != idxs_pos->size() ; ii+=2 ){

      shared_ptr<vector<int>> prev_gamma_idxs_pos = make_shared<vector<int>> ( idxs_pos->begin()+ii, idxs_pos->end());
      prev_gammas_[(ii/2)-1] = get_gamma_name( full_idx_ranges, full_aops_vec, prev_gamma_idxs_pos, Bra_info_->name(), Ket_info_->name());
      prev_sigmas_[(ii/2)-1] = "S_"+prev_gammas_[(ii/2)-1];

      if ( Gamma_map->find( prev_gammas_[(ii/2)-1] )  == Gamma_map->end() ){ //TODO fix Bra_info for rel case
        shared_ptr<GammaInfo> prev_gamma_ = make_shared<GammaInfo>( Bra_info_, Ket_info_, full_aops_vec, full_idx_ranges, prev_gamma_idxs_pos, Gamma_map );
        Gamma_map->emplace( prev_gammas_[(ii/2)-1], prev_gamma_ );
      }

      if(  ii == 2 )
        prev_Bra_info_ = Gamma_map->at( prev_gammas_[0] )->Bra_info() ;

    }
    print_vector( prev_gammas_ , "prev gammas of " + name_ ); cout << endl;
  }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaGenerator::GammaGenerator( shared_ptr<StatesInfo<double>> target_states, int Bra_num, int Ket_num,
                                shared_ptr<const vector<string>> orig_ids, shared_ptr<const vector<bool>> orig_aops,
                                shared_ptr<map<string, shared_ptr<GammaInfo>>> Gamma_map_in,
                                shared_ptr<map<string, shared_ptr<map<string, shared_ptr<AContribInfo> >>>> G_to_A_map_in,
                                double bk_factor_in                                                            ):
                                target_states_(target_states), Bra_num_(Bra_num), Ket_num_(Ket_num),
                                orig_ids_(orig_ids), orig_aops_(orig_aops),
                                G_to_A_map(G_to_A_map_in), Gamma_map(Gamma_map_in), bk_factor(bk_factor_in),
                                projected_bra_(false), projected_ket_(false)  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //needed to keep ordering of contractions consistent, should find better way of defining opname
  auto opname = orig_ids_->at(0)[0]; 

  proj_op_name_ = 'X'; //TODO make it so the projector_op name is taken from expression
 
  op_order = make_shared<map< char, int>>();
  orb_exc_deriv_ = false;
  int ii =0;
  op_order->emplace (opname, ii);

  free_aops_ = make_shared<vector<bool>>();
  free_ids_ = make_shared<vector<string>>();

  int orig_pos = 0;
  int free_pos = 0;
  vector<bool>::const_iterator oa_it = orig_aops_->begin();
  for ( vector<string>::const_iterator oi_it = orig_ids_->begin(); oi_it != orig_ids_->end() ; oi_it++, oa_it++, orig_pos++ ) {

    if ( opname != (*oi_it)[0] ){
      opname = (*oi_it)[0];
      op_order->emplace(opname, ++ii);
    }

    if ( (*oi_it)[0] == proj_op_name_ ) {
      orb_exc_deriv_ = true;
    } else { 
      free_aops_->push_back(*oa_it);
      free_ids_->push_back(*oi_it);
      orig_to_free_pos_.emplace( orig_pos, free_pos++); 
    }     

  }

  ii=0; //TODO This should not be needed here (I think); normal idx order when initialized
        //     Have shuffled idx order specified in range_block_info in add_gamma
  idx_order = make_shared<map< string, int>>();
  for ( string elem : *orig_ids_ )
    idx_order->emplace(elem, ii);

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::add_gamma( shared_ptr<range_block_info> block_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::add_gamma" << endl;

  Bra_names_ = target_states_->civec_names( Bra_num_ );
  Ket_names_ = target_states_->civec_names( Ket_num_ );

  shared_ptr<vector<pair<int,int>>> deltas_pos_in = make_shared<vector<pair<int,int>>>(0);
  int my_sign_in = 1;
  shared_ptr<vector<int>> ids_pos_init =  make_shared<vector<int>>(block_info->orig_idxs()->size());
  iota( ids_pos_init->begin() , ids_pos_init->end(), 0 );
  shared_ptr<vector<string>> id_ranges_in = make_shared<vector<string>>(*block_info->orig_block());

  gamma_vec = make_shared<vector<shared_ptr<GammaIntermediate>>>(1, make_shared<GammaIntermediate>(id_ranges_in, ids_pos_init, deltas_pos_in, my_sign_in));
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediate>>>(0);

  cout << " Leaving GammaGenerator::add_gamma" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::generic_reorderer( string reordering_name, bool first_reordering, bool final_reordering ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::generic_reorderer" << endl; 
  
  int kk = 0;
  bool does_it_contribute = false;
  final_gamma_vec = make_shared<vector<shared_ptr<GammaIntermediate>>>(0);
  
   for ( string bra_name : *Bra_names_ )
     for ( string ket_name : *Ket_names_ )
       does_it_contribute = generic_reorderer_different_sector( reordering_name, bra_name, ket_name, final_reordering);

   return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::generic_reorderer_different_sector( string reordering_name, string bra_name,
                                                         string ket_name, bool final_reordering   ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "GammaGenerator::generic_reorderer_different_sector" << endl;

  shared_ptr<map<char,int>> bra_hole_map = target_states_->hole_range_map(bra_name);;
  shared_ptr<map<char,int>> bra_elec_map = target_states_->elec_range_map(bra_name);;
                                                                                    
  shared_ptr<map<char,int>> ket_hole_map = target_states_->hole_range_map(ket_name);;
  shared_ptr<map<char,int>> ket_elec_map = target_states_->elec_range_map(ket_name);;
  
  if ( reordering_name == "normal order" ) {
    
    int kk = 0;
    while ( kk != gamma_vec->size()) {
      if ( !proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ){
        print_vector( *(gamma_vec->at(kk)->ids_pos), "DIES before normal ordering" ); cout << endl << endl;  
      } else { 
        normal_order(kk);
      }
      kk++;
    }
    kk = 0;
    while ( kk != gamma_vec->size()){
      if ( !proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) { 
        print_vector( *gamma_vec->at(kk)->ids_pos, "DIES after normal order with proj between gamma and proj indexes" ); cout << endl; 
      } else if ( gamma_vec->at(kk)->ids_pos->size() != 0 ) { 
        if ( !all_active_ranges(gamma_vec->at(kk)) ) {  
          Contract_remaining_indexes(kk);
        } else {
          final_gamma_vec->push_back( gamma_vec->at(kk ));
        }
      }
      kk++;
    }

  } else if ( reordering_name == "anti-normal order" ) {
    int kk= 0 ;
    while ( kk != gamma_vec->size()){
      if ( !proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) { 
        print_vector( *gamma_vec->at(kk)->ids_pos, "DIES before anti normal ordering" ); cout << endl; 
      } else { 
        anti_normal_order(kk);
      }
      kk++;
    }
    kk = 0;
    while ( kk != gamma_vec->size()){
      if ( !proj_onto_map( gamma_vec->at(kk), *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ) { 
        print_vector( *gamma_vec->at(kk)->ids_pos, "DIES after anti normal ordering" ); cout << endl;
      } else if ( gamma_vec->at(kk)->ids_pos->size() != 0 ) { 
        if ( !all_active_ranges(gamma_vec->at(kk)) ){ 
          Contract_remaining_indexes(kk);
        } else { 
          final_gamma_vec->push_back( gamma_vec->at(kk ));
        }
      }
      kk++ ;
    } 
  } else if ( reordering_name == "alternating order" ) {

    int kk = 0;
    while ( kk != gamma_vec->size())
      alternating_order(kk++);

    for ( shared_ptr<GammaIntermediate>& gint : *gamma_vec ){
      if ( !proj_onto_map( gint, *bra_hole_map, *bra_elec_map, *ket_hole_map, *ket_elec_map ) ){ 
        print_vector( *(gint->ids_pos), "DIES after alternating" );  cout << endl;
      } else {
        final_gamma_vec->push_back( gint );
      }
    }
  }

  gamma_vec = final_gamma_vec;
  bool does_it_contribute = (gamma_vec->size() > 0 );

  if ( does_it_contribute && final_reordering ) 
    print_gamma_contributions( gamma_vec, reordering_name );

  int kk = 0;
  if ( final_reordering ) { 
    if ( !orb_exc_deriv_ ) {  
      while ( kk != gamma_vec->size()){
        add_Acontrib_to_map( kk, bra_name, ket_name );
        kk++;
      } 
    } else {
      while ( kk != gamma_vec->size()){
        add_Acontrib_to_map_orb_deriv( kk, bra_name, ket_name );
        kk++;
      } 
    }
  }

  return does_it_contribute;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::check_if_same_sector( string bra_name, string ket_name ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::check_if_same_sector" << endl;

  shared_ptr<map<char,int>> ket_hole_map = target_states_->hole_range_map(ket_name); 
  shared_ptr<map<char,int>> ket_elec_map = target_states_->elec_range_map(ket_name);
                           
  shared_ptr<map<char,int>> bra_hole_map = target_states_->hole_range_map(bra_name);
  shared_ptr<map<char,int>> bra_elec_map = target_states_->elec_range_map(bra_name);

  for (auto& elem : *ket_elec_map ){
    auto bhm_loc = bra_hole_map->find(elem.first);
    if ( bhm_loc == bra_hole_map->end() || bhm_loc->second != elem.second )
      return false;
  }

  for (auto& elem : *ket_hole_map ) {
    auto bem_loc = bra_elec_map->find(elem.first);
    if ( bem_loc == bra_elec_map->end() || bem_loc->second != elem.second )
      return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::proj_onto_map( shared_ptr<GammaIntermediate> gint, 
                                    map<char,int> bra_hole_map, map<char,int> bra_elec_map,
                                    map<char,int> ket_hole_map, map<char,int> ket_elec_map  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Grossly inefficient, but totally generic, should write seperate routines for normal and antinormal
//ordering; consecutive operators means can just count.
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::proj_onto_map" << endl;

  shared_ptr<vector<int>> idxs_pos =  gint->ids_pos;

  shared_ptr<vector<string>> id_ranges = make_shared<vector<string>>(idxs_pos->size());
  shared_ptr<vector<bool>>   aops      = make_shared<vector<bool>>(idxs_pos->size());
 
  vector<string>::iterator ir_it = id_ranges->begin();
  vector<bool>::iterator a_it = aops->begin();
  for ( vector<int>::iterator ip_it = idxs_pos->begin(); ip_it != idxs_pos->end(); ip_it++ ){
    *ir_it++ = gint->full_id_ranges->at(*ip_it);
    *a_it++  = orig_aops_->at(*ip_it);
  }

  vector<bool>::reverse_iterator aops_it = aops->rbegin();
  for ( auto rng = id_ranges->rbegin(); rng != id_ranges->rend(); rng++, aops_it++ ) {
    if( !(*aops_it) ){
      auto ket_elec_map_loc = ket_elec_map.find( (*rng)[0] );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        return false;
      } else if ( (ket_elec_map_loc->second -= 1 ) == -1  ) {
        return false;
      }
      auto ket_hole_map_loc = ket_hole_map.find( (*rng)[0] );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        ket_hole_map.emplace( (*rng)[0], 1 );
      } else {
        ket_hole_map_loc->second += 1;
      }
    } else {
      auto ket_hole_map_loc = ket_hole_map.find( (*rng)[0] );
      if ( ket_hole_map_loc == ket_hole_map.end() ) {
        return false;
      } else if ( (ket_hole_map_loc->second -= 1 ) == -1  ) {
        return false;
      }
      auto ket_elec_map_loc = ket_elec_map.find( (*rng)[0] );
      if ( ket_elec_map_loc == ket_elec_map.end() ) {
        ket_elec_map.emplace( (*rng)[0], 1 );
      } else {
        ket_elec_map_loc->second += 1;
      }
    }
  }

//  cout << "ket rng, num_ket_elec, num_bra_elec" << endl;
//  for (auto& elem : ket_elec_map ) { 
//    string rng_name = "";
//    rng_name += elem.first;
//    cout << rng_name << " " << elem.second << " " <<  bra_elec_map.at(elem.first) << endl;
//    if ( bra_elec_map.at(elem.first) != elem.second ) {
//      return false;
//    }
//  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::anti_normal_order( int kk ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::anti_normal_order" << endl;

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
void GammaGenerator::normal_order( int kk ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::normal_order" << endl;

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
void GammaGenerator::alternating_order( int kk ) {  // e.g. +-+-+-+-
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "GammaGenerator::alternating_order" << endl;

  shared_ptr<vector<int>>           ids_pos = gamma_vec->at(kk)->ids_pos;
  shared_ptr<const vector<string>>  full_id_ranges = gamma_vec->at(kk)->full_id_ranges;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos;

  vector<string> unc_ranged_idxs(full_id_ranges->size() - 2*deltas_pos->size());
  vector<string>::iterator unc_ranged_idxs_it = unc_ranged_idxs.begin();
  vector<bool> unc_aops(orig_aops_->size() - 2*deltas_pos->size());

  vector<bool>::iterator unc_aops_it = unc_aops.begin();
  for (int pos : *ids_pos ) {
    *unc_ranged_idxs_it++ = orig_ids_->at(pos)+full_id_ranges->at(pos);
    *unc_aops_it++ = orig_aops_->at(pos);
  }

  vector<int> standard_alt_order = get_standardized_alt_order( unc_ranged_idxs, unc_aops );

  shared_ptr<vector<int>> new_ids_pos = reorder_vector( standard_alt_order, *ids_pos);
  for (int ii = ids_pos->size()-1 ; ii != -1; ii-- ){
    if (ids_pos->at(ii) == new_ids_pos->at(ii))
      continue;

    while( ids_pos->at(ii) != new_ids_pos->at(ii) ){
      for ( int jj = (ii-1); jj != -1 ; jj--) {
        if (ids_pos->at(jj) == new_ids_pos->at(ii)){
          swap( jj, jj+1, kk, gamma_vec);
          break;
        }
      }
    }
  }

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
///////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "void GammaGenerator::add_Acontrib_to_map" << endl;

  shared_ptr<GammaIntermediate> gamma_int = gamma_vec->at(kk);

  double my_sign = bk_factor*gamma_int->my_sign;
  shared_ptr<const vector<string>>  full_id_ranges = gamma_int->full_id_ranges;
  shared_ptr<vector<pair<int,int>>> deltas_pos     = gamma_int->deltas_pos;
  shared_ptr<vector<int>> ids_pos        = gamma_int->ids_pos;

  string Aname_alt = get_Aname( *orig_ids_, *full_id_ranges, *deltas_pos );
  string Gname_alt = get_gamma_name( full_id_ranges, orig_aops_, ids_pos, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo>>>() );

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
                                                         orig_aops_, full_id_ranges, ids_pos, Gamma_map) );
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::add_Acontrib_to_map_orb_deriv( int kk, string bra_name, string ket_name ){  // e.g. ++++----
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGenerator::add_Acontrib_to_map_orb_deriv" << endl;

  shared_ptr<GammaIntermediate> gint = gamma_vec->at(kk);

  cout << "X1" << endl;
 
  double my_sign = bk_factor*gint->my_sign;
  shared_ptr<const vector<string>>  full_id_ranges = gint->full_id_ranges;
  shared_ptr<vector<pair<int,int>>> deltas_pos     = gint->deltas_pos;
  shared_ptr<vector<int>> ids_pos = gint->ids_pos;

  cout << "X2" << endl;
  //TODO this should deal with ordering
  vector<int> aids_pos;          // aid positions (to be contracted )
  vector<int> pids_pos;             // pid positions (never contracted )
  for ( int& pos : *ids_pos ){ 
    if ( orig_ids_->at(pos)[0] == proj_op_name_ ){
      pids_pos.push_back(pos); 
    } else { 
      //aids_pos.push_back(orig_to_free_pos_.at(pos));
      aids_pos.push_back(pos);
    }
  } 
  cout << "X3" << endl;
  shared_ptr<vector<int>> gids_pos_new = make_shared<vector<int>>(vector<int>(aids_pos));
   
  vector<string> pure_exc_contractions;  
  vector<pair<int,int>> deltas_pos_purged;  // Deltas without contractions with exc indexes
  for ( pair<int,int>& delta_pos : *deltas_pos ) {
    if ( orig_ids_->at(delta_pos.first)[0] !=  proj_op_name_ ){
      if ( orig_ids_->at(delta_pos.second)[0] !=  proj_op_name_ ){  
        deltas_pos_purged.push_back(delta_pos);
      } else {  
        //aids_pos.push_back( orig_to_free_pos_.at(delta_pos.first) );
        aids_pos.push_back( delta_pos.first );
        pids_pos.push_back( delta_pos.second );  
      }
    } else if ( orig_ids_->at(delta_pos.second)[0] == proj_op_name_ ) { 
      pure_exc_contractions.push_back(gint->full_id_ranges->at(delta_pos.first)); 
    } else {
      //aids_pos.push_back( orig_to_free_pos_.at(delta_pos.second) ); 
      aids_pos.push_back( delta_pos.second ); 
      pids_pos.push_back( delta_pos.first );
    }
  } 
 
  cout << "X4" << endl;
  vector<int> Aid_order_new = get_Aid_order ( aids_pos ) ;
  vector<int> pid_order_new = get_Aid_order ( pids_pos ) ;

  vector<string> full_aids;
  for ( string id : *orig_ids_ ) 
    if ( id[0] != proj_op_name_ ) 
      full_aids.push_back(id); 
 
  cout << "X4a" << endl;
  string Aname_alt = get_Aname( *orig_ids_, *(gint->full_id_ranges), deltas_pos_purged );
  string Gname_alt = get_gamma_name( gint->full_id_ranges, orig_aops_, gids_pos_new, bra_name, ket_name );

  if ( G_to_A_map->find( Gname_alt ) == G_to_A_map->end() )
    G_to_A_map->emplace( Gname_alt, make_shared<map<string, shared_ptr<AContribInfo>>>() );

  pair<double,double> new_fac = make_pair(my_sign, my_sign); 
  auto AInfo_loc =  G_to_A_map->at( Gname_alt )->find(Aname_alt);
  if ( AInfo_loc == G_to_A_map->at( Gname_alt )->end() ) {
    cout << "X5" << endl;
    auto AInfo = make_shared<AContribInfo_ExcDeriv>( Aname_alt, Aid_order_new, pid_order_new, new_fac );
    G_to_A_map->at( Gname_alt )->emplace(Aname_alt, AInfo);
    cout << "X6" << endl;

  } else {

    cout << "X7a" << endl;
    shared_ptr<AContribInfo> AInfo = AInfo_loc->second;
    cout << "X7aa AInfo->aid_orders().size() = " << AInfo->aid_orders().size() << endl;
    for ( int qq = 0 ; qq != AInfo->aid_orders().size(); qq++ ) {
      cout << "X7aaa qq = " << qq << "  AInfo->aid_orders().size() = " << AInfo->aid_orders().size() << endl;
      cout << "X7b qq = " << qq  << endl;
      if( Aid_order_new == *(AInfo->aid_order(qq)) ){
        cout << "X7bb qq = " << qq  << endl;
        AInfo->remaining_uses_ += 1;
        AInfo->total_uses_ += 1;
        cout << "X7bbb qq = " << qq << " AInfo->pid_orders().size() = " << AInfo->pid_orders().size() << endl;
        bool switcher = false;
        for ( int rr = 0 ; rr != AInfo->aid_pid_orders(qq)->size(); rr++ ) {
          cout << "X7c rr = " << rr  << endl;
          if( pid_order_new == *(AInfo->pid_order(qq, rr )) ){
            cout << "X7d rr = " << rr  << endl;
            AInfo->combine_factors(qq, rr, new_fac ) ;
            cout << "X7dd rr = " << rr  << endl;
            switcher = true;
            break;
          } else if ( rr == AInfo->aid_pid_orders(qq)->size()-1) {
            cout << "X7e rr = " << rr  << endl;
            AInfo->add_pid_order(qq, pid_order_new);
            cout << "X7f rr = " << rr  << endl;
            AInfo->add_factor( qq, new_fac );
            cout << "X7g rr = " << rr  << endl;
            switcher = true;
            break;
          }
        }
       
      //  if ( qq == AInfo->aid_orders().size()-1) {
        if ( switcher ) {
          cout << "X7h qq = " << qq  << endl;
          AInfo->add_aid_order( Aid_order_new );
          cout << "X7hh qq = " << qq << "  AInfo->aid_orders().size() = " << AInfo->aid_orders().size() << endl; cout << "X7i qq = " << qq << endl;
          AInfo->add_pid_order( qq, pid_order_new );
          cout << "X7ii qq = " << qq << "  AInfo->pid_orders().size() = " << AInfo->pid_orders().size() << endl; cout << "X7j qq = " << qq << endl;
          break; 
        } 
      }
    }
    cout << "X8" << endl;
  }

  cout << "X9" << endl;
  Gamma_map->emplace( Gname_alt, make_shared<GammaInfo>( target_states_->civec_info(bra_name), target_states_->civec_info(ket_name),
                                                         orig_aops_, full_id_ranges, gids_pos_new, Gamma_map) );
  cout << "X10" << endl;

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::print_gamma_contributions( shared_ptr<vector<shared_ptr<GammaIntermediate>>> final_gamma_vec, string name) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  print_gamma_contributions( final_gamma_vec, name, Bra_names_->at(0), Ket_names_->at(0) );
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::print_gamma_contributions( shared_ptr<vector<shared_ptr<GammaIntermediate>>> final_gamma_vec,
                                                string name,  string bra_name, string ket_name ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "GammaGenerator::print_gamma_contributions" << endl;
  
  cout << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << "     LIST OF GAMMAS FOLLOWING " << name;  cout << "N = " << final_gamma_vec->size() << endl;
  cout << "-----------------------------------------------------" << endl;

  for ( shared_ptr<GammaIntermediate> gint : *final_gamma_vec ) {
    string Gname_tmp = WickUtils::get_gamma_name( gint->full_id_ranges, orig_aops_,  gint->ids_pos, bra_name, ket_name) ;
    cout <<Gname_tmp <<  "   ("<< gint->my_sign <<","<< gint->my_sign << ")       " ;
    cout << get_Aname( *(orig_ids_), *(gint->full_id_ranges), *(gint->deltas_pos) ) << endl;
  }
  cout << "-----------------------------------------------------" << endl;

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Replace this with something more sophisticated which uses constraint functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::Forbidden_Index( shared_ptr<const vector<string>> full_id_ranges,  int position ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::Forbidden_Index" << endl;

  if ( full_id_ranges->at(position)[0] != 'a' && full_id_ranges->at(position)[0] != 'A'){
    return true;
  } else {
    return false;
  }
}
///////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::Contract_remaining_indexes( int kk ){
//////////////////////////////////////////////////////////////////////////////////////
//cout << " GammaGenerator::Contract_remaining_indexes" << endl;
  
  shared_ptr<vector<int>>           ids_pos = gamma_vec->at(kk)->ids_pos;
  shared_ptr<const vector<string>>  full_id_ranges = gamma_vec->at(kk)->full_id_ranges;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos;
  int my_sign = gamma_vec->at(kk)->my_sign;

  //vector of different index ranges, and vectors containing list of index positions
  //associated with each of these ranges
  vector<string> diff_rngs(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> make_ops_pos =  make_shared<vector<shared_ptr<vector<int>>>>(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> kill_ops_pos =  make_shared<vector<shared_ptr<vector<int>>>>(0);

  int start_pos = 0;
  while( start_pos!= ids_pos->size()  ){
    if ( Forbidden_Index(full_id_ranges, ids_pos->at(start_pos)) )
      break;
    start_pos++;
  }

  // get the positions of the creation and annihilation operators associated with
  // each different range
  for ( int jj = start_pos;  jj != ids_pos->size() ; jj++){
    int ii = 0;
    string rng = full_id_ranges->at(ids_pos->at(jj));
    if ( Forbidden_Index(full_id_ranges, ids_pos->at(jj) ) ) {
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

  // first index  = range.
  // second index = pair vector defining possible way of contracting indexes with that range.
  // third index  = pair defining contraction of indexes with the same range.
  vector<shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>> new_contractions(diff_rngs.size());

  for (int ii =0 ; ii != new_contractions.size(); ii++)
    new_contractions[ii] = get_unique_pairs( make_ops_pos->at(ii), kill_ops_pos->at(ii), make_ops_pos->at(ii)->size() );

  shared_ptr<vector<int>> forvec = make_shared<vector<int>>(diff_rngs.size(),0) ;
  shared_ptr<vector<int>> min = make_shared<vector<int>>(diff_rngs.size(),0) ;
  shared_ptr<vector<int>> max = make_shared<vector<int>>(diff_rngs.size()) ;

  for (int ii = 0; ii != max->size();  ii++)
    max->at(ii) = make_ops_pos->at(ii)->size()-1;

  do {

    shared_ptr<vector<pair<int,int>>> new_deltas_pos_tmp = make_shared<vector<pair<int,int>>>(*deltas_pos);
    for ( int qq = 0; qq != forvec->size(); qq++)
      new_deltas_pos_tmp->insert(new_deltas_pos_tmp->end(), new_contractions[qq]->at(forvec->at(qq))->begin(), new_contractions[qq]->at(forvec->at(qq))->end());

    shared_ptr<vector<pair<int,int>>> new_deltas_pos = Standardize_delta_ordering_generic( new_deltas_pos_tmp ) ;
    shared_ptr<vector<int>> new_ids_pos =  get_unc_ids_from_deltas_ids_comparison( ids_pos , new_deltas_pos );
    final_gamma_vec->push_back(make_shared<GammaIntermediate>(full_id_ranges, new_ids_pos, new_deltas_pos, my_sign));

  } while ( fvec_cycle_skipper(forvec, max, min) ) ;

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Swaps indexes round, flips sign, and if ranges are the same puts new density matrix in the list.
// CAREFUL : always keep creation operator as the left index in the contraction .
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::swap( int ii, int jj, int kk, shared_ptr<vector<shared_ptr<GammaIntermediate>>> gamma_vec  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<const vector<string>> full_id_ranges = gamma_vec->at(kk)->full_id_ranges;
  shared_ptr<vector<int>>  ids_pos             = gamma_vec->at(kk)->ids_pos;
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos;

  int idx_buff = ids_pos->at(ii);
  ids_pos->at(ii) = ids_pos->at(jj);
  ids_pos->at(jj) = idx_buff;

  if ( (full_id_ranges->at(ids_pos->at(jj)) == full_id_ranges->at(ids_pos->at(ii))) &&
       (orig_aops_->at(ids_pos->at(ii)) !=  orig_aops_->at(ids_pos->at(jj))) ){

    shared_ptr<pint_vec> new_deltas_tmp = make_shared<pint_vec>(*deltas_pos);

    pair<int,int> new_delta = orig_aops_->at(ids_pos->at(jj)) ? make_pair(ids_pos->at(jj), ids_pos->at(ii)) : make_pair(ids_pos->at(ii), ids_pos->at(jj));
    new_deltas_tmp->push_back(new_delta);
    shared_ptr<pint_vec> new_deltas = Standardize_delta_ordering_generic( new_deltas_tmp );
    int new_sign = gamma_vec->at(kk)->my_sign;

    shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>();
    for( int qq = 0 ; qq !=ids_pos->size() ; qq++)
      if ( (qq !=ii) && (qq!=jj))
        new_ids_pos->push_back(ids_pos->at(qq));

    shared_ptr<GammaIntermediate> new_gamma = make_shared<GammaIntermediate>(full_id_ranges, new_ids_pos, new_deltas, new_sign );
    gamma_vec->push_back(new_gamma);

  }
  gamma_vec->at(kk)->my_sign *= -1;

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec>  GammaGenerator::Standardize_delta_ordering_generic(shared_ptr<pint_vec> deltas_pos ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::Standardize_delta_ordering_generic" << endl;

 shared_ptr<vector<pair<int,int>>> new_deltas_pos;

 if (deltas_pos->size() > 1 ) {
   new_deltas_pos =  make_shared <vector<pair<int,int>>>( deltas_pos->size());
   vector<int> posvec(deltas_pos->size(),0);
   for (int ii = 0 ; ii != deltas_pos->size() ; ii++)
     for (int jj = 0 ; jj != deltas_pos->size() ; jj++)
       if (deltas_pos->at(ii).first > deltas_pos->at(jj).first )
         posvec[ii]++ ;

   for (int ii = 0 ; ii != deltas_pos->size() ; ii++)
     new_deltas_pos->at(posvec[ii]) = deltas_pos->at(ii);

 } else {
   new_deltas_pos = deltas_pos;
 }
  return new_deltas_pos;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec> GammaGenerator::Standardize_delta_ordering(shared_ptr<pint_vec> deltas_pos ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::Standardize_delta_ordering" << endl;

   bool(GammaGenerator::*orderfunc)(pair<string,string>, pair<string,string>) = &GammaGenerator::ordering;

   auto dtmp = make_shared<pint_vec>(*deltas_pos);

   if ( orig_aops_->at( dtmp->back().second ))
     dtmp->back() = make_pair( deltas_pos->back().second, deltas_pos->back().first );

   //sort(dtmp->begin(), dtmp->end() );
   // note that by construction only the last element of dtmp could be in the wrong place.
   int ii = 0;
   while (ii < dtmp->size()-1) {
     if ( op_order->at(orig_ids_->at(dtmp->back().first)[0])  < op_order->at(orig_ids_->at(dtmp->at(ii).first)[0]) ){

       auto dtmp2 = make_shared<pint_vec>();

       for (int qq= 0 ; qq != ii ; qq++)
          dtmp2->push_back(dtmp->at(qq));

       dtmp2->push_back(dtmp->back());

       for (int qq= ii ; qq !=dtmp->size()-1 ; qq++)
          dtmp2->push_back(dtmp->at(qq));

       dtmp = dtmp2;
       break;
      }
      ii++;
  }
  deltas_pos = dtmp;
  return deltas_pos;
}
////////////////////////////////////////////////////////////////////////////////////
//Returns false if gamma contains anything which isn't active alpha (a) or beta (A) 
////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::all_active_ranges(shared_ptr<GammaIntermediate> gint) {
////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::gamma_survives" << endl;

 for ( vector<int>::iterator ip_it = gint->ids_pos->begin() ; ip_it != gint->ids_pos->end(); ip_it++ )
   if ( gint->full_id_ranges->at(*ip_it)[0] != 'a' && gint->full_id_ranges->at(*ip_it)[0] != 'A')
     return false;

 return true;

}
/////////////////////////////////////////////////////////////////////////////////////
//Most basic constraint of all active
//should instead take function
/////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::gamma_survives(shared_ptr<vector<int>> ids_pos, shared_ptr<const vector<string>> id_ranges) {
/////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::gamma_survives" << endl;

 for (int ii = 0 ; ii != ids_pos->size(); ii++  )
   if ( Forbidden_Index(id_ranges, ids_pos->at(ii) ) )
     return false;

 return true;

}
//////////////////////////////////////////////////////////////////////////////
// This should not be necessary, but keep it for debugging
//////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::RangeCheck(shared_ptr<const vector<string>> full_id_ranges) {
//////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::RangeCheck" << endl;

  vector<string> diff_rngs(1, full_id_ranges->at(0) );
  vector<int> updown(1, ( orig_aops_->at(0) ? 1 : -1)  );

  for ( int jj = 1;  jj !=full_id_ranges->size() ; jj++){
    int ii = 0;

    string rng = full_id_ranges->at(jj);

    do {
      if(rng == diff_rngs[ii]){
        if ( orig_aops_->at(jj)){
          updown[ii]+=1;
        } else {
          updown[ii]-=1;
        }
        break;
      }
      if ( ii == diff_rngs.size()-1){
        diff_rngs.push_back(rng);
        if ( orig_aops_->at(jj)){
          updown.push_back(1);
        } else {
          updown.push_back(-1);

        }
        break;
      }

      ii++;
    } while (true);
  }

  for (int ac : updown )
    if (ac != 0 )
      return false;

  return true;

}

//////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_standard_range_order(const vector<string> &rngs) {
//////////////////////////////////////////////////////////////////////////////

  vector<int> pos(rngs.size());
  iota(pos.begin(), pos.end(), 0);
  sort(pos.begin(), pos.end(), [&rngs](int i1, int i2){
                                   return (bool)( rngs[i1] < rngs[i2] );
});

  return pos;
}

//////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_standard_idx_order(const vector<string>&idxs) {
//////////////////////////////////////////////////////////////////////////////

  vector<int> pos(idxs.size());
  iota(pos.begin(), pos.end(), 0);

  auto op_order_tmp = op_order;
  sort(pos.begin(), pos.end(), [&idxs, &op_order_tmp](int i1, int i2){
                                   return (bool)( op_order_tmp->at(idxs[i1][0]) < op_order_tmp->at(idxs[i2][0]) );
                                 }
                            );

  return pos;
}


///////////////////////////////////////////////////////////////////////////////////////
//Returns ordering vector for ids_pos, e.g., if ids_pos = {6,4,5} , then pos = {1,2,0}
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_position_order(const vector<int> &ids_pos) {
///////////////////////////////////////////////////////////////////////////////////////

  vector<int> pos(ids_pos.size());
  iota(pos.begin(), pos.end(), 0);
  sort(pos.begin(), pos.end(), [&ids_pos](int i1, int i2){return (bool)( ids_pos[i1] < ids_pos[i2] ); });

  return pos;
}

///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_standard_order ( const vector<string>& rngs ) {
///////////////////////////////////////////////////////////////////////////////////////

  vector<string> new_rngs(rngs.size());
  vector<string>::iterator new_rngs_it = new_rngs.begin();
  vector<int> new_order = get_standard_idx_order(rngs) ;

  for ( int pos : new_order ) { *new_rngs_it++ = rngs[pos] ; }
  return new_order;
}
///////////////////////////////////////////////////////////////////////////////////////
// Getting reordering vec to go from Atens uncids to gamma uncids
// Running sort twice to get inverse; seems weird and there's probably a better way...
///////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_Aid_order ( const vector<int>& id_pos ) {
///////////////////////////////////////////////////////////////////////////////////////
//   cout << "GammaGenerator::get_Aid_order " << endl;

  vector<int> new_id_pos(id_pos.size());
  vector<int> tmp_order = get_position_order(id_pos);
  vector<int> new_order = get_position_order(tmp_order) ;

  return new_order;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> GammaGenerator::get_standardized_alt_order ( const vector<string>& rngs ,const vector<bool>& aops ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> standard_order = get_standard_order(rngs) ;
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
