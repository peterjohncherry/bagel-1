#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>
#include <tuple> 
using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_Base::TensOp_Base( std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops ) :
                          name_(name), spinfree_(spinfree), Tsymm_("none"), state_dep_(0), 
                          required_blocks_(std::make_shared<std::set<std::string>>()), sub_tensops_(sub_tensops) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_Base::TensOp_Base (MT constructor) " << endl; 

  sort( sub_tensops_.begin(), sub_tensops_.end(), [](shared_ptr<TensOp_Base> t1, shared_ptr<TensOp_Base> t2){ return (bool)( t1->name() < t2->name() );} );

  cout << "sub_tensops_ = [ " ; cout.flush();  for ( auto& t : sub_tensops_)    cout << t->name() << " ";  cout << "]" << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                  vector<bool>& aops, std::pair<double, double>& factor,
                                  vector<shared_ptr<Transformation>>& symmfuncs,
                                  vector<shared_ptr<Constraint>>& constraints,
                                  string& Tsymm, int state_dep, shared_ptr<map<char, long unsigned int>> range_prime_map ) :
                                  TensOp_Base( name, factor, true ), symmfuncs_(symmfuncs), constraints_(constraints)   {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp" <<   endl;

  num_idxs_ = idxs.size();

  all_ranges_state_specific_ =  make_shared<map<string, shared_ptr<map<const vector<string>, shared_ptr<Range_Block_Info>>>>>();
  CTP_map_ = make_shared< map< string, shared_ptr<CtrTensorPart_Base> >>();

  //TODO switch inputs (and input reading) to use shared_ptrs; needless copying here
  aops_ = make_shared<const vector<bool>>(aops);
  idxs_ = make_shared<const vector<string>>(idxs);
  idx_ranges_ = make_shared<const vector<vector<string>>> (idx_ranges); 
  state_ids_ = make_shared<set<shared_ptr<Op_Info>>>();
  generate_blocks();

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_idx_ranges_map() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << " void TensOp::TensOp<DataType>::generate_idx_ranges_map()" << endl;
 
 idx_ranges_map_ = make_shared< map< char, tuple< pair<double,double>, shared_ptr<vector<vector<string>>> >>>();
 pair<double,double> hconj_fac = make_pair ( 1.0 , -1.0); 
 {
 shared_ptr<vector<vector<string>>> transposed_ranges = make_shared<vector<vector<string>>>(*idx_ranges_);
 reverse( transposed_ranges->begin(), transposed_ranges->end() );

 //TODO should be returned from functions defined from symmetry properties of op
 idx_ranges_map_->emplace( 'h' , tie( hconj_fac, transposed_ranges ) );
 } 
 
 {
 shared_ptr<vector<vector<string>>> spin_flip_ranges = make_shared<vector<vector<string>>>(*idx_ranges_);
 for ( vector<string>& rng_vec : *spin_flip_ranges ) 
   for ( auto& elem : rng_vec ) 
     elem[0] = elem[0] > 'Z' ? toupper(elem[0])  : tolower(elem[0]);

 pair<double,double> time_reversal_fac = make_pair ( 1.0 , -1.0); 
 idx_ranges_map_->emplace( 't' , tie( time_reversal_fac,  spin_flip_ranges ) );
 }
 
 //TODO  should got to function which just returns appropriate idx ranges, not real need for map unless
 //      we have lots of symmetries
 shared_ptr<vector<vector<string>>> idx_ranges_nc = make_shared<vector<vector<string>>>(*idx_ranges_);
 pair<double,double> no_fac = make_pair ( 1.0 , 1.0); 
 idx_ranges_map_->emplace( 'n' , tie( no_fac, idx_ranges_nc) ); 
 idx_ranges_map_->emplace( '0' , tie( no_fac, idx_ranges_nc) ); 
 idx_ranges_map_->emplace( 'i' , tie( hconj_fac, idx_ranges_nc ) );
  
 return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
bool TensOp::TensOp<DataType>::satisfies_constraints( vector<string>& ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout <<  " TensOp::TensOp<DataType>::satisfies_constraints" << endl;
  
  for (auto& cstr : constraints_) 
    if(!(cstr->apply_constraint( ranges )) ) {
      return false;
    }
  
  return true;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and ranges; note this is done by TensOp not TensOp gen as the contractions needed will be 
//dependent on the sparsity associated with the relevant state
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_uncontracted_ctps( std::shared_ptr<Op_Info> op_info ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "TensOp::TensOp<DataType> generate_uncontracted_ctps : "  << name_ << endl;

  shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *idxs_ );
  shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector< pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>>  ReIm_factors = make_shared< vector<pair<int,int>>>(1, factor_); 

  for ( vector<shared_ptr<vector<string>>>::iterator bl_it = block_list_.begin(); bl_it != block_list_.end(); bl_it++ ) {
    shared_ptr<vector<string>> full_ranges = *bl_it;
    shared_ptr<CtrTensorPart<DataType>>  CTP = make_shared< CtrTensorPart<DataType> >( full_idxs, full_ranges, noctrs ); 
    CTP_map_->emplace(CTP->name(), CTP);
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_transformed_ranges( shared_ptr<Op_Info> op_info ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType>::generate_transformed_ranges_new : " <<  name_ <<  endl;

  shared_ptr<map<const vector<string>, shared_ptr<Range_Block_Info>>> all_ranges_ref;

  auto arss_loc = all_ranges_state_specific_->find( op_info->op_state_name_ );
  if (arss_loc == all_ranges_state_specific_->end() ) {
    generate_ranges( op_info );
    all_ranges_ref = all_ranges_state_specific_->at( op_info->op_state_name_ );
  } else {
    all_ranges_ref = arss_loc->second;
  }

  shared_ptr<map<const vector<string>, shared_ptr<Range_Block_Info>>> all_ranges_new =
  make_shared<map<const vector<string>, shared_ptr<Range_Block_Info>>>();

  for (auto& elem : *all_ranges_ref ) {
    shared_ptr<Range_Block_Info>  new_block = get_transformed_range_block( op_info, elem.second );
    all_ranges_new->emplace( *(new_block->orig_rngs() ), new_block);
  } 

  all_ranges_state_specific_->emplace( op_info->op_full_name_, all_ranges_new );

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<Range_Block_Info>  
TensOp::TensOp<DataType>::get_transformed_range_block( std::shared_ptr<Op_Info> op_info, std::shared_ptr<Range_Block_Info>& block ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "TensOp::TensOp<DataType>::get_transformed_range_block : " <<  name_ <<  endl;

  shared_ptr<vector<string>> new_orig_rngs = make_shared<vector<string>>(*(block->unique_block_));
  pair<double,double> new_factors = block->factors();

  pair<double,double> trans_factors;
  char op_transformation = tolower(op_info->transformation());

  shared_ptr<vector<int>> new_idxs_trans = make_shared<vector<int>>(block->num_idxs_); 
  iota( new_idxs_trans->begin(), new_idxs_trans->end(), 0 );

  vector<bool> new_aops = *aops_;
 
  switch ( op_transformation ) { 

    case 'i' : // inverse
      for_each( new_aops.begin(), new_aops.end(), [] ( bool aop ) { return !aop; } ); 
      assert ( hermitian_inverse_ ); // TODO in principal, this should go to an op specific function.
      trans_factors = make_pair(1.0, -1.0);
      break;

    case 'h' : // hconj
      for_each( new_aops.begin(), new_aops.end(), [] ( bool aop ) { return !aop; } ); 
      std::reverse( new_orig_rngs->begin(), new_orig_rngs->end() );
      std::reverse( new_idxs_trans->begin(), new_idxs_trans->end() );
      trans_factors = make_pair(1.0, -1.0);
      break;

    case 't' : // time reversal
      for_each( new_orig_rngs->begin(), new_orig_rngs->end(), [] ( string& rng ) { if ( rng[0] > 'Z' ){ rng[0] -=32;} else { rng[0] += 32;}}); 
      trans_factors = make_pair(-1.0, 1.0);
      break;

    case '0' : // do nothing
      trans_factors = make_pair(1.0, 1.0);
      break;

    case 'n' : // do nothing
      trans_factors = make_pair(1.0, 1.0);
      break;

    default : 
      std::cout << "do not have transformation " << op_transformation << " implemented; please check the braket specification in the input file." << std::endl;
      throw std::logic_error( " Aborting!!" );
  } 

  vector<int> id_trans = *(block->idxs_trans());

  reorder_vector_inplace( id_trans, *new_orig_rngs );
  reorder_vector_inplace( id_trans, *new_idxs_trans );

  shared_ptr<const vector<string>> new_orig_rngs_c = std::const_pointer_cast<const vector<string>>(new_orig_rngs);

  return  make_shared<Range_Block_Info>( new_orig_rngs_c, block->unique_block_, new_idxs_trans, factor_, trans_factors, new_aops, op_info );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_blocks() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "TensOp::TensOp<DataType>::generate_blocks" <<   endl;
  
  vector<int> fvec( num_idxs_, 0 );
  vector<int> mins( num_idxs_, 0 );
  vector<int> maxs( num_idxs_ );
  pair<double,double> fac_new(1.0,1.0);
  int num_blocks = 1;

  for( int ii = 0; ii != idx_ranges_->size(); ii++ ){
    maxs[ii] = ( (*idx_ranges_)[ii].size()-1 );
    num_blocks *= (*idx_ranges_)[ii].size();
  }

  // do all symm test ; loop through ranges , try to transform into another..., if fail,  then add as unique range.
  // generate all possible ranges
  block_list_ = vector<shared_ptr<vector<string>>>(num_blocks);
  vector<shared_ptr<vector<string>>>::iterator bl_it = block_list_.begin();
  do { 

    shared_ptr<vector<string>> new_block = make_shared<vector<string>>(idx_ranges_->size());
    for (auto jj = 0 ; jj != fvec.size() ; jj++ )
      (*new_block)[jj] = (*idx_ranges_)[jj][fvec[jj]];

    if( satisfies_constraints( *new_block ) ){
     *bl_it = new_block; 
      ++bl_it;
    } else {
      --num_blocks; 
    }

  } while( fvec_cycle_skipper( fvec, maxs, mins ) );

  block_list_.resize(num_blocks);

  if ( num_blocks <  0 )
    throw logic_error( " The tensor " +  name_ + " appears to have no valid blocks under the current constraints! Aborting!!" ); 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_ranges( std::shared_ptr<Op_Info> op_info ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "TensOp::TensOp<DataType>::generate_ranges : " << name_ <<   endl;

  all_ranges_tmp_ = make_shared< map< const vector<string>, shared_ptr<Range_Block_Info> > > ();

  vector<shared_ptr<vector<string>>>::iterator bl_it = block_list_.begin();
  for ( bl_it = block_list_.begin(); bl_it != block_list_.end() ; bl_it++ )
    apply_symmetry( *(*bl_it), *aops_, op_info ) ;  

  shared_ptr< map< const vector<string>, shared_ptr<Range_Block_Info> > > all_ranges_new = 
            make_shared< map< const vector<string>, shared_ptr<Range_Block_Info> > > ( *all_ranges_tmp_ );

  all_ranges_state_specific_->emplace( op_info->op_state_name_, all_ranges_new );
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::apply_symmetry( const vector<string>& new_block, const vector<bool>& aops, shared_ptr<Op_Info> op_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   cout << "TensOp::TensOp<DataType>::apply_symmetry : " << op_info->name_ <<  endl; 

  pair<double,double> factor_trans;
  shared_ptr<const vector<string>> new_block_c = apply_direct_range_transformation( new_block, factor_, op_info  );
  for ( auto& func : symmfuncs_ ) {
    auto art_loc = all_ranges_tmp_->find(func->transform( *new_block_c ));
    if( art_loc != all_ranges_tmp_->end() ){
      all_ranges_tmp_->emplace( *new_block_c, make_shared<Range_Block_Info>( new_block_c, art_loc->second->unique_block_, func, factor_, factor_trans, aops , op_info ));
      return;
    }
  }

  shared_ptr<vector<int>> idxs_trans = make_shared<vector<int>>(new_block.size());
  iota( idxs_trans->begin(), idxs_trans->end(), 0 );

  new_block_c =  make_shared<const vector<string>>(*new_block_c);
  factor_trans = make_pair(1.0,1.0);
  all_ranges_tmp_->emplace( *new_block_c, make_shared<Range_Block_Info>( new_block_c, new_block_c, idxs_trans, factor_, factor_trans, aops , op_info ) );  

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<const vector<string>>
TensOp::TensOp<DataType>::apply_direct_range_transformation( const vector<string>& block, pair<double,double>& new_fac, 
                                                             shared_ptr<Op_Info>& op_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "TensOp::TensOp<DataType>::apply_direct_range_transformation : " << op_info->name_ <<  endl;

  // TODO should be in seperate function
  vector<string> block_buff = block ;
  vector<string>::const_iterator b_it = block.begin();
  for ( vector<string>::iterator bb_it = block_buff.begin() ; bb_it != block_buff.end() ; bb_it++ , b_it++ )
    (*bb_it)[0]  = tolower( (*b_it)[0] );

//  WickUtils::pair_fac_mult( op_info->spin_flip_factor_, new_fac );

  return make_shared<const vector<string>>( block_buff );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
TensOp::TensOp<DataType>::transform_aops( const char op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "TensOp::TensOp<DataType>::transform_aops " << name_ << endl;

  vector<bool> aops = *aops_; 
  char bob = tolower(op_trans);
  switch ( bob ) { 
    case 'i' : // inverse
      for_each( aops.begin(), aops.end(), [] ( decltype(aops)::reference aop ) { aop = !aop; } ); 
      break;

    case 'h' : // hconj
      for_each( aops.begin(), aops.end(), [] ( decltype(aops)::reference aop ) { aop = !aop; } ); 
      break;

    case 't' : // time reversal
      break;

    case 'n' : // time reversal
      break;

    default : 
      //throw logic_error( "do not have transformation " + bob +" implemented; please check the braket specification in the input file. TensOp::transform_aops"); 
      assert(false);
      break;
  } 

  return make_shared<vector<bool>>(aops);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
MultiTensOp::MultiTensOp<DataType>::transform_aops( const char trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_aops"  << endl;

  throw logic_error( "do not call this yet ; MultiTensOp::MultiTensOp<DataType>::transform_aops , aborting!!"); 

  shared_ptr<vector<bool>> trans_aops = make_shared<vector<bool>>( num_idxs_, false );
  return trans_aops;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
TensOp::TensOp<DataType>::transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType>::transform_aops"  << endl;

  throw logic_error( "do not call this yet ;  TensOp::TensOp<DataType>::transform_aops  aborting!!"); 
  shared_ptr<vector<bool>> trans_aops =  make_shared<vector<bool>>( num_idxs_, false );
  return trans_aops;

}  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
MultiTensOp::MultiTensOp<DataType>::transform_aops( const vector<int>& op_order , const vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_aops " << endl;

  shared_ptr<vector<bool>> trans_aops = make_shared<vector<bool>>( num_idxs_ );

  vector<bool>::iterator ta_it = trans_aops->begin();
  vector<int>::const_iterator oo_it = op_order.begin(); 

  for ( vector<char>::const_iterator ot_it = op_trans.begin(); ot_it != op_trans.end(); ot_it++, oo_it++  ) {

    shared_ptr<vector<bool>> trans_aops_part = sub_tensops_[*oo_it]->transform_aops( *ot_it );
 
    ta_it = move( trans_aops_part->begin(), trans_aops_part->end(), ta_it);

  } 

  return trans_aops;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
MultiTensOp::MultiTensOp<DataType>::MultiTensOp( std::string name, bool spinfree,
                                                 std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops,
                                                 shared_ptr<map< char , long unsigned int>> range_prime_map ):
                                                 TensOp_Base( name, spinfree, sub_tensops ),
                                                 num_tensors_(sub_tensops_.size()) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::MultiTensOp" << endl;
    
  state_specific_split_ranges_ = make_shared<map<string,shared_ptr<map<const vector<string>, shared_ptr<Split_Range_Block_Info>>>>>(); 

  CTP_map_  = make_shared< map< string, shared_ptr<CtrTensorPart_Base> >>();
  { 
    num_idxs_ = 0;
    factor_ = make_pair(1.0,0.0);
    vector<string> idxs;
    vector<vector<string>> idx_ranges;
    vector<bool> aops;
    vector<int> cmlsizevec(num_tensors_);
    for ( int ii = 0;  ii != sub_tensops_.size() ; ii++ ) {
    
      WickUtils::pair_fac_mult( sub_tensops_[ii]->factor_ , factor_ );
      idx_ranges.insert(  idx_ranges.end(), sub_tensops_[ii]->idx_ranges()->begin(), sub_tensops_[ii]->idx_ranges()->end() );
      idxs.insert( idxs.end(),sub_tensops_[ii]->idxs()->begin(), sub_tensops_[ii]->idxs()->end() );
      aops.insert( aops.end(),sub_tensops_[ii]->aops()->begin(), sub_tensops_[ii]->aops()->end() );
      cmlsizevec[ii]  = num_idxs_;
      num_idxs_ += sub_tensops_[ii]->num_idxs();
    
    } 
    aops_ = make_shared<const vector<bool>>(aops);
    idxs_ = make_shared<const vector<string>>(idxs);
    idx_ranges_ = make_shared<const vector<vector<string>>> (idx_ranges); 
    cmlsizevec_ = make_shared<const vector<int>>( cmlsizevec ); 
  } 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
MultiTensOp::MultiTensOp<DataType>::generate_ranges( shared_ptr<Op_Info> multiop_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;
 
  vector< map< const vector<string>, shared_ptr< Range_Block_Info >>::const_iterator> rng_maps( num_tensors_ );  
  vector< map< const vector<string>, shared_ptr< Range_Block_Info >>::const_iterator> rng_maps_begin( num_tensors_ );  
  split_ranges_ = make_shared<map< const vector<string>, shared_ptr<Split_Range_Block_Info> >>();

  shared_ptr<vector<bool>> trans_aops = transform_aops( *(multiop_info->op_order()),  *(multiop_info->transformations()) );
  if ( num_tensors_ > 1 ) { 
    vector<int> forvec(num_tensors_, 0 ); 
    vector<int> old_forvec(num_tensors_, 0 );
    vector<int> mins(num_tensors_, 0 );  
    vector<int> maxs(num_tensors_ );  
    
    vector<int> trans_cml_sizes(rng_maps.size(),0);
    int cml_size = 0; 
    for ( int ii = 0; ii != sub_tensops_.size(); ii++ ){
      sub_tensops_[multiop_info->op_order(ii)]->generate_transformed_ranges( multiop_info->op_info(ii) );

      maxs[ii] = sub_tensops_[multiop_info->op_order(ii)]->all_ranges_state_specific_->at( multiop_info->op_info(ii)->name_ )->size()-1;
      rng_maps[ii] = sub_tensops_[multiop_info->op_order(ii)]->all_ranges_state_specific_->at( multiop_info->op_info(ii)->name_ )->begin();

      trans_cml_sizes[ii] = (*cmlsizevec_)[multiop_info->op_order(ii)];
    }

    rng_maps_begin = rng_maps;

    do {
      for ( int ii = 0 ; ii != forvec.size() ; ii++ ) {
        if ( old_forvec[ii] != forvec[ii] ) {
          if ( forvec[ii] == 0 ) {
            rng_maps[ii] = rng_maps_begin[ii];
          } else {
            rng_maps[ii]++;
          }
        }
      }
 
      old_forvec = forvec;
      shared_ptr< vector <shared_ptr<Range_Block_Info >>> split_block = make_shared< vector <shared_ptr<Range_Block_Info >>>( num_tensors_ );
      vector<shared_ptr<Range_Block_Info>>::iterator sb_it = split_block->begin();
      for (auto  rm_it = rng_maps.begin(); rm_it != rng_maps.end(); rm_it++, sb_it++ )
        *sb_it = (*rm_it)->second;
       
      //TODO Must obtain from constraint functions 
      shared_ptr<Split_Range_Block_Info> srbi;
      {
        SRBI_Helper helper(*split_block, trans_cml_sizes );
        srbi = make_shared<Split_Range_Block_Info>( *trans_aops, helper, multiop_info );
        split_ranges_->emplace( *(helper.orig_rngs_), srbi ) ;
      }

    } while (fvec_cycle_skipper(forvec, maxs, mins )); 

    state_specific_split_ranges_->emplace( multiop_info->name_, split_ranges_ ); 
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::generate_uncontracted_ctps( shared_ptr<MultiOp_Info> op_info ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_uncontracted_ctps " <<  endl;

  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);

  auto split_ranges = state_specific_split_ranges_->at(op_info->name_); 
  for (auto rng_it = split_ranges->begin(); rng_it != split_ranges->end(); rng_it++) 
    enter_cmtps_into_map(*noctrs, rng_it->first );
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::enter_cmtps_into_map(pint_vec ctr_pos_list, const vector<string>& id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // cout << "MultiTensOp::enter_into_CTP_map" << endl;

  shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart_Base> >> (num_tensors_); 
  vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
  vector<vector<pair<int,int>>> sameT_ctrs_pos( num_tensors_,  pint_vec(0));
  shared_ptr<vector<pair<int,int>>> no_ctrs =  make_shared<vector<pair<int,int>>>(0);

  //seperate contractions into those on the same tensor, and those between different tensors 
  // TODO tidy this up, e.g., use lambda to get cross pos
  for ( pair<int,int> ctr_pos : ctr_pos_list ) {
 
    pair<int,int> ctr1;
    pair<int,int> ctr2;

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.first -= sub_tensops_[ii]->num_idxs()) < 0 ){
         ctr1 = make_pair(ii, ctr_pos.first + sub_tensops_[ii]->num_idxs());
        break;
      }

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.second -= sub_tensops_[ii]->num_idxs()) < 0 ){
         ctr2 = make_pair(ii, ctr_pos.second + sub_tensops_[ii]->num_idxs());
        break;
      } 

    if (ctr1.first == ctr2.first) {
      sameT_ctrs_pos.at(ctr1.first).push_back(make_pair(ctr1.second, ctr2.second));  
    } else {
      diffT_ctrs_pos.push_back(make_pair(ctr1,ctr2));
    }
  } 

  //get_ranges for individual tensors
  for (int ii = 0 ; ii !=num_tensors_; ii++ ){ 
    shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(id_ranges.begin() + (*cmlsizevec_)[ii], id_ranges.begin()+(*cmlsizevec_)[ii]+sub_tensops_[ii]->num_idxs());
    shared_ptr<vector<string>>  TS_idxs = make_shared<vector<string>>( *(sub_tensops_[ii]->idxs()) );
   
    if( sameT_ctrs_pos[ii].size() != 0 ) {
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, make_shared<vector<pair<int,int>>>(sameT_ctrs_pos[ii]) ) ; 
    } else { 
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, no_ctrs ); 
    }
    CTP_map_->emplace(CTP_vec->at(ii)->name(), CTP_vec->at(ii)); 

  }
  if ( CTP_vec->size() < 3  ){ 
    auto new_cmtp = make_shared<CtrMultiTensorPart<DataType>>(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos));
    CTP_map_->emplace(new_cmtp->name(), new_cmtp );
  } else {
    get_cmtp(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos) ); 
  } 

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If the ctp_vec has more than two elements then the multitensor is built up recursively from other multitensors, this is
// to simplify the process of contracting the tensors.  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::get_cmtp( shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>>  ctp_vec,  // ctp : contracted_tensor_parts
                                                   shared_ptr<vector<pair<pair<int,int>, pair<int,int>>>> ccp_vec  ){ // ccp : cross_ctrs_pos
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "MultiTensOp::MultiTensOp<DataType>::get_cmtp " << endl;
  
  auto ctp_vec_buff =  make_shared<vector<shared_ptr<CtrTensorPart_Base>>>(*ctp_vec);  
  auto ccp_vec_buff =  make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(*ccp_vec);
  auto cmtp_orig_ctp_order = make_shared<CtrMultiTensorPart<DataType>>(ctp_vec_buff, ccp_vec_buff);

  vector<pair<pair<int,int>, pair<int,int>>>::iterator ccp_it = ccp_vec->begin();
  int counter = 0 ;
  if ( ccp_vec->size() != 0 ) { 
    do {
      int inside_counter = 0 ;
      int ta_pos = ccp_it->first.first < ccp_it->second.first ? ccp_it->first.first : ccp_it->second.first ;
      int tb_pos = ccp_it->first.first > ccp_it->second.first ? ccp_it->first.first : ccp_it->second.first ;

      auto ctp_vec_tatb = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( vector<shared_ptr<CtrTensorPart_Base>> { ctp_vec->at(ta_pos), ctp_vec->at(tb_pos) } );
      auto ccp_vec_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);
      auto ccp_vec_merged_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);

      for ( auto ccp_ab_it = ccp_vec->begin(); ccp_ab_it != ccp_vec->end();  ccp_ab_it++ ) {
        if (ccp_ab_it->first.first == ta_pos && ccp_ab_it->second.first == tb_pos) {
          ccp_vec_tatb->push_back( make_pair( make_pair( 0, ccp_ab_it->first.second), make_pair(1, ccp_ab_it->second.second) ));
        } else if (ccp_ab_it->first.first == tb_pos && ccp_ab_it->second.first == ta_pos){
          ccp_vec_tatb->push_back( make_pair( make_pair( 0, ccp_ab_it->second.second), make_pair(1, ccp_ab_it->first.second) ));
        } else {
          ccp_vec_merged_tatb->push_back(*ccp_ab_it);
        }
      }

      shared_ptr<CtrMultiTensorPart<DataType>> cmtp_tatb = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec_tatb, ccp_vec_tatb );
      CTP_map_->emplace( cmtp_tatb->name(), cmtp_tatb );

      shift_ccp_and_ctp_vecs( cmtp_tatb, ta_pos, tb_pos, ctp_vec, ccp_vec_merged_tatb );

      ccp_vec = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(*ccp_vec_merged_tatb);
      ccp_it = ccp_vec->begin();

    } while( ccp_vec->size() != 0 ) ;

  } else { 
    while ( ctp_vec->size() > 2 ) { 
      auto ctp_vec_tatb = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( vector<shared_ptr<CtrTensorPart_Base>> { *(ctp_vec->end()-2), *(ctp_vec->end()-1) } );
      auto ccp_vec_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);
      auto cmtp_tatb = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec_tatb, ccp_vec_tatb );
      CTP_map_->emplace( cmtp_tatb->name(), cmtp_tatb );
      ctp_vec->pop_back(); 
      ctp_vec->back() = cmtp_tatb;
    }
  }   

  shared_ptr<CtrMultiTensorPart<DataType>> cmtp_new_ctp_order = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec, ccp_vec );
  std::shared_ptr<std::vector<int>> new_order_to_orig_order =  get_pattern_match_order( cmtp_orig_ctp_order->unc_idxs(), cmtp_new_ctp_order->unc_idxs() ) ; 

  cmtp_orig_ctp_order->use_new_order_compute_list( new_order_to_orig_order, cmtp_new_ctp_order->name() );

  CTP_map_->emplace( cmtp_orig_ctp_order->name(), cmtp_orig_ctp_order ) ;
  CTP_map_->emplace( cmtp_new_ctp_order->name(), cmtp_new_ctp_order );

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Shifts the indexes in the cross_ctrs_vec (necessary as ctps at ta and tb have been merged and put to the end or ctp_vec)
// Note that both ctp vec and ccp vec will be altered in this routine
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
MultiTensOp::MultiTensOp<DataType>::shift_ccp_and_ctp_vecs( shared_ptr<CtrMultiTensorPart<DataType>>& tatb_cmtp,
                                                            int ta, int tb, shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>>& ctp_vec,
                                                            shared_ptr<vector<pair<pair<int,int>, pair<int,int>> >>& ccp_vec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "MultiTensOp::MultiTensOp<DataType>::shift_ccp_ctp_vecs" << endl;

  auto new_ctp_vec = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( ctp_vec->size()-1 );
  auto new_ctp_vec_it = new_ctp_vec->begin();
  map<int,int> shifted_ctp_pos_map;
  int tb_id_shift = ctp_vec->at(ta)->size();

  int new_ctp_pos = 0;
  for (int ii = 0 ; ii != ctp_vec->size() ; ii++ ){
    if ( ii != ta && ii != tb ){
      *new_ctp_vec_it++ = ctp_vec->at(ii);
       shifted_ctp_pos_map.emplace(ii,new_ctp_pos++);
    }
  }
  new_ctp_vec->back() = tatb_cmtp;
  ctp_vec = new_ctp_vec;

  shifted_ctp_pos_map.emplace(ta, new_ctp_vec->size()-1);
  shifted_ctp_pos_map.emplace(tb, new_ctp_vec->size()-1);

  for ( pair<pair<int,int> ,pair<int,int>>& ccp : *ccp_vec ) {
    if (ccp.first.first == tb ) {
      ccp.first.second += tb_id_shift;
    } else if (ccp.second.first == tb ) {
      ccp.second.second += tb_id_shift;
    }
    ccp.first.first = shifted_ctp_pos_map.at(ccp.first.first);
    ccp.second.first = shifted_ctp_pos_map.at(ccp.second.first);
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp::TensOp<double>;
template class TensOp::TensOp<complex<double>>;
template class MultiTensOp::MultiTensOp<double>;
template class MultiTensOp::MultiTensOp<complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
