#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor for use when building Split_Range_Block_Info, note this does not calculate the inverse transformation; that is determined in SRBI_Helper,  and passed 
// directly to Split_Range_Block_Info
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block,
                                    shared_ptr<vector<int>> idxs_trans,  pair<double,double> factors,
                                    pair<double,double> ReIm_factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info ) :
                                    orig_rngs_(orig_block), idxs_trans_(idxs_trans), factors_(factors), ReIm_factors_(ReIm_factors)  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   cout << "Range_Block_Info::Range_Block_Info 1" << endl;

  num_idxs_ = orig_rngs_->size();
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );
  idxs_trans_inverse_ = make_shared<vector<int>>( num_idxs_ );

  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( int ii = 0 ; ii != num_idxs_; ii++, it_it++ ) 
    (*idxs_trans_inverse_)[ *it_it ] = *it_it;
  }

  name_ = op_info->op_state_name_ + "_";
  for ( auto& rng : *orig_rngs_ ) 
    name_ += rng;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  vector<bool>::const_iterator a_it = aops.begin();
  for ( vector<string>::const_iterator or_it = orig_rngs_->begin(); or_it != orig_rngs_->end(); ++or_it, ++a_it ){
    if (*a_it) {
      plus_pnum_ *= range_to_prime( (*or_it)[0] );
    } else {
      kill_pnum_ *= range_to_prime( (*or_it)[0] );
    }
  }
  if ( (plus_pnum_ - kill_pnum_) != 0) {
    ci_sector_transition_ = true;
  } else {
    ci_sector_transition_ = false;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for use only when building from single operators (otherwise we will have problems with the inverse transformation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<Range_Block_Info> unique_block, 
                                    shared_ptr<Transformation> transform,  pair<double,double> factors, pair<double,double> ReIm_factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(transform->idxs_trans(*orig_block)), factors_(factors),
                                    ReIm_factors_(ReIm_factors), transform_(transform) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "Range_Block_Info::Range_Block_Info 2" << endl;

  num_idxs_ = orig_rngs_->size();
  idxs_trans_inverse_ = make_shared<vector<int>>( num_idxs_ );

  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( int ii = 0 ; ii != num_idxs_; ii++, it_it++ ) 
    (*idxs_trans_inverse_)[ *it_it ] = ii;
  }

  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  name_ = op_info->op_state_name_ + "_";
  for ( auto& rng : *(unique_block_->orig_rngs_) ) 
    name_ += rng;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  vector<bool>::const_iterator a_it = aops.begin();
  for ( vector<string>::const_iterator or_it = orig_rngs_->begin(); or_it != orig_rngs_->end(); ++or_it, ++a_it ){
    if (*a_it) {
      plus_pnum_ *= range_to_prime( (*or_it)[0] );
    } else {
      kill_pnum_ *= range_to_prime( (*or_it)[0] );
    }
  }

  if ( (plus_pnum_ - kill_pnum_)  != 0 ) {
    ci_sector_transition_ = true;
  } else {
    ci_sector_transition_ = false;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for use only when building from single operators (otherwise we will have problems with the inverse transformation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<Range_Block_Info> unique_block, shared_ptr<vector<int>> idxs_trans,
                                    pair<double,double> factors, pair<double,double> ReIm_factors, const vector<bool>& aops, shared_ptr<Op_Info>& op_info ) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_( idxs_trans ), factors_(factors),
                                    ReIm_factors_(ReIm_factors) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::Range_Block_Info 3" << endl;

  num_idxs_ = orig_rngs_->size();
  idxs_trans_inverse_ = make_shared<vector<int>>( num_idxs_ );

  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( int ii = 0 ; ii != num_idxs_; ii++, it_it++ ) 
    (*idxs_trans_inverse_)[ *it_it ] = ii;
  }

  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  name_ = op_info->op_state_name_ + "_";
  for ( auto& rng : *(unique_block->orig_rngs_) ) 
    name_ += rng;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  vector<bool>::const_iterator a_it = aops.begin();
  for ( vector<string>::const_iterator or_it = orig_rngs_->begin(); or_it != orig_rngs_->end(); ++or_it, ++a_it ){
    if (*a_it) {
      plus_pnum_ *= range_to_prime( (*or_it)[0] );
    } else {
      kill_pnum_ *= range_to_prime( (*or_it)[0] );
    }
  }

  if ( (plus_pnum_ - kill_pnum_)  != 0 ) {
    ci_sector_transition_ = true;
  } else {
    ci_sector_transition_ = false;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBI_Helper::SRBI_Helper( std::vector<std::shared_ptr<Range_Block_Info>>& range_blocks, std::vector<int>& cml_sizes, 
                          std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<Op_Info> multiop_info,  
                          std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>>& unique_split_ranges ) :
                          range_blocks_(std::make_shared<std::vector<std::shared_ptr<Range_Block_Info>>>( range_blocks )) , factors_(std::make_pair(1.0,0.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "SRBI_Helper::SRBI_Helper" << endl;

  num_idxs_ = 0;
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++) 
    num_idxs_ += (*rb_it)->num_idxs_; 

  string unique_name = ""; 

  idxs_trans_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator it_it = idxs_trans_->begin();

  idxs_trans_inverse_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator iti_it = idxs_trans_inverse_->begin();

  vector<string> unique_rngs(num_idxs_);
  vector<string>::iterator ur_it = unique_rngs.begin();

  vector<string> orig_rngs( num_idxs_ );
  vector<string>::iterator or_it = orig_rngs.begin();

  //DQ : I feel like there should be a nicer way to do this; a lot of iterators, should I add braces to throw the iterators out?
  vector<int>::iterator cs_it = cml_sizes.begin();
  
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++, cs_it++) {
 
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first =  (Re_buff*(*rb_it)->factors_.first   - Im_buff*(*rb_it)->factors_.second);
    factors_.second = (Re_buff*(*rb_it)->factors_.second  + Im_buff*(*rb_it)->factors_.first );
    
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    std::for_each( it_it, it_it+(*rb_it)->num_idxs_, [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    copy( (*rb_it)->idxs_trans_inverse()->begin(), (*rb_it)->idxs_trans_inverse()->end(), iti_it );
    std::for_each( iti_it, iti_it+(*rb_it)->num_idxs_ , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    copy( (*rb_it)->orig_rngs()->begin(), (*rb_it)->orig_rngs()->end(), or_it );

    copy( (*rb_it)->unique_block_->orig_rngs_->begin(), (*rb_it)->unique_block_->orig_rngs_->end(), ur_it + *cs_it );

    unique_name += (*rb_it)->name_;

    it_it += (*rb_it)->num_idxs_;
    iti_it += (*rb_it)->num_idxs_;
    or_it += (*rb_it)->num_idxs_;
  }

  orig_rngs_ = std::make_shared<const std::vector<string>>(orig_rngs);         

  auto  sssr_loc = unique_split_ranges->find( unique_rngs );
  if ( sssr_loc != unique_split_ranges->end() ) {
    unique_block_ = sssr_loc->second; 
  } else  { 
    shared_ptr<vector<int>> no_trans = make_shared<vector<int>>(unique_rngs.size());
    iota( no_trans->begin(), no_trans->end(), 0);
    unique_block_ = std::make_shared<Range_Block_Info>(make_shared<const vector<string>>(unique_rngs), no_trans, factors_, factors_, *aops, multiop_info ) ;
    unique_block_->unique_block_ = unique_block_->shared_from_this();
  } 

      
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Split_Range_Block_Info::is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "Split_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     //if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}
