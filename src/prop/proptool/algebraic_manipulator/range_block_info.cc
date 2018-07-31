#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

//#define __DEBUG_PROPTOOL_RANGE_BLOCK_INFO
//#define __DEBUG_PROPTOOL_SPLIT_RANGE_BLOCK_INFO
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor for use when building Split_Range_Block_Info, note this does not calculate the inverse transformation; that is determined in SRBI_Helper,  and passed 
// directly to Split_Range_Block_Info
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block,
                                    shared_ptr<vector<int>> idxs_trans,  pair<double,double> factors,
                                    pair<double,double> ReIm_factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info ) :
                                    orig_rngs_(orig_block), idxs_trans_(idxs_trans), factors_(factors), ReIm_factors_(ReIm_factors),
                                    op_info_(op_info) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RANGE_BLOCK_INFO
cout << "Range_Block_Info::Range_Block_Info 1" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  num_idxs_ = orig_rngs_->size();
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  name_ = op_info->op_state_name_ + "_";
  for ( auto elem : *orig_block )
    name_ += elem;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  set_transition_vars(aops);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for use only when building from single operators (otherwise we will have problems with the inverse transformation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<Range_Block_Info> unique_block, 
                                    shared_ptr<Transformation> transform,  pair<double,double> factors, pair<double,double> ReIm_factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(transform->idxs_trans(*orig_block)), factors_(factors),
                                    ReIm_factors_(ReIm_factors), transform_(transform), op_info_(op_info)  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RANGE_BLOCK_INFO
cout << "Range_Block_Info::Range_Block_Info 2" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "op_info_->name() = "; cout.flush(); cout << op_info_->name_ << endl;

  num_idxs_ = orig_rngs_->size();
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  name_ = op_info->op_state_name_ + "_";
  for ( auto elem : *orig_block )
    name_ += elem;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  set_transition_vars(aops);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for use only when building from single operators (otherwise we will have problems with the inverse transformation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<Range_Block_Info> unique_block, shared_ptr<vector<int>> idxs_trans,
                                    pair<double,double> factors, pair<double,double> ReIm_factors, const vector<bool>& aops, shared_ptr<Op_Info>& op_info ) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_( idxs_trans ), factors_(factors),
                                    ReIm_factors_(ReIm_factors), op_info_(op_info)  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RANGE_BLOCK_INFO
cout << "Range_Block_Info::Range_Block_Info 3" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  num_idxs_ = orig_rngs_->size();
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  name_ = op_info->op_state_name_ + "_";
  for ( auto elem : *orig_block )
    name_ += elem;

  full_op_name_ = op_info->op_full_name_;
  op_state_name_ = op_info->op_state_name_;

  set_transition_vars(aops);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Range_Block_Info::set_transition_vars( const vector<bool>& aops ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RANGE_BLOCK_INFO
cout << "Range_Block_Info::set_transition_vars" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Split_Range_Block_Info::Split_Range_Block_Info( shared_ptr<vector<shared_ptr<Range_Block_Info>>> range_blocks, vector<int>& cml_sizes, 
                                                shared_ptr<vector<bool>> aops, shared_ptr<Op_Info> op_info,  
                                                shared_ptr<map<const vector<string>, shared_ptr<Range_Block_Info>>>& unique_split_ranges ) :
                                                Range_Block_Info( make_pair(1.0, 0.0), make_pair(1.0, 0.0), op_info ),
                                                range_blocks_( range_blocks ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_SPLIT_RANGE_BLOCK_INFO
cout << "Split_Range_Block_Info::Split_Range_Block_Info" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  num_idxs_ = 0;
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++) 
    num_idxs_ += (*rb_it)->num_idxs_; 

  idxs_trans_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator it_it = idxs_trans_->begin();

  vector<string> unique_rngs(num_idxs_);
  vector<string>::iterator ur_it = unique_rngs.begin();

  vector<string> orig_rngs( num_idxs_ );
  vector<string>::iterator or_it = orig_rngs.begin();

  vector<int>::iterator cs_it = cml_sizes.begin();

  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++, cs_it++) {
  
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first =  (Re_buff*(*rb_it)->factors_.first   - Im_buff*(*rb_it)->factors_.second);
    factors_.second = (Re_buff*(*rb_it)->factors_.second  + Im_buff*(*rb_it)->factors_.first );
    
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    std::for_each( it_it, it_it+(*rb_it)->num_idxs_, [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );
  
    copy( (*rb_it)->orig_rngs()->begin(), (*rb_it)->orig_rngs()->end(), or_it );

    copy( (*rb_it)->unique_block_->orig_rngs_->begin(), (*rb_it)->unique_block_->orig_rngs_->end(), ur_it + *cs_it );
  
    it_it += (*rb_it)->num_idxs_;
    or_it += (*rb_it)->num_idxs_;
  }
  
  orig_rngs_ = std::make_shared<const std::vector<string>>(orig_rngs);//TODO annoying const requires copying 
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  set_transition_vars(*aops);

  name_ = op_info->op_state_name_ + "_";

  for ( auto  elem : *orig_rngs_ch_ )
    name_ += elem;

  canonical_ = true; 
  { 
  int tmp_pos = 0; 
  for ( vector<int>::iterator it_it = idxs_trans_->begin(); it_it !=idxs_trans_->end(); it_it++ ) {
    if ( *it_it != tmp_pos ){
      canonical_ = false;
      break;
    }
    tmp_pos++;
  }
  }
  
  if( !canonical_ ) { 
  
    auto  sssr_loc = unique_split_ranges->find( unique_rngs );
    if ( sssr_loc != unique_split_ranges->end() ) {
      unique_block_ = sssr_loc->second;
      range_blocks_canonical_ = unique_block_->range_blocks_canonical();

    } else { 
      vector<int> cml_sizes_canonical( range_blocks_->size() );
      vector<int>::iterator csc_it  = cml_sizes_canonical.begin();

      range_blocks_canonical_ = make_shared<vector<shared_ptr<Range_Block_Info>>>( range_blocks_->size() );
      vector<shared_ptr<Range_Block_Info>>::iterator rbc_it = range_blocks_canonical_->begin();

      vector<string> unique_ranges_canonical( num_idxs_ );
      vector<string>::iterator urc_it = unique_ranges_canonical.begin();

      int tmp = 0; 
      for ( vector<int>::iterator op_it = op_info_->op_order_->begin(); op_it != op_info_->op_order_->end(); op_it++, rbc_it++, csc_it++, urc_it++ ) {
        vector<shared_ptr<Range_Block_Info>>::iterator rb_it = range_blocks_->begin()+(*op_it);
        *rbc_it = (*rb_it)->unique_block_;
        *csc_it = tmp; 
        tmp += (*rb_it)->num_idxs_;
        copy( (*rb_it)->unique_block_->orig_rngs_->begin(), (*rb_it)->unique_block_->orig_rngs_->end(), urc_it  );
      }
      
      unique_block_ = make_shared<Split_Range_Block_Info>( range_blocks_canonical_, cml_sizes_canonical, aops, op_info_->op_info_canonical(),  unique_split_ranges );
    }
  } else {
    range_blocks_canonical_ =  range_blocks_;

  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Split_Range_Block_Info::is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_SPLIT_RANGE_BLOCK_INFO
cout << "Split_Range_Block_Info::is_sparse" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     //if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}
