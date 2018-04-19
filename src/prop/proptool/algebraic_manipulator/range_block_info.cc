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
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<const vector<string>> unique_block, 
                                    shared_ptr<vector<int>> idxs_trans,  pair<double,double> factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info ) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(idxs_trans), factors_(factors) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Range_Block_Info::Range_Block_Info 1" << endl;

  num_idxs_ = orig_rngs_->size();
  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );
  idxs_trans_inverse_ = make_shared<vector<int>>( num_idxs_ );

  print_vector( *orig_rngs_ , "orig_block " ) ;  print_vector( *unique_block_ , "unique_block " ) ;   print_vector( *idxs_trans_ , "idxs_trans " ) ;  cout << endl;

  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( int ii = 0 ; ii != num_idxs_; ii++, it_it++ ) 
    (*idxs_trans_inverse_)[ *it_it ] = *it_it;
  }


  { // TODO replace this with normal name_, just keep until ctp fixed
  vector<string> idxs(unique_block->size());
  int pos = 0;
  for (auto& elem : idxs ) 
    elem = op_info->op_name_ + to_string(pos++);
  name_ = get_ctp_name( idxs, *unique_block );
  }

  full_op_name_ = op_info->op_full_name_;

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

  if ( plus_pnum_ - kill_pnum_ ) {
    no_transition_ = false;
  } else {
    no_transition_ = true;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor for use only when building from single operators (otherwise we will have problems with the inverse transformation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( shared_ptr<const vector<string>> orig_block, shared_ptr<const vector<string>> unique_block, 
                                    shared_ptr<Transformation> transform,  pair<double,double> factors,
                                    const vector<bool>& aops, shared_ptr<Op_Info>& op_info) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(transform->idxs_trans(*orig_block)), factors_(factors),
                                    transform_(transform) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Range_Block_Info::Range_Block_Info 2" << endl;

  num_idxs_ = orig_rngs_->size();

  idxs_trans_inverse_ = make_shared<vector<int>>( num_idxs_ );

  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( int ii = 0 ; ii != num_idxs_; ii++, it_it++ ) 
    (*idxs_trans_inverse_)[ *it_it ] = ii;
  }

  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  { // TODO replace this with normal name_, just keep until ctp fixed
  vector<string> idxs(unique_block->size());
  int pos = 0;
  for (auto& elem : idxs ) 
    elem = op_info->op_name_ + to_string(pos++);
  name_ = get_ctp_name( idxs, *unique_block );
  }

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

  if ( plus_pnum_ - kill_pnum_ ) {
    no_transition_ = false;
  } else {
    no_transition_ = true;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBI_Helper::SRBI_Helper( std::vector<std::shared_ptr<Range_Block_Info>>& range_blocks, std::vector<int>& cml_sizes ) :
                          range_blocks_(std::make_shared<std::vector<std::shared_ptr<Range_Block_Info>>>( range_blocks )) , factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "SRBI_Helper::SRBI_Helper 2" << endl;

  num_idxs_ = 0;
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++) 
    num_idxs_ += (*rb_it)->num_idxs_; 

  unique_   = true;

  idxs_trans_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator it_it = idxs_trans_->begin();

  idxs_trans_inverse_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator iti_it = idxs_trans_inverse_->begin();

  vector<string> unique_block(num_idxs_);
  vector<string>::iterator ub_it = unique_block.begin();

  vector<string> orig_rngs( num_idxs_ );
  vector<string>::iterator or_it = orig_rngs.begin();

  //DQ : I feel like there should be a nicer way to do this; a lot of iterators, should I add braces to throw the iterators out?
  vector<int>::iterator cs_it = cml_sizes.begin();
  
  cout << endl << endl;
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks_->begin() ; rb_it != range_blocks_->end(); rb_it++, cs_it++) {
 
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first = Re_buff*(*rb_it)->Re_factor() + Im_buff*(*rb_it)->Im_factor();
    factors_.second = Re_buff*(*rb_it)->Im_factor() + Im_buff*(*rb_it)->Re_factor();

    cout <<  (*rb_it)->full_op_name()  << " : " ; cout.flush(); 
    print_vector( *((*rb_it)->idxs_trans()) , "" ) ; 
    print_vector( *((*rb_it)->orig_rngs()) , "" ) ; 
    print_vector( *((*rb_it)->unique_block_) , "" )  ; cout << endl; 
    
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    std::for_each( it_it, it_it+(*rb_it)->num_idxs_, [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    copy( (*rb_it)->idxs_trans_inverse()->begin(), (*rb_it)->idxs_trans_inverse()->end(), iti_it );
    std::for_each( iti_it, iti_it+(*rb_it)->num_idxs_ , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    copy( (*rb_it)->unique_block_->begin(), (*rb_it)->unique_block_->end(), ub_it +(*cs_it) );

    copy( (*rb_it)->orig_rngs()->begin(), (*rb_it)->orig_rngs()->end(), or_it );

    it_it += (*rb_it)->num_idxs_;
    iti_it += (*rb_it)->num_idxs_;
    or_it += (*rb_it)->num_idxs_;
  }
  
  cout << endl << endl;

  unique_block_ = std::make_shared<const std::vector<string>>(unique_block);         
  orig_rngs_ = std::make_shared<const std::vector<string>>(orig_rngs);         
      
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
