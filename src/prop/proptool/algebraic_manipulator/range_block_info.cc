#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_BlockX_Info::Range_BlockX_Info( std::shared_ptr<const std::vector<std::string>> orig_rngs, std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                                      std::shared_ptr<const std::vector<bool>> orig_aops, std::shared_ptr<std::vector<int>> rngs_trans,
                                      std::shared_ptr<std::vector<int>> idxs_trans, std::shared_ptr<std::vector<int>> aops_trans, 
                                      std::pair<double,double> factors  ) :
                                      rngs_trans_(rngs_trans), idxs_trans_(idxs_trans), aops_trans_(aops_trans), 
                                      factors_(factors) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_BlockX_Info::Range_BlockX_Info" << endl;
   
  num_idxs_ = orig_rngs->size();

  name_  = WickUtils::get_ctp_name( *orig_idxs, *orig_rngs );

  vector<bool> trans_aops_(orig_aops->size());
  vector<int>::iterator at_it = aops_trans_->begin();
  for ( vector<bool>::iterator ta_it = trans_aops_.begin(); ta_it != trans_aops_.end(); ta_it++, at_it++ )
    *ta_it = (*orig_aops)[*at_it];

  vector<string> trans_rngs_(num_idxs_);
  vector<int>::iterator rt_it = rngs_trans_->begin();
  for ( vector<string>::iterator tr_it = trans_rngs_.begin(); tr_it != trans_rngs_.end(); tr_it++, rt_it++ )
    *tr_it = (*orig_rngs)[*rt_it];

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  std::vector<bool>::const_iterator ta_it = trans_aops_.begin();
  for ( std::vector<std::string>::const_iterator tr_it = trans_rngs_.begin(); tr_it != trans_rngs_.end(); ++tr_it, ++ta_it ){
    if (*ta_it) {
      plus_pnum_ *= range_to_prime( (*tr_it)[0] );
    } else {
      kill_pnum_ *= range_to_prime( (*tr_it)[0] );
    }
  }

  if ( plus_pnum_ - kill_pnum_ ) {
    no_transition_ = false;
  } else { 
    no_transition_ = true;
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<Range_BlockX_Info> 
Range_BlockX_Info::transform( shared_ptr<const vector<string>> orig_rngs, shared_ptr<const vector<string>> orig_idxs, shared_ptr<const vector<bool>> orig_aops,
                                    vector<int>&  op_order, vector<char> op_trans ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Range_BlockX_Info::transform" << endl;

  if ( (op_order.size() > 1)  || (op_trans.size() > 1) ) {
    cout << "This is a single Range_BlockX_Info, and you are applying transformations for Split_Range_BlockX_Info to it! Aborting!!" << endl;
    assert( (op_order.size() > 1)  || (op_trans.size() > 1) ); 
  } 
  
  char trans = op_trans.front(); 

  vector<int> new_aops_trans  = *aops_trans_; 
  vector<int> new_idxs_trans  = *idxs_trans_; 
  vector<int> new_rngs_trans  = *rngs_trans_; 
  
  
  assert( trans > -1 && trans < 127 ) ;
  
  transform_tens_vec( trans, new_aops_trans ); 
  transform_tens_vec( trans, new_idxs_trans ); 
  transform_tens_vec( trans, new_rngs_trans ); 
  
  pair<double,double> new_factors = factors_; 
  if ( trans == 'H' || trans == 'h' )
    new_factors.second *= -1.0 ; 

  return make_shared<Range_BlockX_Info>( orig_rngs, orig_idxs, orig_aops, make_shared<vector<int>>(new_rngs_trans), make_shared<vector<int>>(new_idxs_trans), make_shared<vector<int>>(new_aops_trans), new_factors  );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBIX_Helper::SRBIX_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks ) :
                            rxnge_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "SRBIX_Helper::SRBIX_Helper" << endl;

  int num_idxs_ = 0;
  unique_   = true;

  vector<int> cml_sizes(range_blocks->size());
  vector<int>::iterator cs_it = cml_sizes.begin();
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  range_blocks->begin(); rb_iter != range_blocks->end();  rb_iter++, cs_it++ ){
    *cs_it += num_idxs_;
     num_idxs_  += (*rb_iter)->num_idxs_;
  } 

  vector<int> idxs_trans(num_idxs_);
  vector<int> aops_trans(num_idxs_);
  vector<int> rngs_trans(num_idxs_);
  vector<int>::iterator it_it = idxs_trans.begin();
  vector<int>::iterator at_it = aops_trans.begin();
  vector<int>::iterator rt_it = rngs_trans.begin();

  //DQ : I feel like there should be a nicer way to do this; a lot of iterators, should I add braces to throw the iterators out?
  cs_it = cml_sizes.begin();
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_it =  range_blocks->begin() ; rb_it != range_blocks->end(); rb_it++, cs_it++) {
 
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first = Re_buff*(*rb_it)->Re_factor() + Im_buff*(*rb_it)->Im_factor();
    factors_.second = Re_buff*(*rb_it)->Im_factor() + Im_buff*(*rb_it)->Re_factor();
    
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    std::for_each( it_it, it_it+(*rb_it)->idxs_trans()->size() , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } ); // these are lambdas 
    
    copy( (*rb_it)->aops_trans()->begin(), (*rb_it)->aops_trans()->end(), at_it );
    std::for_each( at_it, at_it+(*rb_it)->idxs_trans()->size() , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    copy( (*rb_it)->rngs_trans()->begin(), (*rb_it)->rngs_trans()->end(), rt_it );
    std::for_each( rt_it, rt_it+(*rb_it)->idxs_trans()->size() , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } );

    it_it += (*rb_it)->num_idxs_;
    at_it += (*rb_it)->num_idxs_;
    rt_it += (*rb_it)->num_idxs_;
  }

  idxs_trans_      = std::make_shared<std::vector<int>>(idxs_trans);         
  rngs_trans_      = std::make_shared<std::vector<int>>(rngs_trans);        
  aops_trans_      = std::make_shared<std::vector<int>>(aops_trans);         
      
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<Range_BlockX_Info> 
SplitX_Range_Block_Info::transform( shared_ptr<const vector<string>> orig_rngs, shared_ptr<const vector<string>> orig_idxs, shared_ptr<const vector<bool>> orig_aops,
                                    vector<int>&  op_order, vector<char> op_trans ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "SplitX_Range_Block_Info::transform" << endl;

  vector<int> cml_sizes(op_order.size());
  vector<int>::iterator cs_it = cml_sizes.begin(); 
  int cml_size = 0 ;
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++, cs_it++ ){
    *cs_it += cml_size;
     cml_size += (*rb_iter)->num_idxs_;
  }

  vector<int> idxs_trans(num_idxs_);
  vector<int> aops_trans(num_idxs_);
  vector<int> rngs_trans(num_idxs_);
  vector<int>::iterator it_it = idxs_trans.begin();
  vector<int>::iterator at_it = aops_trans.begin();
  vector<int>::iterator rt_it = rngs_trans.begin();

  vector<char>::iterator ot_it = op_trans.begin();
   
  // DQ : is it better to pass by reference to a void function, or have a function return a reference? 
  // note, transformations impact are handled differently depending on the blocks, hence shifting the characters
  for( vector<int>::iterator oo_it = op_order.begin(); oo_it != op_order.end(); oo_it++, ot_it++ ) {

    std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_it =  range_blocks_->begin() + *oo_it;
     
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    transform_tens_vec( (*ot_it), it_it, it_it + (*rb_it)->num_idxs_ ); 
    std::for_each( it_it, it_it+(*rb_it)->num_idxs_, [  &cml_sizes, &oo_it ] ( int &pos ) { pos += cml_sizes[*oo_it] ; } );

    copy( (*rb_it)->rngs_trans()->begin(), (*rb_it)->rngs_trans()->end(), rt_it );
    transform_tens_vec( (*ot_it)+1, rt_it, rt_it + (*rb_it)->num_idxs_ ); 
    std::for_each( rt_it, rt_it+(*rb_it)->num_idxs_, [  &cml_sizes, &oo_it ] ( int &pos ) { pos += cml_sizes[*oo_it] ; } );

    copy( (*rb_it)->aops_trans()->begin(), (*rb_it)->aops_trans()->end(), at_it );
    transform_tens_vec( (*ot_it)+2, at_it,  at_it + (*rb_it)->num_idxs_ ); 
    std::for_each( at_it, at_it+(*rb_it)->num_idxs_, [  &cml_sizes, &oo_it ] ( int &pos ) { pos += cml_sizes[*oo_it] ; } );

    it_it += (*rb_it)->num_idxs_;
    rt_it += (*rb_it)->num_idxs_;
    at_it += (*rb_it)->num_idxs_;
  }

  return make_shared<Range_BlockX_Info>( orig_rngs, orig_idxs, orig_aops, make_shared<vector<int>>(rngs_trans), make_shared<vector<int>>(idxs_trans), make_shared<vector<int>>(aops_trans), factors_  );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool SplitX_Range_Block_Info::is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "SplitX_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     //if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}
