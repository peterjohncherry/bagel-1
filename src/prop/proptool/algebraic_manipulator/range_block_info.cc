#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( bool is_unique, bool survives, std::pair<double,double> factors, 
                                    std::shared_ptr<const std::vector<std::string>> orig_block, std::shared_ptr<const std::vector<std::string>> unique_block, 
                                    std::shared_ptr<const std::vector<std::string>> orig_idxs, std::shared_ptr<const std::vector<std::string>> transformed_idxs,
                                    std::shared_ptr<const std::vector<bool>> orig_aops,  std::shared_ptr< std::map < char, long unsigned int>> range_prime_map  ) : 
                                    is_unique_(is_unique), survives_(survives),factors_(factors), orig_block_(orig_block), unique_block_(unique_block),
                                    orig_idxs_(orig_idxs), transformed_idxs_(transformed_idxs), orig_aops_(orig_aops), num_idxs_(orig_block_->size()), 
                                    orig_name_(WickUtils::get_Aname(*orig_idxs, *orig_block)), transformed_name_(WickUtils::get_Aname(*transformed_idxs, *unique_block )),
                                    TensOp_name_(orig_idxs->at(0).substr(0,1)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::Range_Block_Info" << endl;
  plus_pnum_ = 1;   
  kill_pnum_ = 1;  

  vector<unsigned int> kill_pos(orig_aops_->size()/2); 
  vector<unsigned int> plus_pos(orig_aops_->size()/2); 
  vector<unsigned int>::iterator kp_it = kill_pos.begin();
  vector<unsigned int>::iterator pp_it = plus_pos.begin(); 

  int pos = 0;
  std::vector<bool>::const_iterator oa_it = orig_aops_->begin(); 
  for ( std::vector<std::string>::const_iterator ob_it = orig_block_->begin(); ob_it != orig_block_->end(); ++ob_it, ++oa_it, ++pos ){ 
    if (*oa_it) { 
      plus_pnum_ *= range_prime_map->at( (*ob_it)[0] );   
      *pp_it = range_prime_map->at( (*ob_it)[0] );
       ++pp_it;
    } else {  
      kill_pnum_ *= range_prime_map->at( (*ob_it)[0] );   
      *kp_it = range_prime_map->at( (*ob_it)[0] );
       ++kp_it;
    }
  }

  // TODO this generates the contractions, should probably be replaced with arithmetical version
  allowed_contractions_ = vector<bool>(kill_pos.size() * plus_pos.size() ); 
  pp_it = plus_pos.begin(); 
  for ( vector<bool>::iterator ac_it = allowed_contractions_.begin(); ac_it != allowed_contractions_.end(); ++pp_it ){ 
    kp_it = kill_pos.begin();
    for( ; kp_it != kill_pos.end(); ++kp_it, ++ac_it ) 
      *ac_it = ( *kp_it == *pp_it );
  } 
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBI_Helper::SRBI_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks ) :
                            range_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "SRBI_Helper::SRBI_Helper" << endl;

  num_idxs_ = 0;
  unique_   = true;

  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){
    num_idxs_  += (*rb_iter)->num_idxs(); } 

  std::vector<std::string> orig_idxs(num_idxs_);              std::vector<std::string>::iterator oi_it = orig_idxs.begin();  
  std::vector<std::string> orig_block(num_idxs_);             std::vector<std::string>::iterator ob_it = orig_block.begin();

  std::vector<std::string> unique_block(num_idxs_);           std::vector<std::string>::iterator ub_it = unique_block.begin();
  std::vector<std::string> transformed_idxs(num_idxs_);       std::vector<std::string>::iterator ti_it = transformed_idxs.begin();

  std::vector<bool> orig_aops(num_idxs_);       std::vector<bool>::iterator oa_it = orig_aops.begin();

  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){

    if ( unique_ && !(*rb_iter)->is_unique() )
      unique_ = false;
 
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first = Re_buff*(*rb_iter)->Re_factor() + Im_buff*(*rb_iter)->Im_factor();
    factors_.second = Re_buff*(*rb_iter)->Im_factor() + Im_buff*(*rb_iter)->Re_factor();

    copy_n( (*rb_iter)->orig_idxs()->begin(), (*rb_iter)->num_idxs(), oi_it );  
    copy_n( (*rb_iter)->unique_block()->begin(), (*rb_iter)->num_idxs(), ub_it );  

    copy_n( (*rb_iter)->orig_block()->begin(), (*rb_iter)->num_idxs(), ob_it );  
    copy_n( (*rb_iter)->transformed_idxs()->begin(), (*rb_iter)->num_idxs(), ti_it );  

    copy_n( (*rb_iter)->orig_aops()->begin(), (*rb_iter)->num_idxs(), oa_it );  

    oi_it += (*rb_iter)->num_idxs();
    ob_it += (*rb_iter)->num_idxs();
    ub_it += (*rb_iter)->num_idxs();
    ti_it += (*rb_iter)->num_idxs();
    oa_it += (*rb_iter)->num_idxs();
 
  }
 
  survives_ =  WickUtils::RangeCheck( orig_block,  orig_aops );

  orig_idxs_        = std::make_shared<const std::vector<std::string>>(orig_idxs);         
  orig_block_       = std::make_shared<const std::vector<std::string>>(orig_block);        
  unique_block_     = std::make_shared<const std::vector<std::string>>(unique_block);      
  transformed_idxs_ = std::make_shared<const std::vector<std::string>>(transformed_idxs);  
  orig_aops_        = std::make_shared<const std::vector<bool>>(orig_aops);         
      
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Split_Range_Block_Info::is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "Split_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}

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
    transform_tens_vec( (*ot_it+2), it_it, it_it + (*rb_it)->num_idxs_ ); 
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
  //cout << "Split_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     //if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}
