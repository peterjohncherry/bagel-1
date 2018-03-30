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
                                      orig_rngs_(orig_rngs), orig_idxs_(orig_idxs), orig_aops_(orig_aops),
                                      rngs_trans_(rngs_trans), idxs_trans_(idxs_trans), aops_trans_(aops_trans), 
                                      factors_(factors) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Range_BlockX_Info::Range_BlockX_Info" << endl;

   
  num_idxs_ = orig_rngs->size();;
 
  vector<bool> trans_aops_(orig_aops_->size());
  vector<int>::iterator at_it = aops_trans_->begin();
  for ( vector<bool>::iterator ta_it = trans_aops_.begin(); ta_it != trans_aops_.end(); ta_it++, at_it++ )
    *ta_it = (*orig_aops_)[*at_it];

  vector<string> trans_rngs_(orig_rngs_->size());
  vector<int>::iterator rt_it = rngs_trans_->begin();
  for ( vector<string>::iterator tr_it = trans_rngs_.begin(); tr_it != trans_rngs_.end(); tr_it++, rt_it++ )
    *tr_it = (*orig_rngs_)[*rt_it];

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  vector<unsigned int> kill_pos(orig_aops_->size()/2);
  vector<unsigned int> plus_pos(orig_aops_->size()/2);
  vector<unsigned int>::iterator kp_it = kill_pos.begin();
  vector<unsigned int>::iterator pp_it = plus_pos.begin(); 

  std::vector<bool>::const_iterator ta_it = trans_aops_.begin();
  for ( std::vector<std::string>::const_iterator tr_it = trans_rngs_.begin(); tr_it != trans_rngs_.end(); ++tr_it, ++ta_it ){
    if (*ta_it) {
      plus_pnum_ *= range_to_prime( (*tr_it)[0] );
      *pp_it = range_to_prime( (*tr_it)[0] );
      ++pp_it;
    } else {
      kill_pnum_ *= range_to_prime( (*tr_it)[0] );
      *kp_it = range_to_prime( (*tr_it)[0] );
      ++kp_it;
    }
  }

  if ( plus_pnum_ - kill_pnum_ ) {
    no_transition_ = false;
  } else { 
    no_transition_ = true;
  }

  // TODO this generates the list of allowed contractions, should probably be replaced with arithmetical version
  allowed_contractions_ = vector<bool>( kill_pos.size() * plus_pos.size() );
  pp_it = plus_pos.begin();
  for ( vector<bool>::iterator ac_it = allowed_contractions_.begin(); ac_it != allowed_contractions_.end(); ++pp_it ){
    kp_it = kill_pos.begin();
    for( ; kp_it != kill_pos.end(); ++kp_it, ++ac_it )
      *ac_it = ( *kp_it == *pp_it );
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBIX_Helper::SRBIX_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks ) :
                            rxnge_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "SRBIX_Helper::SRBIX_Helper" << endl;

  num_idxs_ = 0;
  unique_   = true;

  int counter = 0;
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  rxnge_blocks_->begin(); rb_iter != rxnge_blocks_->end();  rb_iter++ ){
    cout << "num_idxs_ = " << num_idxs_ <<  "    : " << counter++ <<  endl;
    num_idxs_  += (*rb_iter)->num_idxs();
  } 

  cout << "T1" << endl;
  cout << "num_idxs = " << num_idxs_ << endl;
  std::vector<std::string> orig_idxs(num_idxs_);
  cout << "num_idxs = " << num_idxs_ << endl;
  std::vector<std::string>::iterator oi_it = orig_idxs.begin();
  cout << "num_idxs = " << num_idxs_ << endl;

  cout << "T2" << endl;
  std::vector<std::string> orig_rngs(num_idxs_);
  std::vector<std::string>::iterator or_it = orig_rngs.begin();

  cout << "T3" << endl;
  std::vector<bool> orig_aops(num_idxs_);
  std::vector<bool>::iterator oa_it = orig_aops.begin();

  cout << "T4" << endl;
  std::vector<int> idxs_trans(num_idxs_);
  std::vector<int>::iterator it_it = idxs_trans.begin();

  cout << "T5" << endl;
  std::vector<int> rngs_trans(num_idxs_);
  std::vector<int>::iterator rt_it = rngs_trans.begin();

  cout << "T6" << endl;
  std::vector<int> aops_trans(num_idxs_);
  std::vector<int> ::iterator at_it = aops_trans.begin();
  cout << "T7" << endl;
  
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_it =  rxnge_blocks_->begin(); rb_it != rxnge_blocks_->end();  rb_it++ ){
 
    cout << " counter = " << counter++ <<  endl;
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first = Re_buff*(*rb_it)->Re_factor() + Im_buff*(*rb_it)->Im_factor();
    factors_.second = Re_buff*(*rb_it)->Im_factor() + Im_buff*(*rb_it)->Re_factor();

    copy_n( (*rb_it)->orig_idxs()->begin(), (*rb_it)->num_idxs(), oi_it );
    copy_n( (*rb_it)->orig_rngs()->begin(), (*rb_it)->num_idxs(), or_it );
    copy_n( (*rb_it)->orig_aops()->begin(), (*rb_it)->num_idxs(), oa_it );

    copy_n( (*rb_it)->idxs_trans()->begin(), (*rb_it)->num_idxs(), it_it );
    copy_n( (*rb_it)->rngs_trans()->begin(), (*rb_it)->num_idxs(), rt_it );
    copy_n( (*rb_it)->aops_trans()->begin(), (*rb_it)->num_idxs(), at_it );

    oi_it += (*rb_it)->num_idxs();
    or_it += (*rb_it)->num_idxs();
    oa_it += (*rb_it)->num_idxs();
 
    it_it += (*rb_it)->num_idxs();
    rt_it += (*rb_it)->num_idxs();
    at_it += (*rb_it)->num_idxs();

  }
 
  orig_idxs_      = std::make_shared<const std::vector<std::string>>(orig_idxs);         
  orig_rngs_      = std::make_shared<const std::vector<std::string>>(orig_rngs);        
  orig_aops_      = std::make_shared<const std::vector<bool>>(orig_aops);         

  idxs_trans_      = std::make_shared<std::vector<int>>(idxs_trans);         
  rngs_trans_      = std::make_shared<std::vector<int>>(rngs_trans);        
  aops_trans_      = std::make_shared<std::vector<int>>(aops_trans);         
      
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
SRBIX_Helper::add_trans( std::shared_ptr<SplitX_Range_Block_Info> srbi, vector<int>&  op_order, vector<char> op_trans ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "SRBIX_Helper::add_trans" << endl;

  num_idxs_ = 0;
  unique_   = true;
  factors_  = make_pair(1.0,1.0);

  rxnge_blocks_ = srbi->range_blocks() ; 

  std::vector<std::string> orig_idxs = *(srbi->orig_idxs_) ;

  int counter = 0;
  vector<int> cml_sizes(op_order.size());
  vector<int>::iterator cs_it = cml_sizes.begin();
  for ( std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  rxnge_blocks_->begin(); rb_iter != rxnge_blocks_->end();  rb_iter++, cs_it++ ){
    *cs_it += num_idxs_;
     num_idxs_  += (*rb_iter)->num_idxs();
  } 

  vector<int> idxs_trans(num_idxs_);
  vector<int> aops_trans(num_idxs_);
  vector<int> rngs_trans(num_idxs_);
  vector<int>::iterator it_it = idxs_trans.begin();

  for( vector<int>::iterator oo_it = op_order.begin(); oo_it != op_order.end(); oo_it++  ){

    std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_it =  rxnge_blocks_->begin() + *oo_it;
     
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    
    transform_tens_vec( op_trans[*oo_it], it_it, (*rb_it)->idxs_trans()->end() ); 

    std::for_each( it_it, it_it+(*rb_it)->idxs_trans()->size() , [  &cml_sizes, &oo_it ] ( int &pos ) { pos += cml_sizes[*oo_it] ; } );

    it_it += (*rb_it)->idxs_trans()->size();
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool SplitX_Range_Block_Info::is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "Split_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_BlockX_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}


