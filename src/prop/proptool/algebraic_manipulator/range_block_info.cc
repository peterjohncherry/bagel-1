#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>

using namespace std;
using namespace WickUtils;

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
  cout << "range_blcok_info::Range_Block_Info" << endl;
  plus_pnum_ = 1;   
  kill_pnum_ = 1;  
  std::vector<bool>::const_iterator oa_it = orig_aops_->begin(); 
  for ( std::vector<std::string>::const_iterator ob_it = orig_block_->begin(); ob_it != orig_block_->end(); ob_it++, oa_it++ ){ 
    if (*oa_it) { 
      plus_pnum_ *= range_prime_map->at( (*ob_it)[0] );   
    } else {  
      kill_pnum_ *= range_prime_map->at( (*ob_it)[0] );   
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBI_Helper::SRBI_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks ) :
                            range_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "SRBI_Helper::SRBI_Helper" << endl;

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
  cout << "Split_Range_Block_Info::is_sparse" << endl;

  std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks_->begin();
  for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
     if ( (*rb_iter)->is_sparse(*si_iter) ) 
       return true;      
     rb_iter++;
  }
  return false; 
}
