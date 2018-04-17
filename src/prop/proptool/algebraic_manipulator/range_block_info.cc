#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Range_Block_Info::Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_block, std::shared_ptr<const std::vector<std::string>> unique_block, 
                                    std::shared_ptr<std::vector<int>> idxs_trans,  std::pair<double,double> factors,
                                    const std::vector<bool>& aops ) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(idxs_trans), factors_(factors) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::Range_Block_Info SSSS" << endl;

  num_idxs_ = orig_rngs_->size();

  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  unsigned int block_hash = WickUtils::get_block_hash( *orig_block ); // TODO block hashes should be used in original construction.. 

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  std::vector<bool>::const_iterator a_it = aops.begin();
  for ( std::vector<std::string>::const_iterator or_it = orig_rngs_->begin(); or_it != orig_rngs_->end(); ++or_it, ++a_it ){
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
Range_Block_Info::Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_block, std::shared_ptr<const std::vector<std::string>> unique_block, 
                                    std::shared_ptr<Transformation> transform,  std::pair<double,double> factors,
                                    const std::vector<bool>& aops, std::shared_ptr<Op_Info>& op_info) :
                                    orig_rngs_(orig_block), unique_block_(unique_block), idxs_trans_(transform->idxs_trans(*orig_block)), factors_(factors),
                                    transform_(transform) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::Range_Block_Info SSSS" << endl;

  num_idxs_ = orig_rngs_->size();

  orig_rngs_ch_ = make_shared< vector<char>> ( strvec_to_chrvec ( *orig_rngs_ ) );

  unsigned int block_hash = WickUtils::get_block_hash( *orig_block ); // TODO block hashes should be used in original construction.. 

  { // TODO replace this with normal name_, just keep until ctp fixed
  vector<string> idxs(unique_block->size());
  int pos = 0;
  for (auto& elem : idxs ) 
    elem = op_info->op_name_ + to_string(pos++);
  name_ = get_ctp_name( idxs, *unique_block );
  }

  plus_pnum_ = 1;
  kill_pnum_ = 1;

  std::vector<bool>::const_iterator a_it = aops.begin();
  for ( std::vector<std::string>::const_iterator or_it = orig_rngs_->begin(); or_it != orig_rngs_->end(); ++or_it, ++a_it ){
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
// This is for dealing with time reversal symmetry as applied to the Ket or Bra; we want to keep the relevant MO integrals the same, but the ranges
// of the creation and annihilation operators need to be transformed.
// Note that the aops_rngs corresponds to the ranges on which the creation and annhiliation operators act, not the block.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Range_Block_Info::transform_aops_rngs( std::vector<bool>& aops, std::vector<char>& aops_rngs,  std::pair<double,double>& factors , char transformation ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::transform_aops_rngs " << endl;

  transformation = tolower(transformation); 
 
  switch ( transformation ) { 

    case 'i' : // inverse
      for_each( aops.begin(), aops.end(), [] ( bool aop ) { return !aop; } ); 
      factors.second *= -1.0;
      return;

    case 'h' : // hconj
      for_each( aops.begin(), aops.end(), [] ( bool aop ) { return !aop; } ); 
      std::reverse( aops_rngs.begin(), aops_rngs.end() ) ;
      factors.second *= -1.0;
      return;

    case 't' : // time reversal
      for_each( aops_rngs.begin(), aops_rngs.end(), [] ( char rng ) { if ( rng > 'Z' ){ rng -=32;} else { rng += 32;}}); 
      factors.first *= -1.0;
      return;

    default : 
      string trans_str = "" ; trans_str += transformation;
      std::cout << "do not have transformation " << trans_str << "implemented; please check the braket specification in the input file." << std::endl;
      throw std::logic_error( " Aborting!!" );

  } 

  throw std::logic_error( " transform_aops_rngs in rng block Aborting!!" );
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<Range_Block_Info> 
Range_Block_Info::get_transformed_block( shared_ptr<const vector<string> > trans_block, shared_ptr<const vector<string>> unique_block,
                                         char op_trans, vector<bool>& aops ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Range_Block_Info::get_transformed_block" << endl;

  char trans = tolower(op_trans) ;

  shared_ptr<vector<int>> new_idxs_trans  = make_shared<vector<int>>(*idxs_trans_); 
  transform_tens_vec( trans, *new_idxs_trans ); 

  pair<double,double> new_factors = factors_; 
  if ( trans == 'H' || trans == 'h' )
    new_factors.second *= -1.0 ; 

  return make_shared<Range_Block_Info>( trans_block, unique_block, new_idxs_trans, new_factors, aops );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SRBI_Helper::SRBI_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks ) :
                            rxnge_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "SRBI_Helper::SRBI_Helper" << endl;

  int num_idxs_ = 0;
  unique_   = true;

  vector<int> cml_sizes(range_blocks->size());
  vector<int>::iterator cs_it = cml_sizes.begin();
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_iter =  range_blocks->begin(); rb_iter != range_blocks->end();  rb_iter++, cs_it++ ){
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
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  range_blocks->begin() ; rb_it != range_blocks->end(); rb_it++, cs_it++) {
 
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
SRBI_Helper::SRBI_Helper( std::vector<std::shared_ptr<Range_Block_Info>>& range_blocks, std::vector<int>& cml_sizes ) :
                          rxnge_blocks_(std::make_shared<std::vector<std::shared_ptr<Range_Block_Info>>>( range_blocks )) , factors_(std::make_pair(1.0,1.0)) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "SRBI_Helper::SRBI_Helper" << endl;

  num_idxs_ = 0;
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  rxnge_blocks_->begin() ; rb_it != rxnge_blocks_->end(); rb_it++) 
    num_idxs_ += (*rb_it)->num_idxs_; 

  unique_   = true;

  idxs_trans_ = make_shared<vector<int>>(num_idxs_); 
  vector<int>::iterator it_it = idxs_trans_->begin();

  vector<string> unique_block(num_idxs_);
  vector<string>::iterator ub_it = unique_block.begin();

  vector<string> orig_rngs( num_idxs_ );
  vector<string>::iterator or_it = orig_rngs.begin();

  //DQ : I feel like there should be a nicer way to do this; a lot of iterators, should I add braces to throw the iterators out?
  vector<int>::iterator cs_it = cml_sizes.begin();
  
  for ( std::vector<std::shared_ptr<Range_Block_Info>>::iterator rb_it =  rxnge_blocks_->begin() ; rb_it != rxnge_blocks_->end(); rb_it++, cs_it++) {
 
    double Re_buff = factors_.first;
    double Im_buff = factors_.second;
    factors_.first = Re_buff*(*rb_it)->Re_factor() + Im_buff*(*rb_it)->Im_factor();
    factors_.second = Re_buff*(*rb_it)->Im_factor() + Im_buff*(*rb_it)->Re_factor();
    
    copy( (*rb_it)->idxs_trans()->begin(), (*rb_it)->idxs_trans()->end(), it_it );
    std::for_each( it_it, it_it+(*rb_it)->idxs_trans()->size() , [  &cs_it ] ( int &pos ) { pos += *cs_it ; } ); // these are lambdas 

    copy( (*rb_it)->unique_block_->begin(), (*rb_it)->unique_block_->end(), ub_it );

    copy( (*rb_it)->orig_rngs()->begin(), (*rb_it)->orig_rngs()->end(), or_it );

    it_it += (*rb_it)->num_idxs_;
    ub_it += (*rb_it)->num_idxs_;
    or_it += (*rb_it)->num_idxs_;
  }

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
