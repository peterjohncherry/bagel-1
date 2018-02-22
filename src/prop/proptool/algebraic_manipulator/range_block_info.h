#ifndef __SRC_PROP_PROPTOOL_range_block_info_H
#define __SRC_PROP_PROPTOOL_range_block_info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
#include <src/prop/proptool/proputils.h>

class range_block_info : public std::enable_shared_from_this<range_block_info> {

  // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
  protected :
      const  bool is_unique_;
      const  bool survives_;
      const  std::pair<double,double> factors_; 
      const  std::shared_ptr<const std::vector<std::string>> orig_block_;
      const  std::shared_ptr<const std::vector<std::string>> unique_block_;
      const  std::shared_ptr<const std::vector<std::string>> orig_idxs_;
      const  std::shared_ptr<const std::vector<std::string>> transformed_idxs_;
      const  std::shared_ptr<const std::vector<bool>> orig_aops_; 
      const  int num_idxs_;                      
      const  std::string orig_name_;
      const  std::string transformed_name_;
      const  std::string TensOp_name_;

      std::set<std::vector<int>> sparsity_ ;

  public :

    range_block_info( bool is_unique,
                      bool survives,
                      std::pair<double,double> factors, 
                      std::shared_ptr<const std::vector<std::string>> orig_block,   
                      std::shared_ptr<const std::vector<std::string>> unique_block, 
                      std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                      std::shared_ptr<const std::vector<std::string>> transformed_idxs,
                      std::shared_ptr<const std::vector<bool>> orig_aops ) : 
                      is_unique_(is_unique), survives_(survives),factors_(factors),
                      orig_block_(orig_block), unique_block_(unique_block),
                      orig_idxs_(orig_idxs), transformed_idxs_(transformed_idxs),
                      orig_aops_(orig_aops), num_idxs_(orig_block_->size()), 
                      orig_name_(WickUtils::get_Aname(*orig_idxs, *orig_block)),
                      transformed_name_(WickUtils::get_Aname(*transformed_idxs, *unique_block )),
                      TensOp_name_(orig_idxs->at(0).substr(0,1)) {};

    ~range_block_info(){};
    
    bool is_unique() const { return is_unique_ ; } 
    bool survives() const { return survives_ ; }
    
    std::pair<double,double> factors() const { return factors_; } 
    double Re_factor() const { return factors_.first; } 
    double Im_factor() const { return factors_.second; } 
    
    std::shared_ptr<const std::vector<std::string>> orig_block() const { return orig_block_; }
    std::shared_ptr<const std::vector<std::string>> unique_block() const { return unique_block_; }

    std::shared_ptr<const std::vector<std::string>> orig_idxs()const { return orig_idxs_; }
    std::shared_ptr<const std::vector<std::string>> transformed_idxs() const { return transformed_idxs_; }

    std::shared_ptr<const std::vector<bool>> orig_aops() const { return orig_aops_; }
    
    void add_sparse( std::vector<int>& state_idxs ) { sparsity_.emplace(state_idxs); return; } // defines this range block as being sparse for input states 

    int num_idxs() const { return num_idxs_; } 

    std::string TensOp_name() const { return TensOp_name_; }
    std::string orig_name() const { return orig_name_; }
    std::string transformed_name() const { return transformed_name_; }

    // returns true if this block is sparse for input states
    virtual bool is_sparse( std::vector<int>& state_idxs ) { return ( sparsity_.find(state_idxs) != sparsity_.end() ); } 
    virtual bool is_sparse( const std::vector<int>& state_idxs ) { return ( sparsity_.find(state_idxs) != sparsity_.end() ); }  
    virtual std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks(){
        throw std::logic_error("Not a split_range_block; cannot get range_block_vec! Aborting!!" );
          return std::make_shared<std::vector<std::shared_ptr<range_block_info>>>(1, shared_from_this()); }  


};

//This is an absurd hack....
class srbi_helper { 

  public :
    bool unique_;
    bool survives_;
    std::pair<double,double> factors_; 
    int num_idxs_; 
    std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks_;
 
    std::shared_ptr<const std::vector<std::string>> orig_block_;
    std::shared_ptr<const std::vector<std::string>> unique_block_;

    std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    std::shared_ptr<const std::vector<std::string>> transformed_idxs_;

    std::shared_ptr<const std::vector<bool>> orig_aops_;
                              
    srbi_helper( std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks ) :
                            range_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {
      num_idxs_ = 0;
      unique_   = true;

      for ( std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){
        num_idxs_  += (*rb_iter)->num_idxs(); } 

      std::vector<std::string> orig_idxs(num_idxs_);              std::vector<std::string>::iterator oi_it = orig_idxs.begin();  
      std::vector<std::string> orig_block(num_idxs_);             std::vector<std::string>::iterator ob_it = orig_block.begin();

      std::vector<std::string> unique_block(num_idxs_);           std::vector<std::string>::iterator ub_it = unique_block.begin();
      std::vector<std::string> transformed_idxs(num_idxs_);       std::vector<std::string>::iterator ti_it = transformed_idxs.begin();

      std::vector<bool> orig_aops(num_idxs_);       std::vector<bool>::iterator oa_it = orig_aops.begin();

      for ( std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){

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
 
   ~srbi_helper(){};

};

class split_range_block_info : public  range_block_info { 

  private : 
    std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks_;

  public :

    split_range_block_info( srbi_helper& helper ) :
                      range_block_info( helper.unique_, helper.survives_, helper.factors_, helper.orig_block_,
                                        helper.unique_block_, helper.orig_idxs_, helper.transformed_idxs_, helper.orig_aops_ ),
                                        range_blocks_(helper.range_blocks_) {}  
   ~split_range_block_info(){};

    int num_idxs() { return num_idxs_ ; } 

    std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks(){ return range_blocks_ ;} 

    //TODO can't really need three is_sparse functions...
    bool is_sparse( std::vector<std::vector<int>>& state_idxs ) { 
      std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin();
      for ( std::vector<std::vector<int>>::iterator si_iter = state_idxs.begin(); si_iter != state_idxs.end(); si_iter++ ){
         if ( (*rb_iter)->is_sparse(*si_iter) ) 
           return true;      
         rb_iter++;
      }
      return false; 
    }

    bool is_sparse( const std::vector<std::vector<int>>& state_idxs ) { 
      std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin();
      for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs.begin(); si_iter != state_idxs.end(); si_iter++ ){
         if ( (*rb_iter)->is_sparse(*si_iter) ) 
           return true;      
         rb_iter++;
      }
      return false; 
    }

    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs ) { 
      std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin();
      for ( std::vector<std::vector<int>>::const_iterator si_iter = state_idxs->begin(); si_iter != state_idxs->end(); si_iter++ ){
         if ( (*rb_iter)->is_sparse(*si_iter) ) 
           return true;      
         rb_iter++;
      }
      return false; 
    }

};


#endif
