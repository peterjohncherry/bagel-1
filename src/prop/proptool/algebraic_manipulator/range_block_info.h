#ifndef __SRC_PROP_PROPTOOL_range_block_info_H
#define __SRC_PROP_PROPTOOL_range_block_info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
class range_block_info {

  // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
  protected :
    const bool is_unique_;
    const bool survives_;
    const std::pair<double,double> factors_; 
    const std::shared_ptr<const std::vector<std::string>> orig_block_;
    const std::shared_ptr<const std::vector<std::string>> unique_block_;
    const std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    const std::shared_ptr<const std::vector<std::string>> transformed_idxs_;
                      
    std::set<std::vector<int>> sparsity_ ;

  public :

    const int num_idxs_;

    range_block_info( bool is_unique,
                      bool survives,
                      std::pair<double,double> factors, 
                      std::shared_ptr<const std::vector<std::string>> orig_block,   
                      std::shared_ptr<const std::vector<std::string>> unique_block, 
                      std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                      std::shared_ptr<const std::vector<std::string>> transformed_idxs) : 
                      is_unique_(is_unique), survives_(survives),factors_(factors),
                      orig_block_(orig_block), unique_block_(unique_block),
                      orig_idxs_(orig_idxs), transformed_idxs_(transformed_idxs),
                      num_idxs_(orig_block_->size()) {};

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
    
    void add_sparse( std::vector<int>& state_idxs ) { sparsity_.emplace(state_idxs); return; } // defines this range block as being sparse for input states 

    int num_idxs() const { return num_idxs_; } 

    virtual bool is_sparse( std::vector<int>& state_idxs ) { return ( sparsity_.find(state_idxs) != sparsity_.end() ); }  // returns true if this block is sparse for input states


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

                              
    srbi_helper( std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks ) :
                            range_blocks_(range_blocks), factors_(std::make_pair(1.0,1.0)) {

      num_idxs_ = 0;
      for ( std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){
        num_idxs_  += (*rb_iter)->num_idxs(); } 

      std::vector<std::string> orig_idxs(num_idxs_);              std::vector<std::string>::iterator oi_it = orig_idxs.begin();  
      std::vector<std::string> orig_block(num_idxs_);             std::vector<std::string>::iterator ob_it = orig_block.begin();
      std::vector<std::string> unique_block(num_idxs_);           std::vector<std::string>::iterator ub_it = unique_block.begin();
      std::vector<std::string> transformed_idxs(num_idxs_);       std::vector<std::string>::iterator ti_it = transformed_idxs.begin();

      for ( std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin(); rb_iter != range_blocks_->end();  rb_iter++ ){

        if ( unique_ && !(*rb_iter)->is_unique() )
          unique_ = false;
 
        if ( survives_ && !(*rb_iter)->survives() )
          survives_ = false;
 
        double Re_buff = factors_.first;
        double Im_buff = factors_.second;
        factors_.first = Re_buff*(*rb_iter)->Re_factor() + Im_buff*(*rb_iter)->Im_factor();
        factors_.second = Re_buff*(*rb_iter)->Im_factor() + Im_buff*(*rb_iter)->Re_factor();

        copy_n( (*rb_iter)->orig_idxs()->begin(), (*rb_iter)->num_idxs(), oi_it );  
        copy_n( (*rb_iter)->orig_block()->begin(), (*rb_iter)->num_idxs(), ob_it );  
        copy_n( (*rb_iter)->unique_block()->begin(), (*rb_iter)->num_idxs(), ub_it );  
        copy_n( (*rb_iter)->transformed_idxs()->begin(), (*rb_iter)->num_idxs(), ti_it );  
 
      }

      auto orig_idxs_       = std::make_shared<const std::vector<std::string>>(orig_idxs);         
      auto orig_block_     = std::make_shared<const std::vector<std::string>>(orig_block);        
      auto unique_block_        = std::make_shared<const std::vector<std::string>>(unique_block);      
      auto transformed_idxs_ = std::make_shared<const std::vector<std::string>>(transformed_idxs);  

    }
 
   ~srbi_helper(){};

};

class split_range_block_info : public  range_block_info { 

  private : 
    srbi_helper tmp;
    std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks_;
  public :

     const int num_idxs_;
                              
    split_range_block_info( std::shared_ptr<std::vector<std::shared_ptr<range_block_info>>> range_blocks ) :
                      range_blocks_(range_blocks), tmp(srbi_helper(range_blocks_)), num_idxs_(tmp.num_idxs_),  
                      range_block_info( tmp.unique_, tmp.survives_, tmp.factors_, tmp.orig_block_,
                                        tmp.unique_block_, tmp.orig_idxs_, tmp.transformed_idxs_ ) {}  

 
   ~split_range_block_info(){};

    int num_idxs() { return num_idxs_ ; } 

    bool is_sparse( std::vector<std::vector<int>>& state_idxs ) { 
      std::vector<std::shared_ptr<range_block_info>>::iterator rb_iter =  range_blocks_->begin();
      for ( std::vector<std::vector<int>>::iterator si_iter = state_idxs.begin(); si_iter != state_idxs.end(); si_iter++ ){
         if ( (*rb_iter)->is_sparse(*si_iter) ) 
           return true;      
         rb_iter++;
      }
      return false; 
    }

};


#endif
