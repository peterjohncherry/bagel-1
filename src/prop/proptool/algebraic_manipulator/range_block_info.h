#ifndef __SRC_PROP_PROPTOOL_Range_Block_Info_H
#define __SRC_PROP_PROPTOOL_Range_Block_Info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
#include <src/prop/proptool/proputils.h>
class Range_Block_Info : public std::enable_shared_from_this<Range_Block_Info> {

  // TODO Ultimately, this should contain a vector defining the reordering of the aops from the standard, and a vector defining the reordering of the idxs,
  //      rather than the reordered vetors themselves.
  //      all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
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

    std::vector<bool> allowed_contractions_;
    long unsigned int plus_pnum_;
    long unsigned int kill_pnum_;

    Range_Block_Info( bool is_unique,
                      bool survives,
                      std::pair<double,double> factors, 
                      std::shared_ptr<const std::vector<std::string>> orig_block,   
                      std::shared_ptr<const std::vector<std::string>> unique_block, 
                      std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                      std::shared_ptr<const std::vector<std::string>> transformed_idxs,
                      std::shared_ptr<const std::vector<bool>> orig_aops, 
                      std::shared_ptr< std::map < char, long unsigned int>> range_prime_map  ); 

    ~Range_Block_Info(){};
    
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
    virtual std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){
        throw std::logic_error("Not a split_range_block; cannot get range_block_vec! Aborting!!" );
          return std::make_shared<std::vector<std::shared_ptr<Range_Block_Info>>>(1, shared_from_this()); }  


};

//This is an absurd hack....
class SRBI_Helper { 

  public :
    bool unique_;
    bool survives_;
    std::pair<double,double> factors_; 
    int num_idxs_; 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;
 
    std::shared_ptr<const std::vector<std::string>> orig_block_;
    std::shared_ptr<const std::vector<std::string>> unique_block_;

    std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    std::shared_ptr<const std::vector<std::string>> transformed_idxs_;

    std::shared_ptr<const std::vector<bool>> orig_aops_;
                              
    SRBI_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks );
   ~SRBI_Helper(){};
 
};
 
class Split_Range_Block_Info : public  Range_Block_Info { 

  private : 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;

  public :

    Split_Range_Block_Info( SRBI_Helper& helper, std::shared_ptr<std::map< char, long unsigned int>> range_prime_map  ) :
                      Range_Block_Info( helper.unique_, helper.survives_, helper.factors_, helper.orig_block_,
                                        helper.unique_block_, helper.orig_idxs_, helper.transformed_idxs_, helper.orig_aops_, range_prime_map ),
                                        range_blocks_(helper.range_blocks_) {}  
   ~Split_Range_Block_Info(){};

    int num_idxs() { return num_idxs_ ; } 

    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){ return range_blocks_ ;} 
    std::shared_ptr<Range_Block_Info> range_blocks(int ii){ return range_blocks_->at(ii) ;} 

    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs );
};
 
class Range_BlockX_Info : public std::enable_shared_from_this<Range_BlockX_Info> {

  protected :
     std::pair<double,double> factors_; 
     std::shared_ptr<const std::vector<std::string>> orig_idxs_;
     std::shared_ptr<const std::vector<std::string>> orig_rngs_;
     std::shared_ptr<const std::vector<bool>>   orig_aops_;
     std::shared_ptr<std::vector<int>> idxs_trans_;
     std::shared_ptr<std::vector<int>> aops_trans_;
     std::shared_ptr<std::vector<int>> rngs_trans_;

     int num_idxs_;                      
     std::string orig_name_;
     std::string transformed_name_;
     std::string TensOp_name_;

     std::set<std::vector<int>> sparsity_ ;

  public :

    std::vector<bool> allowed_contractions_;
    long unsigned int plus_pnum_;
    long unsigned int kill_pnum_;
    bool no_transition_;

    Range_BlockX_Info( std::shared_ptr<const std::vector<std::string>> orig_block,   
                       std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                       std::shared_ptr<const std::vector<bool>> orig_aops, 
                       std::shared_ptr<std::vector<int>> idxs_trans,
                       std::shared_ptr<std::vector<int>> aops_trans,
                       std::shared_ptr<std::vector<int>> rngs_trans,
                       std::pair<double,double> factors  ); 

    ~Range_BlockX_Info(){};
    
    std::pair<double,double> factors() const { return factors_; } 
    double Re_factor() const { return factors_.first; } 
    double Im_factor() const { return factors_.second; } 
    
    std::shared_ptr<const std::vector<std::string>> orig_rngs() const { return orig_rngs_; }
    std::shared_ptr<const std::vector<std::string>> orig_idxs() const { return orig_idxs_; }
    std::shared_ptr<const std::vector<bool>> orig_aops() const { return orig_aops_; }
   
    std::shared_ptr<std::vector<int>> idxs_trans() const { return idxs_trans_; }
    std::shared_ptr<std::vector<int>> aops_trans() const { return aops_trans_; }
    std::shared_ptr<std::vector<int>> rngs_trans() const { return rngs_trans_; }

    void add_sparse( std::vector<int>& state_idxs ) { sparsity_.emplace(state_idxs); return; } // defines this range block as being sparse for input states 

    int num_idxs() const { return num_idxs_; } 

    std::string orig_name() const { return orig_name_; }


    // returns true if this block is sparse for input states
    virtual bool is_sparse( std::vector<int>& state_idxs ) { return ( sparsity_.find(state_idxs) != sparsity_.end() ); } 
    virtual bool is_sparse( const std::vector<int>& state_idxs ) { return ( sparsity_.find(state_idxs) != sparsity_.end() ); }  
    virtual std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks(){
        throw std::logic_error("Not a split_range_block; cannot get range_block_vec! Aborting!!" );
          return std::make_shared<std::vector<std::shared_ptr<Range_BlockX_Info>>>(1, shared_from_this()); }  

};


class SRBIX_Helper { 

  public :
    bool unique_;
    bool survives_;
    std::pair<double,double> factors_; 
    int num_idxs_; 
    std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> rxnge_blocks_;
 
    std::shared_ptr<const std::vector<std::string>> orig_rngs_;
    std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    std::shared_ptr<const std::vector<bool>> orig_aops_;

    std::shared_ptr<std::vector<int>> idxs_trans_;
    std::shared_ptr<std::vector<int>> aops_trans_;
    std::shared_ptr<std::vector<int>> rngs_trans_;
                              
    SRBIX_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks );
   ~SRBIX_Helper(){};
 
};
 

class SplitX_Range_Block_Info : public  Range_BlockX_Info { 

  private : 
    std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks_;

  public :

    SplitX_Range_Block_Info( SRBIX_Helper& helper ) :
		             Range_BlockX_Info( helper.orig_rngs_, helper.orig_idxs_, helper.orig_aops_,
                                                helper.idxs_trans_,  helper.aops_trans_, helper.rngs_trans_, helper.factors_ ),
                                                range_blocks_(helper.rxnge_blocks_) {}  
   ~SplitX_Range_Block_Info(){};

    int num_idxs() { return num_idxs_ ; } 

    std::shared_ptr<std::vector<std::shared_ptr<Range_BlockX_Info>>> range_blocks(){ return range_blocks_ ;} 
    std::shared_ptr<Range_BlockX_Info> range_blocks(int ii){ return range_blocks_->at(ii) ;} 

    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs );
};



#endif
