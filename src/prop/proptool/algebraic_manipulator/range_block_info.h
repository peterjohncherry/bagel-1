#ifndef __SRC_PROP_PROPTOOL_Range_Block_Info_H
#define __SRC_PROP_PROPTOOL_Range_Block_Info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
#include <src/prop/proptool/proputils.h>
 
class SRBI_Helper; 
class Split_Range_Block_Info; 
class Range_Block_Info : public std::enable_shared_from_this<Range_Block_Info> {
 
  friend SRBI_Helper;
  friend Split_Range_Block_Info;

  protected :
     std::pair<double,double> factors_; 
     std::shared_ptr<std::vector<int>> idxs_trans_;
     std::shared_ptr<std::vector<int>> aops_trans_;
     std::shared_ptr<std::vector<int>> rngs_trans_;

     std::shared_ptr<const std::vector<std::string>> orig_rngs_;
     std::shared_ptr<const std::vector<bool>> orig_aops_;

     std::set<std::vector<int>> sparsity_ ;
     std::string name_;

  public :

    std::vector<bool> allowed_contractions_;
    long unsigned int plus_pnum_;
    long unsigned int kill_pnum_;
    bool no_transition_;
    int num_idxs_;
    
    std::shared_ptr<const std::vector<std::string>> unique_block_;

    Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_rngs,   
                      std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                      std::shared_ptr<const std::vector<bool>> orig_aops, 
                      std::shared_ptr<std::vector<int>> idxs_trans,
                      std::shared_ptr<std::vector<int>> aops_trans,
                      std::shared_ptr<std::vector<int>> rngs_trans,
                      std::pair<double,double> factors  ); 

    Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_rngs,   
                      std::shared_ptr<const std::vector<bool>> orig_aops, 
                      std::shared_ptr<std::vector<int>> rngs_trans,
                      std::shared_ptr<std::vector<int>> aops_trans,
                      std::pair<double,double> factors  ); 

    ~Range_Block_Info(){};
    
    std::pair<double,double> factors() const { return factors_; } 
    double Re_factor() const { return factors_.first; } 
    double Im_factor() const { return factors_.second; } 
  
    std::string name() const { return name_; } 
     
    std::shared_ptr<const std::vector<std::string>> orig_rngs() { return orig_rngs_; } 
 
    std::shared_ptr<std::vector<int>> idxs_trans() const { return idxs_trans_; }
    std::shared_ptr<std::vector<int>> aops_trans() const { return aops_trans_; }
    std::shared_ptr<std::vector<int>> rngs_trans() const { return rngs_trans_; }

    virtual
    std::shared_ptr<Range_Block_Info> 
    transform( std::shared_ptr<const std::vector<std::string>> orig_rngs, std::shared_ptr<const std::vector<std::string>> orig_idxs, std::shared_ptr<const std::vector<bool>> orig_aops,
               std::vector<int>&  op_order, std::vector<char> op_trans ); 

    void transform_aops_rngs ( std::vector<bool>& aops, std::vector<char>& aops_rngs,  std::pair<double,double>& factors , char transformation );

    
    std::shared_ptr<Range_Block_Info> get_transformed_block( char op_trans, std::shared_ptr< const std::vector<bool> > aops );
};


class SRBI_Helper { 

  public :
    bool unique_;
    bool survives_;
    std::pair<double,double> factors_; 
    int num_idxs_; 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> rxnge_blocks_;
 
    std::shared_ptr<const std::vector<std::string>> orig_rngs_;
    std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    std::shared_ptr<const std::vector<bool>> orig_aops_;

    std::shared_ptr<std::vector<int>> idxs_trans_;
    std::shared_ptr<std::vector<int>> aops_trans_;
    std::shared_ptr<std::vector<int>> rngs_trans_;
                              
    SRBI_Helper( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks );

    void add_trans( std::shared_ptr<Split_Range_Block_Info> srbi, std::vector<int>&  op_order, std::vector<char> op_trans );
   ~SRBI_Helper(){};
 
};
 

class Split_Range_Block_Info : public  Range_Block_Info { 

  private : 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;

  public :

    Split_Range_Block_Info(  std::shared_ptr<std::vector<std::string>> orig_rngs, std::shared_ptr<const std::vector<std::string>> orig_idxs,
                              std::shared_ptr<const std::vector<bool>> orig_aops, SRBI_Helper& helper ) : 
                              Range_Block_Info( orig_rngs, orig_idxs, orig_aops, helper.idxs_trans_, helper.aops_trans_, helper.rngs_trans_, helper.factors_ ), 
                              range_blocks_(helper.rxnge_blocks_) {} 
                     
    Split_Range_Block_Info(  std::shared_ptr<std::vector<std::string>> orig_rngs, std::shared_ptr<const std::vector<bool>> orig_aops, SRBI_Helper& helper ) : 
                             Range_Block_Info( orig_rngs, orig_aops, helper.aops_trans_, helper.rngs_trans_, helper.factors_ ), 
                             range_blocks_(helper.rxnge_blocks_) {} 
                     
   ~Split_Range_Block_Info(){};

    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){ return range_blocks_ ;} 
    std::shared_ptr<Range_Block_Info> range_blocks(int ii){ return range_blocks_->at(ii) ;} 

    std::shared_ptr<Range_Block_Info> 
    transform( std::shared_ptr<const std::vector<std::string>> orig_rngs, std::shared_ptr<const std::vector<std::string>> orig_idxs, std::shared_ptr<const std::vector<bool>> orig_aops,
               std::vector<int>&  op_order, std::vector<char> op_trans );

    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs );
};

#endif
