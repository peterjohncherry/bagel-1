#ifndef __SRC_PROP_PROPTOOL_Range_Block_Info_H
#define __SRC_PROP_PROPTOOL_Range_Block_Info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/op_info.h>
 
class SRBI_Helper; 
class Split_Range_Block_Info; 
class Range_Block_Info : public std::enable_shared_from_this<Range_Block_Info> {
 
  friend SRBI_Helper;
  friend Split_Range_Block_Info;

  private :
    std::shared_ptr<std::vector<int>> state_ids_; 

  protected :
    
    std::pair<double,double> ReIm_factors_;  // Factor arising from parts due to hermitian conjugation etc.,

    std::pair<double,double> factors_; //  Original block factor

    std::shared_ptr<std::vector<int>> idxs_trans_;         
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_;  
    std::shared_ptr<std::vector<int>> aops_trans_;

    std::shared_ptr<std::vector<char>> orig_rngs_ch_;
    std::shared_ptr<const std::vector<bool>> orig_aops_;

    std::set<std::vector<int>> sparsity_ ;
    std::string full_op_name_;
    std::string name_;

  public :

    long unsigned int plus_pnum_;
    long unsigned int kill_pnum_;
    bool ci_sector_transition_;
    int num_idxs_;
    std::shared_ptr<Transformation> transform_;
    
    std::shared_ptr<const std::vector<std::string>> orig_rngs_;
    std::shared_ptr<Range_Block_Info> unique_block_;

    Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_block, 
                      std::shared_ptr<std::vector<int>> idxs_trans,  std::pair<double,double> factors, std::pair<double,double> ReIm_factors,
                      const std::vector<bool>& aops,  std::shared_ptr<Op_Info>& op_info );

    Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_block, std::shared_ptr<Range_Block_Info> unique_block, 
                      std::shared_ptr<Transformation> transform,  std::pair<double,double> factors, std::pair<double,double> ReIm_factors,
                      const std::vector<bool>& aops, std::shared_ptr<Op_Info>& op_info );

    Range_Block_Info( std::shared_ptr<const std::vector<std::string>> orig_block,  std::shared_ptr<Range_Block_Info> unique_block, 
                      std::shared_ptr<std::vector<int>> idxs_trans,  std::pair<double,double> factors, std::pair<double,double> ReIm_factors,
                      const std::vector<bool>& aops,  std::shared_ptr<Op_Info>& op_info );

    ~Range_Block_Info(){};
    
    std::pair<double,double> factors() const { return factors_; } 
    double Re_factor() const { return factors_.first; } 
    double Im_factor() const { return factors_.second; } 
  
    std::string name() const { return name_; } 
    std::string full_op_name() const { return full_op_name_; } 


    virtual
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){ 
      throw std::logic_error(" Should only call range_blocks() from Split_Range_Block" ); 
      std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> dummy; 
      return dummy;
    } 
     
    std::shared_ptr<const std::vector<std::string>> orig_rngs() { return orig_rngs_; } 
    std::shared_ptr< std::vector<char>> orig_rngs_ch() { return orig_rngs_ch_; } 
 
    std::shared_ptr<std::vector<int>> idxs_trans() const { return idxs_trans_; }
    std::shared_ptr<std::vector<int>> aops_trans() const { return aops_trans_; }

    std::shared_ptr<std::vector<int>> idxs_trans_inverse() { return idxs_trans_inverse_;  }  
};

class SRBI_Helper { 

  public :
    bool unique_;
    bool survives_;
    std::pair<double,double> factors_; 
    std::pair<double,double> ReIm_factors_; 
    int num_idxs_; 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;
 
    std::shared_ptr<const std::vector<std::string>> orig_rngs_;
    std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    std::shared_ptr<const std::vector<bool>> orig_aops_;

    std::shared_ptr<std::vector<int>> idxs_trans_;
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_;
    std::shared_ptr<std::vector<int>> aops_trans_;
                              
    std::shared_ptr<Range_Block_Info> unique_block_;

    SRBI_Helper( std::vector<std::shared_ptr<Range_Block_Info>>& range_blocks, std::vector<int>& cml_sizes, 
                 std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<Op_Info> multiop_info,  
                 std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>>& state_specific_split_ranges ); 
   ~SRBI_Helper(){};

    void add_trans( std::shared_ptr<Split_Range_Block_Info> srbi, std::vector<int>&  op_order, std::vector<char> op_trans );
 
};

class Split_Range_Block_Info : public  Range_Block_Info, std::enable_shared_from_this<Split_Range_Block_Info> { 

  private : 
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;

  public :

    Split_Range_Block_Info( const std::vector<bool>& aops, SRBI_Helper& helper, std::shared_ptr<Op_Info> op_info ) : 
                            Range_Block_Info( helper.orig_rngs_, helper.unique_block_, helper.idxs_trans_, helper.factors_, helper.ReIm_factors_, aops , op_info ),
                            range_blocks_(helper.range_blocks_) {}// this->unique_blockXX_ = helper.unique_blockXX_; } 
                     
   ~Split_Range_Block_Info(){};

    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){ return range_blocks_ ;} 
    std::shared_ptr<Range_Block_Info> range_blocks(int ii){ return range_blocks_->at(ii) ;} 

    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs );

};

#endif
