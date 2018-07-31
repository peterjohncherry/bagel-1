#ifndef __SRC_PROP_PROPTOOL_Range_Block_Info_H
#define __SRC_PROP_PROPTOOL_Range_Block_Info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
#include <set> 
#include <vector>
#include <memory>
#include <src/prop/proptool/algebraic_manipulator/op_info.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
 
class Split_Range_Block_Info; 
class Range_Block_Info : public std::enable_shared_from_this<Range_Block_Info> {
 
  friend Split_Range_Block_Info;

  private :
    std::shared_ptr<std::vector<int>> state_ids_; 

  protected :
    
    //NOTE : These factors correspond to very different things, and must be handled seperately  
    std::pair<double,double> ReIm_factors_;  // Factor arising from parts due to hermitian conjugation etc.,
    std::pair<double,double> factors_;       // Original block factor

    std::shared_ptr<std::vector<int>> idxs_trans_;         

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
    
    std::string op_state_name_;

    std::shared_ptr<Op_Info> op_info_;

    std::shared_ptr<const std::vector<std::string>> orig_rngs_;
    std::shared_ptr<Range_Block_Info> unique_block_;

    Range_Block_Info( std::pair<double,double> factors, std::pair<double,double> ReIm_factors, std::shared_ptr<Op_Info>& op_info ) :
                      factors_(factors), ReIm_factors_(ReIm_factors), op_info_(op_info) {} 

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

    void set_transition_vars(const std::vector<bool>& aops );

    std::pair<double,double> factors() const { return factors_; }
    double Re_factor() const { return factors_.first; }
    double Im_factor() const { return factors_.second; }

    std::string name() const { return name_; }
    std::string full_op_name() const { return full_op_name_; }

    std::shared_ptr<const std::vector<std::string>> orig_rngs() { return orig_rngs_; }
    std::shared_ptr< std::vector<char>> orig_rngs_ch() { return orig_rngs_ch_; }
 
    std::shared_ptr<std::vector<int>> idxs_trans() const { return idxs_trans_; }

    std::shared_ptr<Op_Info> op_info() { return op_info_; } 

    virtual std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){
      return std::make_shared<std::vector<std::shared_ptr<Range_Block_Info>>> ( 1, shared_from_this() );
    }

    virtual std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_canonical(){
      return range_blocks();
    }

};

class Split_Range_Block_Info : public  Range_Block_Info { 

  //private : 

  public :
    bool canonical_;
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_;
    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_canonical_;

    Split_Range_Block_Info( std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks, std::vector<int>& cml_sizes, 
                            std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<Op_Info> op_info,  
                            std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>>& split_ranges );
 
   ~Split_Range_Block_Info(){};

    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks(){ return range_blocks_ ;} 
    std::shared_ptr<Range_Block_Info> range_blocks(int ii){ return range_blocks_->at(ii) ;} 

    std::shared_ptr<std::vector<std::shared_ptr<Range_Block_Info>>> range_blocks_canonical(){ return range_blocks_canonical_; }
    bool is_sparse( const std::shared_ptr<std::vector<std::vector<int>>> state_idxs );
    
};

#endif
