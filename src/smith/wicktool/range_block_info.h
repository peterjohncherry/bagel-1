#ifndef __SRC_SMITH_range_block_info_H
#define __SRC_SMITH_range_block_info_H

// On reflection, all this information is not needed, only the original ranges (maybe also idxs and aops), and the transformation vector.
// It could be done this way, but a lot will need to be fed into the gamma generator.
// Factor is still necessary though.
class range_block_info {

  // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
  private :
    const bool is_unique_;
    const bool survives_;
    const std::pair<double,double> factors_; 
    const std::shared_ptr<const std::vector<std::string>> orig_block_;
    const std::shared_ptr<const std::vector<std::string>> unique_block_;
    const std::shared_ptr<const std::vector<std::string>> orig_idxs_;
    const std::shared_ptr<const std::vector<std::string>> transformed_idxs_;
                              
  public :
    range_block_info( bool is_unique,
                      bool survives,
                      std::pair<double,double> factors, 
                      std::shared_ptr<const std::vector<std::string>> orig_block,   
                      std::shared_ptr<const std::vector<std::string>> unique_block, 
                      std::shared_ptr<const std::vector<std::string>> orig_idxs,   
                      std::shared_ptr<const std::vector<std::string>> transformed_idxs) : 
                      is_unique_(is_unique), survives_(survives),factors_(factors),
                      orig_block_(orig_block), unique_block_(unique_block),
                      orig_idxs_(orig_idxs), transformed_idxs_(transformed_idxs) {};
    ~range_block_info(){};
    
    bool is_unique() const { return is_unique_ ; } //get<0>
    bool survives() const { return survives_ ; }
    
    std::pair<double,double> factors() const { return factors_; } // get<3>
    double Re_factor() const { return factors_.first; } 
    double Im_factor() const { return factors_.second; } 
    
    std::shared_ptr<const std::vector<std::string>> orig_idxs() const      { return orig_idxs_;      } 
    std::shared_ptr<const std::vector<std::string>> orig_block() const      { return orig_block_;      }// key 
    std::shared_ptr<const std::vector<std::string>> unique_block() const    { return unique_block_;    }// get<1>
    std::shared_ptr<const std::vector<std::string>> transformed_idxs() const { return transformed_idxs_;}// get<2>

};
class split_range_block_info {

  // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
  public :
    std::shared_ptr<std::vector<std::string>> op_names_;
    std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>> range_blocks_;
    std::shared_ptr<std::vector<std::pair<double,double>>> factors_;
    std::shared_ptr<std::vector<std::shared_ptr<const std::vector<int>>>> transformations_;
                              
    split_range_block_info( std::shared_ptr<std::vector<std::string>> op_names,   
                            std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> range_blocks,
                            std::shared_ptr<std::vector<std::pair<double,double>>> factors,
                            std::shared_ptr<std::vector< std::shared_ptr<const std::vector<int>>>> transformations):
                            op_names_(op_names), range_blocks_(range_blocks),
                            factors_(factors), transformations_(transformations) {};
    ~split_range_block_info(){};

};

#endif
