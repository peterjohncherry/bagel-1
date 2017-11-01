#ifndef __SRC_SMITH_Equation_Computer_Utils_H
#define __SRC_SMITH_Equation_Computer_Utils_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/smith/wicktool/WickUtils.h>
namespace bagel {
namespace SMITH { 


namespace Equation_Computer_Utils {  

     std::shared_ptr<std::vector<int>> get_CTens_strides( std::shared_ptr<std::vector<int>> range_sizes, int ctr1 , int ctr2 ) ;
     std::shared_ptr<std::vector<int>> get_CTens_strides( std::vector<int>& range_sizes, int ctr1 , int ctr2 ) ;
     std::shared_ptr<std::vector<int>> get_CTens_strides( std::vector<int>& range_sizes, std::vector<int>& ctr_idxs_pos );

     std::shared_ptr<std::vector<int>> get_Tens_strides ( std::vector<int>& range_sizes) ;
 
     std::shared_ptr<std::vector<Index>> get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> old_ids) ;
     std::shared_ptr<std::vector<Index>> get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::shared_ptr<std::vector<std::shared_ptr<IndexRange>>> old_ids) ;
     std::shared_ptr<std::vector<Index>> get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::vector<IndexRange>& id_ranges);
 
     std::shared_ptr<std::vector<int>> get_num_index_blocks_vec(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec) ;
     std::shared_ptr<std::vector<int>> get_num_index_blocks_vec(std::shared_ptr<std::vector<IndexRange>> rngvec) ;
     std::vector<int>                  get_num_index_blocks_vec(std::vector<IndexRange>& rngvec);

     std::shared_ptr<std::vector<int>>     get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec, int skip_id) ;
     std::shared_ptr<std::vector<int>>     get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec) ;
 
     std::shared_ptr<std::vector<size_t>>  get_sizes(std::shared_ptr<std::vector<Index>> Idvec, int skip_id);
     std::shared_ptr<std::vector<size_t>>  get_sizes(std::shared_ptr<std::vector<Index>> Idvec);
     std::vector<int>                      get_sizes(std::vector<Index>& Idvec) ;
 
     std::shared_ptr<std::vector<int>>  get_range_lengths( std::shared_ptr<std::vector<IndexRange>> Id_ranges ) ;
     std::shared_ptr<std::vector<int>>  get_range_lengths( std::vector<IndexRange>& Id_ranges ) ;
 
     std::shared_ptr<std::vector<int>>  put_ctr_at_front(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);
     std::shared_ptr<std::vector<int>>  put_ctr_at_back(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);

     void Print_Tensor(std::shared_ptr<Tensor_<double>> Tens) ;
     
     size_t get_block_size(std::vector<Index>::iterator beginpos, std::vector<Index>::iterator endpos  ); 

     size_t get_unc_block_size( std::vector<Index>& idvec, std::pair<int,int> ctr ) ;

     void check_contracted_indexes( std::vector<IndexRange>&  idx_block, std::vector<int>& contracted_index_positions );
}
}
}
#endif 
