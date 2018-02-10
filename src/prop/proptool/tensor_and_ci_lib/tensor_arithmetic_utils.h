#ifndef __SRC_PROP_PROPTOOL_Tensor_Arithmetic_Utils_H
#define __SRC_PROP_PROPTOOL_Tensor_Arithmetic_Utils_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/proputils.h>
namespace bagel {

namespace Tensor_Arithmetic_Utils {  

     std::shared_ptr<std::vector<int>> get_CTens_strides( std::shared_ptr<std::vector<int>> range_sizes, int ctr1 , int ctr2 ) ;
     std::shared_ptr<std::vector<int>> get_CTens_strides( std::vector<int>& range_sizes, int ctr1 , int ctr2 ) ;
     std::shared_ptr<std::vector<int>> get_CTens_strides( std::vector<int>& range_sizes, std::vector<int>& ctr_idxs_pos );

     std::shared_ptr<std::vector<int>> get_Tens_strides ( std::vector<int>& range_sizes);
     std::shared_ptr<std::vector<int>> get_Tens_strides_column_major( std::vector<int>& range_sizes); 
     std::shared_ptr<std::vector<int>> get_Tens_strides_row_major( std::vector<int>& range_sizes);
 
     std::shared_ptr<std::vector<SMITH::Index>>
     get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>> old_ids) ;

     std::shared_ptr<std::vector<SMITH::Index>>
     get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::shared_ptr<std::vector<std::shared_ptr<SMITH::IndexRange>>> old_ids) ;

     std::shared_ptr<std::vector<SMITH::Index>>
     get_rng_blocks(std::shared_ptr<std::vector<int>> block_pos, std::vector<SMITH::IndexRange>& id_ranges);
 
     std::shared_ptr<std::vector<int>> get_num_index_blocks_vec(std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>> rngvec) ;
     std::shared_ptr<std::vector<int>> get_num_index_blocks_vec(std::shared_ptr<std::vector<SMITH::IndexRange>> rngvec) ;
     std::vector<int>                  get_num_index_blocks_vec(std::vector<SMITH::IndexRange>& rngvec);

     std::shared_ptr<std::vector<int>>     get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>> rngvec, int skip_id) ;
     std::shared_ptr<std::vector<int>>     get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>> rngvec) ;
 
     std::shared_ptr<std::vector<size_t>>  get_sizes(std::shared_ptr<std::vector<SMITH::Index>> Idvec, int skip_id);
     std::shared_ptr<std::vector<size_t>>  get_sizes(std::shared_ptr<std::vector<SMITH::Index>> Idvec);
     std::vector<int>                      get_sizes(std::vector<SMITH::Index>& Idvec) ;
 
     std::shared_ptr<std::vector<int>>  get_range_lengths( std::shared_ptr<std::vector<SMITH::IndexRange>> Id_ranges ) ;
     std::shared_ptr<std::vector<int>>  get_range_lengths( std::vector<SMITH::IndexRange>& Id_ranges ) ;
 
     std::shared_ptr<std::vector<int>>  put_ctr_at_front(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);
     std::shared_ptr<std::vector<int>>  put_ctr_at_back(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);
     
     void put_ctrs_at_front( std::vector<int>& id_pos, std::vector<int>& ctr_todo);
     void put_ctrs_at_front( std::vector<int>& id_pos, std::pair<int,int>& ctr_todo);

     void put_reversed_ctrs_at_front( std::vector<int>& id_pos, std::vector<int> ctr_todo);

     void put_ctrs_at_back( std::vector<int>& id_pos, std::vector<int>& ctr_todo);
     void put_ctrs_at_back( std::vector<int>& id_pos, std::pair<int,int>& ctr_todo);

     void Print_Tensor(std::shared_ptr<SMITH::Tensor_<double>> Tens, std::string name = "") ;
     void Print_Tensor_row_major( std::shared_ptr<SMITH::Tensor_<double>> Tens, std::string name = "") ;

     void Print_Vector_Tensor_Format( std::shared_ptr<SMITH::Tensor_<double>> VecIn, std::string name );

     size_t get_block_size(std::vector<SMITH::Index>::iterator beginpos, std::vector<SMITH::Index>::iterator endpos  ); 

     size_t get_unc_block_size( std::vector<SMITH::Index>& idvec, std::pair<int,int> ctr ) ;

     std::shared_ptr<std::vector<std::pair<size_t,size_t>>> get_block_start( std::shared_ptr<std::vector<SMITH::IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos ) ;

     std::shared_ptr<std::vector<std::vector<int>>> get_block_offsets(std::vector<SMITH::IndexRange>&  ranges ) ;
     std::shared_ptr<std::vector<std::vector<int>>> get_block_offsets(std::shared_ptr<std::vector<SMITH::IndexRange>> ranges ) ;

     void check_contracted_indexes( std::vector<SMITH::IndexRange>&  idx_block, std::vector<int>& contracted_index_positions );

     template<typename DataType> 
     std::shared_ptr<SMITH::Tensor_<DataType>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::vector<SMITH::IndexRange>& range_names );

     template<> 
     std::shared_ptr<SMITH::Tensor_<double>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<double>> Tens_in,  std::vector<SMITH::IndexRange>& range_names );


     template<typename DataType> 
     std::shared_ptr<SMITH::Tensor_<DataType>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::vector<std::string>& range_names,
                                                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange> >> range_conversion_map );
     template<>
     std::shared_ptr<SMITH::Tensor_<double>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<double>> Tens_in,  std::vector<std::string>& range_names,
                                                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange> >> range_conversion_map );


     //messy, but I can't see a better way
     template<class vtype>
     std::shared_ptr<std::vector<vtype>> inverse_reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) {
     auto newvec = std::make_shared<std::vector<vtype>>(origvec->size());
     for( int ii = 0; ii != origvec->size(); ii++ )
         newvec->at(neworder->at(ii)) =  origvec->at(ii);
     return newvec;
     };

     template<class vtype>
     std::shared_ptr<std::vector<vtype>> inverse_reorder_vector(std::shared_ptr<std::vector<int>> neworder, std::vector<vtype>& origvec ) {
     auto newvec = std::make_shared<std::vector<vtype>>(origvec.size());
     for( int ii = 0; ii != origvec->size(); ii++ )
       newvec->at(neworder->at(ii)) =  origvec->at(ii);
     return newvec;
     };

     template<class vtype>
     std::shared_ptr<std::vector<vtype>> reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::vector<vtype>& origvec ) {
     auto newvec = std::make_shared<std::vector<vtype>>(origvec.size());
     for( int ii = 0; ii != origvec.size(); ii++ )
       newvec->at(ii) = origvec[neworder->at(ii)];
     return newvec;
     };

     template<class vtype>
     std::shared_ptr<std::vector<vtype>> reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) {
       auto newvec = std::make_shared<std::vector<vtype>>(origvec->size());
       for( int ii = 0; ii != origvec->size(); ii++ )
          newvec->at(ii) = origvec->at(neworder->at(ii));
       return newvec;
     };

}
}
#endif 
