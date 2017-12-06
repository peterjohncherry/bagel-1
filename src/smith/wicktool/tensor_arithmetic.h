#ifndef __SRC_SMITH_WICKTOOL_TENSOR_ARITH_H
#define __SRC_SMITH_WICKTOOL_TENSOR_ARITH_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic_utils.h>

//Seems ridiculous as all members are static. However, it is simpler to template a class than a name space.
//TODO At the moment I assume that all indexes running similar ranges are defined using IndexRange objects 
//     with similar range members. If this is not true, it will impact the id_offsets, and also the block pos
//     in the contractions.
//     Whilst this assumption hold true for the current method it may well not if I start to play with the
//     block size (e.g. set active blocksize = 1 for sigmas), and so it must be fixed.

namespace bagel {
namespace SMITH { 

namespace Tensor_Arithmetic { 

template<class DataType>
class Tensor_Arithmetic { 

    public: 
    Tensor_Arithmetic(){};
    ~Tensor_Arithmetic(){};

     static DataType
     sum_tensor_elems( std::shared_ptr<Tensor_<DataType>> Tens_in) ;

     static DataType
     trace_tensor__number_return( std::shared_ptr<Tensor_<DataType>> Tens_in) ;

     static std::shared_ptr<Tensor_<DataType>>
     trace_tensor__tensor_return( std::shared_ptr<Tensor_<DataType>> Tens_in) ;

     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor( std::shared_ptr<Tensor_<DataType>> Tens_in,  std::vector<int>& ctrs ); 

     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor( std::shared_ptr<Tensor_<DataType>> Tens_in,  std::pair<int,int> ctrs_pair ); 
     
     static std::shared_ptr<Tensor_<DataType>>
     contract_different_tensors( std::shared_ptr<Tensor_<DataType>> Tens1_in,
                                 std::shared_ptr<Tensor_<DataType>> Tens2_in,
                                 std::pair<int,int> ctr_todo                  );

     static std::shared_ptr<Tensor_<DataType>>
     contract_different_tensors( std::shared_ptr<Tensor_<DataType>> Tens1_in,
                                 std::shared_ptr<Tensor_<DataType>> Tens2_in,
                                 std::pair<std::vector<int>,std::vector<int>>& ctr_todo   );

     static std::shared_ptr<Tensor_<DataType>>
     direct_tensor_product( std::shared_ptr<Tensor_<DataType>> Tens1, std::shared_ptr<Tensor_<DataType>> Tens2 );

     static std::shared_ptr<Tensor_<DataType>>
     contract_tensor_with_vector( std::shared_ptr<Tensor_<DataType>> Tens1_in,
                                  std::shared_ptr<Tensor_<DataType>> Tens2_in,  
                                  int ctr_todo                                );

     static DataType
     contract_vectors( std::shared_ptr<Tensor_<DataType>> Tens1_in,  std::shared_ptr<Tensor_<DataType>> Tens2_in);

     static std::shared_ptr<Tensor_<DataType>>
     reorder_block_Tensor(std::shared_ptr<Tensor_<DataType>> T_in_name, std::shared_ptr<std::vector<int>> new_order);

     //Note the reordering assums column major ordering
     static std::unique_ptr<DataType[]>
     reorder_tensor_data( const DataType* orig_data,
                          std::shared_ptr<std::vector<int>> new_order_vec,
                          std::shared_ptr<std::vector<Index>> orig_index_blocks ) ;

     static std::unique_ptr<DataType[]>
     get_block_of_data( DataType* data_ptr , std::shared_ptr<std::vector<IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos) ;
     
     static DataType
     get_tensor_element( std::shared_ptr<Tensor_<DataType>> Tens, std::vector<int>& id_pos); 

     static void
     set_tensor_elems(std::shared_ptr<Tensor_<DataType>> Tens, DataType elem_val );
 
     static std::shared_ptr<Tensor_<DataType>>
     get_uniform_Tensor(std::shared_ptr<std::vector<IndexRange>> T_id_ranges, DataType XX );

     static std::shared_ptr<Tensor_<DataType>>
     get_test_Tensor_row_major( std::shared_ptr<std::vector<IndexRange>> T_id_ranges );

     static std::shared_ptr<Tensor_<DataType>>
     get_test_Tensor_column_major( std::shared_ptr<std::vector<IndexRange>> T_id_ranges );

}; 
 
}
}
}
#endif
