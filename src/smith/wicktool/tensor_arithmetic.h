#ifndef __SRC_SMITH_WICKTOOL_TENSOR_ARITH_H
#define __SRC_SMITH_WICKTOOL_TENSOR_ARITH_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic_utils.h>

//Seems ridiculous as all members are static. However, eventually it will be templated, and it is simpler template a class than a name space.

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

     //only for two contracted indexes
     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor(std::shared_ptr<Tensor_<DataType>> Tens_in, std::pair<int,int> ctr_todo) ;
     
     //for an arbitrary number of contracted indexes
     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor(std::shared_ptr<Tensor_<DataType>> Tens_in , std::vector<int>& contracted_index_positions) ;
     
     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor_new( std::shared_ptr<Tensor_<DataType>> Tens_in,  std::vector<int>& ctrs ); 

     static std::shared_ptr<Tensor_<DataType>>
     contract_on_same_tensor_new( std::shared_ptr<Tensor_<DataType>> Tens_in,  std::pair<int,int>& ctrs_pair ); 

     static std::shared_ptr<Tensor_<DataType>>
     contract_different_tensors( std::shared_ptr<Tensor_<DataType>> Tens1_in, std::shared_ptr<Tensor_<DataType>> Tens2_in, std::pair<int,int> ctr_todo );
     
     static std::shared_ptr<Tensor_<DataType>>
     reorder_block_Tensor(std::shared_ptr<Tensor_<DataType>> T_in_name, std::shared_ptr<std::vector<int>> new_order);

     static std::unique_ptr<DataType[]>
     get_block_of_data( DataType* data_ptr , std::shared_ptr<std::vector<IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos) ;
     
     static void set_tensor_elems(std::shared_ptr<Tensor_<DataType>> Tens, DataType elem_val );
 
     static std::shared_ptr<Tensor_<DataType>>
     get_uniform_Tensor(std::shared_ptr<std::vector<IndexRange>> T_id_ranges, DataType XX );

     static DataType contract_vectors( std::shared_ptr<Tensor_<DataType>> Tens1_in,  std::shared_ptr<Tensor_<DataType>> Tens2_in);

     static std::unique_ptr<DataType[]>
     reorder_tensor_data( const DataType* orig_data, std::shared_ptr<std::vector<int>> new_order_vec, std::shared_ptr<std::vector<Index>> orig_index_blocks ) ;

     static std::shared_ptr<Tensor_<DataType>> contract_tensor_with_vector( std::shared_ptr<Tensor_<DataType>> Tens1_in, std::shared_ptr<Tensor_<DataType>> Tens2_in,  
                                                                            std::pair<int,int> ctr_todo);

     static std::shared_ptr<Tensor_<DataType>>
     get_test_Tensor_row_major( std::shared_ptr<std::vector<IndexRange>> T_id_ranges );

     static std::shared_ptr<Tensor_<DataType>>
     get_test_Tensor_column_major( std::shared_ptr<std::vector<IndexRange>> T_id_ranges );


     static std::shared_ptr<Tensor_<DataType>>
     contract_different_tensors_column_major( std::shared_ptr<Tensor_<DataType>> Tens1_in, std::shared_ptr<Tensor_<DataType>> Tens2_in, std::pair<int,int> ctr_todo );

}; 
 
}
}
}
#endif
