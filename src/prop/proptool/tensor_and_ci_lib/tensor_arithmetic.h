#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_H

#include <src/prop/proptool/proputils.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>

//Seems ridiculous as all members are static. However, it is simpler to template a class than a name space.
//TODO At the moment I assume that all indexes running similar ranges are defined using IndexRange objects 
//     with similar range members. If this is not true, it will impact the id_offsets, and also the block pos
//     in the contractions.
//     Whilst this assumption hold true for the current method it may well not if I start to play with the
//     block size (e.g. set active blocksize = 1 for sigmas), and so it must be fixed.

namespace bagel {

namespace Tensor_Arithmetic { 

template<class DataType>
class Tensor_Arithmetic { 

    public: 
    Tensor_Arithmetic(){};
    ~Tensor_Arithmetic(){};

     static DataType
     sum_tensor_elems( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in) ;

     static DataType
     trace_tensor__number_return( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in) ;

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     trace_tensor__tensor_return( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in) ;

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     sum_over_idxs( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in, std::vector<int>& summed_idxs_pos); 

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     contract_on_same_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::vector<int>& ctrs ); 

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     contract_on_same_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::pair<int,int> ctrs_pair ); 
     
     static std::shared_ptr<SMITH::Tensor_<DataType>>
     contract_different_tensors( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1_in,
                                 std::shared_ptr<SMITH::Tensor_<DataType>> Tens2_in,
                                 std::pair<int,int> ctr_todo                  );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     contract_different_tensors( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1_in,
                                 std::shared_ptr<SMITH::Tensor_<DataType>> Tens2_in,
                                 std::pair<std::vector<int>,std::vector<int>>& ctr_todo   );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     direct_tensor_product( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1, std::shared_ptr<SMITH::Tensor_<DataType>> Tens2 );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     contract_tensor_with_vector( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1_in,
                                  std::shared_ptr<SMITH::Tensor_<DataType>> Tens2_in,  
                                  int ctr_todo                                );

     static DataType
     contract_vectors( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1_in,  std::shared_ptr<SMITH::Tensor_<DataType>> Tens2_in);

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     reorder_block_Tensor(std::shared_ptr<SMITH::Tensor_<DataType>> T_in_name, std::shared_ptr<std::vector<int>> new_order);

     //Note the reordering assums column major ordering
     static std::unique_ptr<DataType[]>
     reorder_tensor_data( const DataType* orig_data,
                          std::shared_ptr<std::vector<int>> new_order_vec,
                          std::shared_ptr<std::vector<SMITH::Index>> orig_index_blocks ) ;


     static std::unique_ptr<DataType[]>
     get_block_of_data( DataType* data_ptr , std::shared_ptr<std::vector<SMITH::IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos) ;
     
     static DataType
     get_tensor_element( std::shared_ptr<SMITH::Tensor_<DataType>> Tens, std::vector<int>& id_pos); 

     static void
     set_tensor_elems(std::shared_ptr<SMITH::Tensor_<DataType>> Tens, DataType elem_val );
 
     static void
     put_sub_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1, std::shared_ptr<SMITH::Tensor_<DataType>> Tens2 );

     static void
     put_tensor_range_block( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1, std::shared_ptr<SMITH::Tensor_<DataType>> Tens2, std::vector<SMITH::IndexRange>& id_ranges );

     static void
     put_reordered_range_block( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1, std::shared_ptr<SMITH::Tensor_<DataType>> Tens2,
                                std::vector<SMITH::IndexRange>& id_ranges, std::shared_ptr<std::vector<int>> new_order     );


     static void 
     put_reordered_range_block( std::shared_ptr<SMITH::Tensor_<DataType>> T1, std::vector<SMITH::IndexRange>& id_ranges_T1,
                                std::shared_ptr<SMITH::Tensor_<DataType>> T2, std::vector<SMITH::IndexRange>& id_ranges_T2,
                                std::shared_ptr<std::vector<int>> new_order );

     static void
     set_tensor_elems(std::shared_ptr<SMITH::Tensor_<DataType>> Tens, std::vector<SMITH::IndexRange>& id_ranges, DataType elem_val );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_uniform_Tensor(std::shared_ptr<std::vector<SMITH::IndexRange>> T_id_ranges, DataType XX );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_test_Tensor_row_major( std::shared_ptr<std::vector<SMITH::IndexRange>> T_id_ranges );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_test_Tensor_column_major( std::shared_ptr<std::vector<SMITH::IndexRange>> T_id_ranges );

     static void 
     gemm( char op1, char op2, int size_i, int size_l, int size_j, DataType* A_data, DataType* B_data, DataType* C_data, DataType alpha, DataType beta ); 

     static void
     gemv( char op1, int size_i, int size_j, DataType* A_data, DataType* B_data, DataType* C_data, DataType alpha, DataType beta ); 
 
     static void 
     scaler(int T1_block_size, DataType T2_data_ptr, DataType* Tout_data_ptr);  

     static DataType
     dot_product( size_t vec_size, DataType* v1, DataType* v2); 
}; 
 
template<> void Tensor_Arithmetic<double>::gemm ( char op1, char op2, int size_i, int size_l, int size_j, 
                                                  double* A_data, double* B_data, double* C_data, double alpha, double beta );

template<> void Tensor_Arithmetic<std::complex<double>>::gemm( char op1, char op2, int size_i, int size_l, int size_j, 
                                                          std::complex<double>* A_data, std::complex<double>* B_data, std::complex<double>* C_data,
                                                          std::complex<double> alpha, std::complex<double> beta );

template<> void Tensor_Arithmetic<double>::gemv( char op1, int size_i, int size_j, double* A_data, double* B_data, double* C_data, double alpha, double beta ); 

template<> void Tensor_Arithmetic<std::complex<double>>::gemv( char op1, int size_i, int size_j,
                                                               std::complex<double>* A_data, std::complex<double>* B_data, std::complex<double>* C_data, 
                                                               std::complex<double> alpha, std::complex<double> beta );

template<> void Tensor_Arithmetic<double>::scaler( int T1_block_size, double T2_data_ptr, double* Tout_data_ptr);  

template<> void Tensor_Arithmetic<std::complex<double>>::scaler( int T1_block_size, std::complex<double> T2_data_ptr, std::complex<double>* Tout_data_ptr);  

template<> double Tensor_Arithmetic<double>::dot_product( size_t vec_size, double* v1, double* v2);

template<> std::complex<double> Tensor_Arithmetic<std::complex<double>>::dot_product( size_t vec_size, std::complex<double>* v1, std::complex<double>* v2); 
}
}
#endif
