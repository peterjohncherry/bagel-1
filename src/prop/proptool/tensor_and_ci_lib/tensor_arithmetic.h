#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>

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

     static void
     add_tensors( std::shared_ptr<SMITH::Tensor_<DataType>> tens_target, std::shared_ptr<SMITH::Tensor_<DataType>> tens_summand, DataType factor );

     static void
     add_list_of_reordered_tensors( std::shared_ptr<SMITH::Tensor_<DataType>>& target,
                                    std::shared_ptr<SMITH::Tensor_<DataType>>& summand_orig_order,
                                    std::vector<std::vector<int>>& summand_reorderings,
                                    std::vector<DataType>& summand_factors                       );

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
     reorder_block_Tensor(std::shared_ptr<SMITH::Tensor_<DataType>> T_in_name, std::vector<int>& new_order);

     //Note the reordering assums column major ordering
     static std::unique_ptr<DataType[]>
     reorder_tensor_data( const DataType* orig_data,  std::vector<int>& new_order_vec, std::vector<SMITH::Index>& orig_index_blocks ) ;

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
     set_tensor_elems(std::shared_ptr<SMITH::Tensor_<DataType>> Tens, const std::vector<SMITH::IndexRange>& id_ranges, DataType elem_val );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_uniform_Tensor( const std::vector<SMITH::IndexRange>& T_id_ranges, DataType XX );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_test_tensor_row_major( const std::vector<SMITH::IndexRange>& T_id_ranges );

     static std::shared_ptr<SMITH::Tensor_<DataType>>
     get_test_tensor_column_major( const std::vector<SMITH::IndexRange>& T_id_ranges );

     static void 
     gemm( char op1, char op2, int size_i, int size_l, int size_j, DataType* A_data, DataType* B_data, DataType* C_data, DataType alpha, DataType beta ); 

     static void
     gemv( char op1, int size_i, int size_j, DataType* A_data, DataType* B_data, DataType* C_data, DataType alpha, DataType beta ); 
 
     static void 
     scaler(int T1_block_size, DataType T2_data_ptr, DataType* Tout_data_ptr);  

     static DataType
     dot_product( size_t vec_size, DataType* v1, DataType* v2); 

     static std::shared_ptr<SMITH::Tensor_<DataType>> 
     divide_tensors( std::shared_ptr<SMITH::Tensor_<DataType>> T1, std::shared_ptr<SMITH::Tensor_<DataType>> T2 );

     static void
     divide_tensors_in_place( std::shared_ptr<SMITH::Tensor_<DataType>> T1, std::shared_ptr<SMITH::Tensor_<DataType>> T2 );
    
     static void
     add_tensor_along_trace( std::shared_ptr<SMITH::Tensor_<DataType>>& t_target, std::shared_ptr<SMITH::Tensor_<DataType>>& t_summand,
                             std::vector<int>& summand_pos, DataType factor );

     static void
     ax_plus_y( int array_length, DataType factor,  DataType* target_ptr, DataType* summand_ptr );

     static void
     invert_matrix_general( int nrows, int ncols, DataType* data_ptr );

     static void 
     half_inverse_matrix_hermitian( int nrows, std::unique_ptr<DataType[]>& orig_data_ptr );

     static void
     diagonalize_matrix_hermitian( int nrows, DataType* orig_data_ptr, double* eigenvalues_ptr ); 

     static void
     zero_all_but_block( std::shared_ptr<SMITH::Tensor_<DataType>>& tens, const std::vector<SMITH::IndexRange>& nonzero_block );


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

template<> void Tensor_Arithmetic<std::complex<double>>::ax_plus_y( int array_length, std::complex<double> factor, std::complex<double>* target_ptr, std::complex<double>* summand_ptr );

template<> void Tensor_Arithmetic<double>::ax_plus_y( int array_length, double factor, double* target_ptr, double* summand_ptr );

template<> void Tensor_Arithmetic<double>::diagonalize_matrix_hermitian( int nrows, double* orig_data_ptr, double* eigenvalues_ptr ); 

template<> void Tensor_Arithmetic<std::complex<double>>::diagonalize_matrix_hermitian( int nrows, std::complex<double>* orig_data_ptr, double* eigenvalues_ptr ); 

template<> void Tensor_Arithmetic<double>::invert_matrix_general( int nrows, int ncols, double* data_ptr );

template<> void Tensor_Arithmetic<double>::half_inverse_matrix_hermitian( int nrows, std::unique_ptr<double[]>& orig_data_ptr ); 

//template<> void Tensor_Arithmetic<double>::put_sub_tensor( std::shared_ptr<SMITH::Tensor_<double>> Tens1, std::shared_ptr<SMITH::Tensor_<double>> Tens2 );

//template<> void Tensor_Arithmetic<std::complex<double>>::put_sub_tensor( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens1, std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens2 );
}
}
#endif
