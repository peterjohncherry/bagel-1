#ifndef __SRC_PROP_PROPTOOL_Tensor_Arithmetic_Utils_H
#define __SRC_PROP_PROPTOOL_Tensor_Arithmetic_Utils_H
// #define __DEBUG_TENSOR_ARITHMETIC_UTILS
// #define __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE 
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
namespace bagel {

namespace Tensor_Arithmetic_Utils {  

     std::vector<size_t> get_strides_row_major( const std::vector<SMITH::Index>& block ); // C ordering

     std::vector<size_t> get_strides_column_major( const std::vector<SMITH::Index>& block ); // Fortran (SMITH::Tensor ordering)

     std::pair<std::vector<std::vector<size_t>>,std::vector<std::vector<size_t>>>
     get_block_start_ends( const std::vector<SMITH::IndexRange>& ranges );

     double sum_elems( std::unique_ptr<double[]>& some_data, size_t length );
     std::complex<double> sum_elems( std::unique_ptr<std::complex<double>[]>& some_data, size_t length );
     double sum_elem_norms( std::unique_ptr<double[]>& some_data, size_t length );
     double sum_elem_norms( std::unique_ptr<std::complex<double>[]>& some_data, size_t length );
 
     std::vector<int> get_Tens_strides ( const std::vector<int>& range_sizes);
     std::vector<int> get_Tens_strides_column_major( const std::vector<int>& range_sizes); 
 
     std::vector<SMITH::Index>
     get_rng_blocks(const std::vector<int>& block_pos, const std::vector<SMITH::IndexRange>& id_ranges);

     std::vector<int>  get_num_index_blocks_vec(std::vector<SMITH::IndexRange>& rngvec);

     std::shared_ptr<std::vector<size_t>> get_sizes_sp(const std::vector<SMITH::Index>& Idvec);

     std::vector<size_t> get_sizes(const std::vector<SMITH::Index>& Idvec);

     template< typename T = size_t >
     std::vector<T> get_sizes_m1(const std::vector<SMITH::Index>& Idvec); // routine for fvec_cycle (maxs = size-1)

     template<>
     std::vector<int> get_sizes_m1(const std::vector<SMITH::Index>& Idvec); // routine for fvec_cycle (maxs = size-1)
 
     std::vector<int>  get_range_lengths( const std::vector<SMITH::IndexRange>& id_ranges ) ;
 
     std::vector<int>  put_ctr_at_front(const std::vector<int>& orig_pos , int ctr_pos);
     std::vector<int>  put_ctr_at_back( const std::vector<int>& orig_pos , int ctr_pos);
     
     void put_ctrs_at_front( std::vector<int>& id_pos, std::vector<int>& ctr_todo);
     void put_ctrs_at_front( std::vector<int>& id_pos, std::pair<int,int>& ctr_todo);

     void put_reversed_ctrs_at_front( std::vector<int>& id_pos, std::vector<int> ctr_todo);

     void put_ctrs_at_back( std::vector<int>& id_pos, std::vector<int>& ctr_todo);

     size_t get_block_size(std::vector<SMITH::Index>::iterator beginpos, std::vector<SMITH::Index>::iterator endpos  ); 

     size_t get_unc_block_size( std::vector<SMITH::Index>& idvec, std::pair<int,int> ctr ) ;

     std::shared_ptr<std::vector<std::pair<size_t,size_t>>> get_block_start( std::shared_ptr<std::vector<SMITH::IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos ) ;

     std::shared_ptr<std::vector<std::vector<int>>> get_block_offsets_sp(const std::vector<SMITH::IndexRange>&  ranges ) ;

     std::vector<std::vector<int>> get_block_offsets(const std::vector<SMITH::IndexRange>& ranges ) ;

     void check_contracted_indexes( std::vector<SMITH::IndexRange>&  idx_block, std::vector<int>& contracted_index_positions );

     std::vector<SMITH::IndexRange> get_indexrange_from_dimension_vector( std::vector<size_t> dimensions, std::vector<size_t> maxblock  );

     template<typename DataType>
     void Print_Tensor(std::shared_ptr<SMITH::Tensor_<DataType>> Tens, std::string name = "") ;
     template<>
     void Print_Tensor(std::shared_ptr<SMITH::Tensor_<double>> Tens, std::string name) ;
     template<>
     void Print_Tensor(std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens, std::string name) ;

     template<typename DataType>
     void Print_Tensor_row_major( std::shared_ptr<SMITH::Tensor_<DataType>> Tens, std::string name = "") ;
     template<>
     void Print_Tensor_row_major( std::shared_ptr<SMITH::Tensor_<double>> Tens, std::string name ) ;
     template<>
     void Print_Tensor_row_major( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens,  std::string name ) ;

     template<typename DataType>
     void Print_Vector_Tensor_Format( std::shared_ptr<SMITH::Tensor_<DataType>> VecIn, std::string name );
     template<>
     void Print_Vector_Tensor_Format( std::shared_ptr<SMITH::Tensor_<double>> VecIn, std::string name );
     template<>
     void Print_Vector_Tensor_Format( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> VecIn, std::string name );

     template<typename DataType> 
     std::shared_ptr<SMITH::Tensor_<DataType>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in, const std::vector<SMITH::IndexRange>& ranges );
     template<> 
     std::shared_ptr<SMITH::Tensor_<double>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<double>> Tens_in,  const std::vector<SMITH::IndexRange>& ranges );
     template<> 
     std::shared_ptr<SMITH::Tensor_<std::complex<double>>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens_in, const std::vector<SMITH::IndexRange>& ranges );

     template<typename DataType> 
     std::shared_ptr<SMITH::Tensor_<DataType>> get_sub_tensor_symm( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in, const std::vector<SMITH::IndexRange>& ranges,
                                                                    const std::vector<std::vector<int>>& transforms, const std::vector<DataType>& factors );
     template<> 
     std::shared_ptr<SMITH::Tensor_<double>> get_sub_tensor_symm( std::shared_ptr<SMITH::Tensor_<double>> Tens_in,  const std::vector<SMITH::IndexRange>& ranges,
                                                                  const std::vector<std::vector<int>>& transforms, const std::vector<double>& factors );
     template<> 
     std::shared_ptr<SMITH::Tensor_<std::complex<double>>> get_sub_tensor_symm( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens_in, const std::vector<SMITH::IndexRange>& ranges,
                                                                                const std::vector<std::vector<int>>& transforms, const std::vector<std::complex<double>>& factors );


     template<typename DataType> 
     std::shared_ptr<SMITH::Tensor_<DataType>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::vector<std::string>& range_names,
                                                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange> >> range_conversion_map );
     template<>
     std::shared_ptr<SMITH::Tensor_<double>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<double>> Tens_in,  std::vector<std::string>& range_names,
                                                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange> >> range_conversion_map );
     template<>
     std::shared_ptr<SMITH::Tensor_<std::complex<double>>> get_sub_tensor( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens_in,  std::vector<std::string>& range_names,
                                                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange> >> range_conversion_map );


     template<typename DataType>     
     void print_tensor_with_indexes( std::shared_ptr<SMITH::Tensor_<DataType>> tens, std::string name, bool print_zero_blocks = false  );
     template<>
     void print_tensor_with_indexes( std::shared_ptr<SMITH::Tensor_<double>> tens, std::string name, bool print_zero_blocks  );
     template<>
     void print_tensor_with_indexes( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> tens, std::string name, bool print_zero_blocks );

     template<typename DataType>     
     void set_test_elems( std::shared_ptr<SMITH::Tensor_<DataType>> Tens, std::string name  );
     template<>
     void set_test_elems( std::shared_ptr<SMITH::Tensor_<double>> Tens, std::string name  );
     template<>
     void set_test_elems( std::shared_ptr<SMITH::Tensor_<std::complex<double>>> Tens, std::string name  );

     template<typename DataType>     
     bool check_idx_rng_existence( const std::shared_ptr<SMITH::Tensor_<DataType>>& Tens, const std::vector<SMITH::IndexRange>& id_ranges  );
     template<>     
     bool check_idx_rng_existence( const std::shared_ptr<SMITH::Tensor_<double>>& Tens, const std::vector<SMITH::IndexRange>& id_ranges  );
     template<>     
     bool check_idx_rng_existence( const std::shared_ptr<SMITH::Tensor_<std::complex<double>>>& Tens, const std::vector<SMITH::IndexRange>& id_ranges  );
};
};
#endif 
