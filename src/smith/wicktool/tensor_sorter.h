#ifndef __SRC_SMITH_WICKTOOL_TENSOR_SORTER_H
#define __SRC_SMITH_WICKTOOL_TENSOR_SORTER_H

#include <src/util/prim_op.h>
#include <iostream>
namespace bagel {
namespace SMITH { 

namespace Tensor_Sorter { 

template<class DataType>
class Tensor_Sorter { 

    public: 
      Tensor_Sorter(){};
      ~Tensor_Sorter(){};
  
     std::string get_comb_str(std::shared_ptr<std::vector<size_t>> invec);

     std::string get_comb_str(std::array<int,4>&  sort_options);

     size_t get_perm_idnum(std::shared_ptr<std::vector<size_t>> seq) ;
  
     size_t get_sort_options_idnum(std::array<int,4>& sort_options) ;

     void sort_indices_2( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_3( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_4( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_5( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_6( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_7( std::shared_ptr<std::vector<size_t>> new_id_order,
                   std::shared_ptr<std::vector<size_t>> rlen,
                   std::array<int,4>& sort_options, 
                   const DataType* in, DataType* out ) ;
    
     void sort_indices_8( std::shared_ptr<std::vector<size_t>> new_id_order,
                          std::shared_ptr<std::vector<size_t>> rlen,
                          std::array<int,4>& sort_options, 
                          const DataType* in, DataType* out ) ;
    
};
}
}
}
#endif
