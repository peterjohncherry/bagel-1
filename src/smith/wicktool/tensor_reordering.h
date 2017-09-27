#ifndef __SRC_SMITH_WICKTOOL_TENSOR_SORTER_H
#define __SRC_SMITH_WICKTOOL_TENSOR_SORTER_H

#include <src/smith/wicktool/equation_tools.h>
namespace bagel {
namespace SMITH { 

namespace Tensor_Sorter { 

template<class DataType>
class Tensor_Sorter { 

    public: 
      Tensor_Sorter(){};
      ~Tensor_Sorter(){};
  
     std::unique_ptr<DataType[]> reorder_tensor_data(const DataType* orig_data,  size_t data_size, std::vector<int>  new_order_vec, std::vector<size_t> new_sizes_vec );

};
}
}
}
#endif
