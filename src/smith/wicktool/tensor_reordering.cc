#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_tools.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This would be best put in another part, probably CTP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Tensor_Sorter::Tensor_Sorter::reorder_tensor_data(const DataType* orig_data,  size_t data_size, vector<int>  new_order_vec, vector<size_t> new_sizes_vec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const double fac1= 1.0;
  const double fac2= 1.0;

  int nbase = new_order_vec->size();
  int id_num = 0;
  for (int ii = 0; ii!=new_order_vec.size(); ii++)
    id_num+= pow(nbase, ii)*new_order_vec[ii]; 

  unique_ptr<DataType[]> reordered_data(new DataType[data_size]);
 
  return reordered_data;
}

#endif
