#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/indexrange.h>
#include <src/smith/wicktool/equation_tools.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This would be best put in another part, probably CTP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
 Equation_Computer::Equation_Computer::reorder_tensor_data(const DataType* orig_data,  size_t data_size, vector<int>  new_order_vec, vector<size_t> new_sizes_vec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const double fac1= 1.0;
  const double fac2= 1.0;

  unique_ptr<DataType[]> reordered_data(new DataType[data_size]);
  // Converts to arrays; sort_indices fixed size (defined in template) arrays, so I now have this ugly mess,
  if (new_order_vec.size() == 2 ){
    array<int,2> new_order_arr; 
    array<int,2> new_sizes_arr; 
    for (auto ii = 0; ii != new_order_vec.size(); ii++)
      new_order_arr[ii] = new_order_vec[ii];
    for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
      new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;

  }
  else if (new_order_vec.size() == 3 ){
      array<int,3> new_order_arr;  
      array<int,3> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
        new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(), new_order_arr) ;
 
  } else if (new_order_vec.size() == 4 ){
      array<int,4> new_order_arr;  
      array<int,4> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
        new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;
 
  } else if (new_order_vec.size() == 5 ){
      array<int,5> new_order_arr;  
      array<int,5> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
        new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;
 
  } else if (new_order_vec.size() == 6 ){
      array<int,6> new_order_arr;  
      array<int,6> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
        new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;
 
  } else if (new_order_vec.size() == 7 ){
      array<int,7> new_order_arr;  
      array<int,7> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
        new_sizes_arr[ii] = new_sizes_vec[ii];
    sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;
 
  } else if (new_order_vec.size() == 8 ){
      array<int,8> new_order_arr;   
      array<int,8> new_sizes_arr; 
      for (auto ii = 0; ii != new_order_vec.size(); ii++)
        new_order_arr[ii] = new_order_vec[ii];
      for (auto ii = 0; ii != new_sizes_vec.size(); ii++)
       new_sizes_arr[ii] = new_sizes_vec[ii];
   sort_indices( new_sizes_arr, fac1, fac2, orig_data, reordered_data.get(),  new_order_arr) ;
  }
 
  return reordered_data;

}

#endif
