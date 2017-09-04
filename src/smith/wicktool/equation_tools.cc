#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/CtrTensOp.h>
#include <src/smith/indexrange.h>
#include <src/smith/wicktool/equation_tools.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace equation_tools;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> inverse_reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
    newvec->at(neworder->at(ii)) =  origvec->at(ii);

  return newvec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
     newvec->at(ii) = origvec->at(neworder->at(ii));

  return newvec;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Contracts tensors T1 and T2 over two specified indexes
//T1_org_rg T2_org_rg are the original index ranges for the tensors (not necessarily normal ordered).
//T2_new_rg T1_new_rg are the new ranges, with the contracted index at the end, and the rest in normal ordering.
//T1_new_order and T2_new_order are the new order of indexes, and are used for rearranging the tensor data.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<DType> contract_different_tensors( string T1name, string T2name, pair<int,int> ctr_todo,
                                              shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto T1 = Tmap->at(T1name);
  auto T2 = Tmap->at(T2name);
 
  auto T1_org_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);//SHOULD BE TAKEN FROM CTP 
  auto T2_org_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);//SHOULD BE TAKEN FROM CTP

  auto T_out = make_shared<DType>(); 

  auto T1_new_order = make_shared<vector<int>>(0);
  for (int ii = 0; ii !=T1->unc_pos->size(); ii++)
    if (T1->unc_pos->at(ii) != ctr_todo.first)
      T1_new_order->push_back(T1->unc_pos->at(ii));
    
  T1_new_order->push_back(ctr_todo.first);
  auto T1_new_rngs = reorder_vector(T1_new_order, T1_org_rngs);
  auto maxs1 = get_sizes(T1_new_rngs);
  int T1_num_total_blocks = accumulate(maxs1->begin(), maxs1->end(), 1, multiplies<int>());
 
  auto T2_new_order = make_shared<vector<int>>(1,ctr_todo.second);
  for (int ii = 0; ii !=T2->unc_pos->size(); ii++)
    if (T2->unc_pos->at(ii) != ctr_todo.second)
      T2_new_order->push_back(T2->unc_pos->at(ii));

  T2_new_order->push_back(ctr_todo.second);
  auto T2_new_rngs = reorder_vector(T2_new_order, T2_org_rngs);
  auto maxs2 = get_sizes(T2_new_rngs);
  maxs2->pop_back();
  int T2_num_unc_blocks = accumulate(maxs1->begin(), maxs1->end(), 1, multiplies<int>());
  

  auto Tout_unc_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T2_new_rngs->begin()+1, T2_new_rngs->end());
 

  //loops over all index blocks of T1 and  T2; inner loop for T2 has same final index as T2 due to contraction
  auto T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);
  for (int ii = 0 ; ii != T1_num_total_blocks; ii++) { 

    std::unique_ptr<double[]> T1_data_new; 
    size_t ctr_block_size;
    size_t T1_unc_block_size;
    auto T_out_rng_block = make_shared<vector<Index>>();

    {
    auto T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, T1_new_rngs); 
    auto T1_org_rng_blocks = inverse_reorder_vector(T1_new_order, T1_org_rng_blocks); 

    auto T1_new_rng_block_sizes = get_sizes(T1_new_rng_blocks);
    T_out_rng_block->insert(T_out_rng_block->end(), T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
    
    ctr_block_size = T1_new_rng_blocks->back().size(); 
    T1_unc_block_size = get_block_size( T1_new_rng_blocks, 0, T1_new_rng_blocks->size()); 
    
    auto T1_data_org = T1->get_block(*T1_org_rng_blocks);
    T1_data_new = reorder_tensor_data( T1_data_org.get(), get_block_size(T1_new_rng_blocks, 0, T1_new_rng_blocks->size()),  *T1_new_order, *T1_new_rng_block_sizes); 
    }


    auto T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size()-1, 0);
    for (int jj = 0 ; jj != T2_num_unc_blocks; jj++) { 
      
      std::unique_ptr<double[]> T2_data_new; 
      size_t T2_unc_block_size;
      T2_rng_block_pos->push_back(T1_rng_block_pos->back());

      {
      auto T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, T2_new_rngs); 
      auto T2_org_rng_blocks = inverse_reorder_vector(T2_new_order, T2_org_rngs); 

      auto T2_new_rng_block_sizes = get_sizes(T2_new_rng_blocks);
      T_out_rng_block->insert(T_out_rng_block->end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      
      T2_unc_block_size = get_block_size( T2_new_rng_blocks, 0, T2_new_rng_blocks->size()); 
      
      auto T2_data_org = T2->get_block(*T2_org_rng_blocks);
      T2_data_new = reorder_tensor_data(T2_data_org.get(), get_block_size(T2_new_rng_blocks, 0, T2_new_rng_blocks->size()), *T2_new_order, *T2_new_rng_block_sizes); 
      }
      
      std::unique_ptr<double[]> T_out_data(new double[T1_unc_block_size*T2_unc_block_size]);
     
      //should not use transpose; instead build T2_new_order backwards... 
      dgemm_("N", "T", T1_unc_block_size, ctr_block_size,  T1_unc_block_size, 1.0, T1_data_new.get(), T1_unc_block_size,
              T2_data_new.get(), ctr_block_size, 1.0, T_out_data.get(), T1_unc_block_size);
      

      T_out->put_block( T_out_data, *T_out_rng_block );
      T2_rng_block_pos->pop_back(); //remove last index; contracted index is cycled in T1 loop
      fvec_cycle(T2_rng_block_pos, maxs2 );
    }
    fvec_cycle(T1_rng_block_pos, maxs1 );
  }                                                                                                                                              
  return T_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets a block of data from the original indexes, and reorders it so the indexes are in T_new_rng order
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType, class DType>
unique_ptr<DataType[]>
get_reordered_Tensor_data(shared_ptr<vector<int>> rng_block_pos, shared_ptr<vector<const IndexRange>> T_org_rng,
                                                                 shared_ptr<vector<const IndexRange>> T_new_rng, shared_ptr<DType> Tens )  { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// auto blocksize = [](shared_ptr<vector<size_t>> rb_sizes, int startpos, int endpos ){
//      size_t tmp_sz = 0;
//      for ( int ii = startpos ; ii!=endpos+1; ii++)
//         tmp_sz += *rb_sizes[ii];
//      return tmp_sz;
//  };
//
// shared_ptr<vector<Index>> T_new_rng_blocks = get_rng_blocks( rng_block_pos, T_new_rng); 
//
// shared_ptr<vector<size_t>> T_new_rng_block_sizes = get_sizes(T_new_rng_blocks);
//                                                                                                                          
// shared_ptr<vector<Index>> T_org_rng_blocks = get_rng_blocks( fvec, T_org_rng); 
//
// std::unique_ptr<DataType[]> T_data_org = Tens->get_block(*T_org_idxs);
// 
// return  reorder_tensor_data(T_data_org.get(),  get_block_size(T_new_rng_blocks, T_new_rng_blocks->size() -1)  , *T_new_order, *T_idx_sizes); 
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>> get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));}
  return new_ids;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> get_sizes(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<int>>(); 
  for( auto elem : *rngvec ) 
     size_vec->push_back(elem->size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>> get_sizes(shared_ptr<vector<Index>> Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( auto elem : *Idvec ) 
     size_vec->push_back(elem.size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
//already a function for this.... so replace
////////////////////////////////////////////////////////////////////////////////////////
size_t get_block_size(shared_ptr<vector<Index>> Idvec, int startpos, int endpos) {
////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
//  for( int ii = startpos ; ii!=(endpos+1); ii++ ) 
//   block_size *= Idvec->at(ii).size();
  return  block_size;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This would be best put in another part, probably CTP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
reorder_tensor_data(const DataType* orig_data,  size_t data_size, vector<int>  new_order_vec, vector<size_t> new_sizes_vec ) {
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void build_range_conversion_map(std::shared_ptr<std::vector<pair<std::string, std::shared_ptr<IndexRange>>>> range_conversion_pairs ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  range_conversion_map = make_shared<map< string, shared_ptr<IndexRange>>>();
  for (auto rng_pair : *range_conversion_pairs)
    range_conversion_map->emplace(rng_pair.first, rng_pair.second);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> convert_str_to_Bagel_Index(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}

////////////////////////////////////////////////////
//template class CtrTensorPart<std::vector<double>>;
//template class CtrTensorPart<Tensor_<double>>;
///////////////////////////////////////////////////
#endif
