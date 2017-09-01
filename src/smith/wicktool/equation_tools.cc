#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/WickUtils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Contracts tensors T1 and T2 over two specified indexes
//T1_org_rg T2_org_rg are the original index ranges for the tensors (not necessarily normal ordered).
//T2_new_rg T1_new_rg are the new ranges, with the contracted index at the end, and the rest in normal ordering.
//T1_new_order and T2_new_order are the new order of indexes, and are used for rearranging the tensor data.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<DType> CtrMultiTensorPart<DType>::contract_different_tensors( string T1name, string T2name,  pair<int,int> ctr_todo,
                                                                         shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  auto toi = [](char cc ){ return((int)(cc) - (int)('0')); } ; //there must be a better way.....
  
  auto CTP_new->data = make_shared<DType>(); 

  //Binary_contract_arbitrary(T1_org_rgs, T2_org_rgs, T1_new_pos, T2_new_pos, CTP1->data, CTP2->data, CTP_new->data);
  std::complex<double> cx_one = (1.0,0.0);
  std::complex<double> cx_zero = (0.0,0.0);

  auto T1_bagel_rngs = make_shared<vector<const IndexRange>>(0);//MUST BE SET IN CTP WHEN INITIALIZED
  auto T2_bagel_rngs = make_shared<vector<const IndexRange>>(0);//MUST BE SET IN CTP WHEN INITIALIZED

  cout << "generating new order for T1 indexes with contracted indexes at the end " << endl;
  auto T1_new_order = make_shared<vector<int>>(T1->unc_pos);
  auto T1_new_rng = make_shared<vector<shared_ptr<const IndexRange>>>();
  for (int ii = 0; ii !=unc_pos->size(); ii++){
    if (pos != ctr_todo.first){
        T1_new_order->push_back(pos);
        T1_new_rng->push_back(T1_bagel_rngs->at(pos)); 
    }
  }

  cout << "generating new order for T2 indexes with contracted indexes at the end " << endl;
  auto T2_new_order = make_shared<vector<int>>(T2->unc_pos);
  auto T2_new_rng = make_shared<vector<shared_ptr<const IndexRange>>>();
  for (int ii = 0; ii !=unc_pos->size(); ii++){
    if (pos != ctr_todo.first){
        T2_new_order->push_back(pos);
        T2_new_rng->push_back(T2_bagel_rngs->at(pos)); 
    }
  }

  cout << "defining flen and forvec for looping through blocks " << endl;
  int flen =1;
  auto maxs  = make_shared<vector<int>>();
  for (auto elem : *out_range) {
     maxs->push_back(elem->size());
     flen *= elem->size();
  }
  auto fvec = make_shared<vector<int>>(T1_new_order->size()+T2_new_order->size(),0);

  for (int ii = 0 ; ii != flen; ii++) { 
    for (int ll = 0 ; ll!= T1_new_order->back()->size(); ll++) { 

      cout << "fvec = [ "; for (auto elem : *fvec ) {cout << elem  << " " ;} cout << " ]" << endl; 

      std::unique_ptr<complex<double>[]> T1_data_new = get_reordered_Tensor_data(fvec, T1_org_rngs, T1_new_order, CTP1->data ); 
      std::unique_ptr<complex<double>[]> T2_data_new = get_reordered_Tensor_data(fvec, T2_org_rngs, T2_new_order, CTP1->data ); 

      cout << "creating T_out data "  << endl;
      std::unique_ptr<complex<double>[]> T_out_data(new complex<double>[T1_unc_block_size*T2_unc_block_size]);
                                                                                                                 
      cout << " calculating T_out data"  << endl;
      zgemv_("N", T1_unc_block_size*T2_unc_block_size,  ctr_block_size, cx_one,                 
                  T1_data_new.get(), T1_unc_block_size*T2_unc_block_size,  T2_data_new.get(),                  
                  1, cx_one, T_out_data.get(), 1);                                                                   

      cout << " getting T_out idxs"  << endl;
      shared_ptr<vector<Index>> T1_out_idxs = get_idxs( fvec, out_range); 

      cout << "putting T_out_data  in block"  << endl;
      T_out->put_block( T_out_data, *T1_out_idxs  );       
      
      cout << "cycling fvec" << endl;
      fvec_cycle(fvec, maxs );
    }
  }                                                                                                                                              
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets a block of data from the original indexes, and reorders it so the indexes are in T_new_rng order
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
unique_ptr<double[]>
CtrMultiTensOp<DType>::CtrMultiTensOp::get_reordered_Tensor_data(shared_ptr<vector<int>> rng_block_pos, shared_ptr<vector<const IndexRange>> T_org_rng,
                                                                 shared_ptr<vector<const IndexRange>> T_new_rng, shared_ptr<DType> Tens )  { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  auto blocksize = [](shared_ptr<vector<size_t>> rb_sizes, int startpos, int endpos ){
       size_t tmp_sz = 0;
       for ( int ii = startpos ; ii!=endpos+1; ii++)
          tmp_sz += *rb_sizes[ii];
       return tmp_sz;
   };

  shared_ptr<vector<Index>> T_new_rng_blocks = get_rng_blocks( rng_block_pos, T_new_rng); 

  shared_ptr<vector<size_t>> T_new_rng_block_sizes = get_sizes(T_new_rng_blocks);
                                                                                                                           
  shared_ptr<vector<Index>> T_org_rng_blocks = get_rng_blocks( fvec, T_org_rng); 

  std::unique_ptr<double[]> T_data_org = Tens->get_block(*T_org_idxs);
  
  return  reorder_tensor_data(T_data_org.get(),  get_block_size(T_new_rng_blocks, T_new_rng_blocks->size() -1)  , *T_new_order, *T_idx_sizes); 
}
////////////////////////////////////////////////////////////////////////////////////////
template <class DType>
shared_ptr<vector<Index>> CtrMultiTensOp<DType>::CtrMultiTensOp::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));}
  return new_ids;
}
////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<vector<size_t>> CtrMultiTensOp<DType>::CtrMultiTensOp::get_sizes(shared_ptr<vector<Index>> Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( auto elem : *Idvec ) 
     size_vec->push_back(elem.size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
//already a function for this.... so replace
////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
size_t CtrMultiTensOp<DType>::CtrMultiTensOp::get_block_size(shared_ptr<vector<Index>> Idvec, int startpos, int endpos) {
////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( int ii = startpos ; ii!=(endpos+1); ii++ ) 
    block_size *= Idvec->at(ii).size();
  return  block_size;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This would be best put in another part, probably CTP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
unique_ptr<complex<double>[]>
CtrMultiTensOp<DType>::CtrMultiTensOp::reorder_tensor_data(const complex<double>* orig_data,  size_t data_size, vector<int>  new_order_vec, vector<size_t> new_sizes_vec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const double fac1= 1.0;
  const double fac2= 1.0;

  unique_ptr<complex<double>[]> reordered_data(new complex<double>[data_size]);
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
vector<DataType> MSRelCASPT2::MSRelCASPT2::reorder_vector(vector<DataType> unordered_vec, vector<int> new_idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << " MSRelCASPT2::reorder" << endl;
  
  vector<DataType> ordered_vec;
  for (int elem : new_idxs) 
    ordered_vec.push_back(unordered_vec[elem]);
  return ordered_vec;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void build_index_conversion_map(std::shared_ptr<std::vector<pair<std::string, std::shared_ptr<IndexRange>>>> range_conversion_pairs ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (auto rng_pair : *range_conversion_pairs)
    map->emplace(rng_pair.first, rng_pair.second);
 
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> convert_str_to_Bagel_Index(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(ranges_str->size());
  for ( auto rng : *id_ranges_str) 
    *ranges_Bagel->push_back(idx_conversion_map->at(rng));

  return;
}
#endif
