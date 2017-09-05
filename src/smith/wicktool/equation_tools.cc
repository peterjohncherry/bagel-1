#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/indexrange.h>
#include <src/smith/wicktool/equation_tools.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Equation_Computer::Equation_Computer::Equation_Computer(std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<Tensor_<double>>> eqn_info_in,
                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map_in ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  eqn_info =  eqn_info_in;

  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  cc_ = ref->ciwfn()->civectors();
  det_ = ref->ciwfn()->civectors()->det();

  range_conversion_map = make_shared<map<string, shared_ptr<IndexRange>>>();
  
  const int max = ref->maxtile();
  auto closed_rng  =  make_shared<IndexRange>(IndexRange(ref->nclosed()-ref->ncore(), max, 0, ref->ncore()));
  auto active_rng  =  make_shared<IndexRange>(IndexRange(ref->nact(), min(10,max), closed_rng->nblock(), ref->ncore()+closed_rng->size()));
  auto virtual_rng =  make_shared<IndexRange>(IndexRange(ref->nvirt(), max, closed_rng->nblock()+active_rng->nblock(), ref->ncore()+closed_rng->size()+active_rng->size()));

  range_conversion_map->emplace("cor", closed_rng);//change the naming of the ranges from cor to clo... 
  range_conversion_map->emplace("act", active_rng);
  range_conversion_map->emplace("vir", virtual_rng);

  CTP_map = eqn_info_in->CTP_map;
  CTP_data_map = CTP_data_map_in;

}  


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_block_Tensor(string Tname){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   auto unc_ranges = make_shared<vector<string>>(0);
   for (auto pos : *(CTP_map->at(Tname)->unc_pos))
     unc_ranges->push_back( CTP_map->at(Tname)->id_ranges->at(pos));  

   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   cout << " A" << endl;      
   auto Bagel_id_ranges_csptr = make_shared<vector<shared_ptr<const IndexRange>>>(0) ; 
   auto range_lengths = make_shared<vector<int>>(0);
   for (auto idrng : *Bagel_id_ranges){
      Bagel_id_ranges_csptr->push_back(make_shared<const IndexRange>(idrng));
      range_lengths->push_back(idrng.size());
   }

   cout << " B" << endl;      
   auto block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   auto mins = make_shared<vector<int>>(unc_ranges->size(),0);  
    
   cout << " C" << endl;      
   auto block_tensor = make_shared<Tensor_<double>>(*Bagel_id_ranges);
   do {
     
     cout << "fvec = " ;cout.flush();  for (auto elem : *block_pos) { cout <<  elem <<  " "  ; } cout << endl;
     auto T_id_blocks = get_rng_blocks( block_pos, Bagel_id_ranges_csptr ); 
     cout << " D" << endl;      
     auto T_block_data = CTP_data_map->at(Tname)->get_block(*T_id_blocks);
     cout << " E" << endl;      
     block_tensor->put_block(T_block_data, *T_id_blocks);
     cout << " F" << endl;      

   } while (fvec_cycle(block_pos, range_lengths, mins ));
 
   //unique_ptr<double[]> Tblock_data = CTP_data_map->at(to_string(Tname[0]))->get_block(Bagel_id_ranges);

   return block_tensor;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_same_tensors( pair<int,int> ctr_todo, std::string Tname) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  auto CTP = CTP_map->at(Tname);
  auto CTP_data = CTP_data_map->at(Tname);
  shared_ptr<Tensor_<double>>  T_out;  

  auto T1_new_order = make_shared<vector<int>>(0);
  for (int ii = 0; ii !=CTP->unc_pos->size(); ii++)
    if (CTP->unc_pos->at(ii) != ctr_todo.first)
      T1_new_order->push_back(CTP->unc_pos->at(ii));
    
  T1_new_order->push_back(ctr_todo.first);
  T1_new_order->push_back(ctr_todo.second);

  auto T1_org_rngs = Get_Bagel_const_IndexRanges(CTP->id_ranges);

  auto T1_new_rngs = reorder_vector(T1_new_order, T1_org_rngs);
  auto maxs1 = get_sizes(T1_new_rngs);
  int T1_num_total_blocks = accumulate(maxs1->begin(), maxs1->end(), 1, multiplies<int>());

  auto Tout_unc_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T1_new_rngs->begin(), T1_new_rngs->end()-2);

  {
  auto Tout_unc_rngs_raw = make_shared<vector<IndexRange>>(0); 
  for (auto id_rng : *Tout_unc_rngs)
    Tout_unc_rngs_raw->push_back(*id_rng);
  
  T_out = make_shared<Tensor_<double>>(*Tout_unc_rngs_raw); 
  }

  //loops over all index blocks of T1 and  T2; inner loop for T2 has same final index as T2 due to contraction
  auto T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);
  for (int ii = 0 ; ii != T1_num_total_blocks; ii++) { 

    std::unique_ptr<double[]> T1_data_new; 
    size_t ctr_block_size;
    size_t T1_unc_block_size;
    auto T_out_rng_block = make_shared<vector<Index>>();

    {
    auto T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, T1_new_rngs); 
    auto T1_org_rng_blocks = inverse_reorder_vector(T1_new_order, T1_new_rng_blocks); 

    auto T1_new_rng_block_sizes = get_sizes(T1_new_rng_blocks);
    T_out_rng_block->insert(T_out_rng_block->end(), T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-2);
    
    ctr_block_size = T1_new_rng_blocks->back().size(); 
    T1_unc_block_size = get_block_size( T1_new_rng_blocks, 0, T1_new_rng_blocks->size()); 
    
    auto T1_data_org = CTP_data->get_block(*T1_org_rng_blocks);
    T1_data_new = reorder_tensor_data( T1_data_org.get(), get_block_size(T1_new_rng_blocks, 0, T1_new_rng_blocks->size()),  *T1_new_order, *T1_new_rng_block_sizes); 

    const int one = 1;
    const double one_d = 1.0;
    const int T1_ubs_int = T1_unc_block_size;
    std::unique_ptr<double[]> T_out_data(new double[T1_unc_block_size]);

    for (int kk = 0; kk != ctr_block_size; kk++ ) 
      cout <<" fix lapack problem" << endl;
//       daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
       //daxpy_(&T1_ubs_int, &one_d, T1_data_new.get(), &one, T_out_data.get()+(kk*T1_unc_block_size), 1.0, &one);
//void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
    
    T_out->put_block( T_out_data, *T1_new_rng_blocks  );       
    }

    fvec_cycle(T1_rng_block_pos, maxs1 );
  }                                                                                                                                              
  return T_out;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Contracts tensors T1 and T2 over two specified indexes
//T1_org_rg T2_org_rg are the original index ranges for the tensors (not necessarily normal ordered).
//T2_new_rg T1_new_rg are the new ranges, with the contracted index at the end, and the rest in normal ordering.
//T1_new_order and T2_new_order are the new order of indexes, and are used for rearranging the tensor data.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_different_tensors( pair<int,int> ctr_todo, std::string T1name, std::string T2name){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  auto CTP1 = CTP_map->at(T1name);
  auto CTP2 = CTP_map->at(T2name);
  std::shared_ptr<Tensor_<double>> CTP1_data = CTP_data_map->at(T1name);
  std::shared_ptr<Tensor_<double>> CTP2_data = CTP_data_map->at(T2name);

  auto T1_org_rngs = Get_Bagel_const_IndexRanges(CTP1->id_ranges) ;
  auto T2_org_rngs = Get_Bagel_const_IndexRanges(CTP2->id_ranges) ;

  shared_ptr<Tensor_<double>> T_out;  

  auto T1_new_order = make_shared<vector<int>>(0);
  for (int ii = 0; ii !=CTP1->unc_pos->size(); ii++)
    if (CTP1->unc_pos->at(ii) != ctr_todo.first)
      T1_new_order->push_back(CTP1->unc_pos->at(ii));
    
  T1_new_order->push_back(ctr_todo.first);
  auto T1_new_rngs = reorder_vector(T1_new_order, T1_org_rngs);
  auto maxs1 = get_sizes(T1_new_rngs);
  int T1_num_total_blocks = accumulate(maxs1->begin(), maxs1->end(), 1, multiplies<int>());
 
  auto T2_new_order = make_shared<vector<int>>(1,ctr_todo.second);
  for (int ii = 0; ii !=CTP2->unc_pos->size(); ii++)
    if (CTP2->unc_pos->at(ii) != ctr_todo.second)
      T2_new_order->push_back(CTP2->unc_pos->at(ii));

  T2_new_order->push_back(ctr_todo.second);
  auto T2_new_rngs = reorder_vector(T2_new_order, T2_org_rngs);
  auto maxs2 = get_sizes(T2_new_rngs);
  maxs2->pop_back();
  int T2_num_unc_blocks = accumulate(maxs1->begin(), maxs1->end(), 1, multiplies<int>());

  auto Tout_unc_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T2_new_rngs->begin()+1, T2_new_rngs->end());

  {
  auto Tout_unc_rngs_raw = make_shared<vector<IndexRange>>(0); 
  for (auto id_rng : *Tout_unc_rngs)
    Tout_unc_rngs_raw->push_back(*id_rng);
  
  T_out = make_shared<Tensor_<double>>(*Tout_unc_rngs_raw); 
  }

  //loops over all index blocks of T1 and  T2; inner loop for T2 has same final index as T2 due to contraction
  auto T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);
  for (int ii = 0 ; ii != T1_num_total_blocks; ii++) { 

    std::unique_ptr<double[]> T1_data_new; 
    size_t ctr_block_size;
    size_t T1_unc_block_size;
    auto T_out_rng_block = make_shared<vector<Index>>();

    {
    auto T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, T1_new_rngs); 
    auto T1_org_rng_blocks = inverse_reorder_vector(T1_new_order, T1_new_rng_blocks); 

    auto T1_new_rng_block_sizes = get_sizes(T1_new_rng_blocks);
    T_out_rng_block->insert(T_out_rng_block->end(), T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
    
    ctr_block_size = T1_new_rng_blocks->back().size(); 
    T1_unc_block_size = get_block_size( T1_new_rng_blocks, 0, T1_new_rng_blocks->size()); 
    
    auto T1_data_org = CTP1_data->get_block(*T1_org_rng_blocks);
    T1_data_new = reorder_tensor_data( T1_data_org.get(), get_block_size(T1_new_rng_blocks, 0, T1_new_rng_blocks->size()),  *T1_new_order, *T1_new_rng_block_sizes); 
    }


    auto T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size()-1, 0);
    for (int jj = 0 ; jj != T2_num_unc_blocks; jj++) { 
      
      std::unique_ptr<double[]> T2_data_new; 
      size_t T2_unc_block_size;
      T2_rng_block_pos->push_back(T1_rng_block_pos->back());

      {
      auto T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, T2_new_rngs); 
      auto T2_org_rng_blocks = inverse_reorder_vector(T2_new_order, T2_new_rng_blocks); 

      auto T2_new_rng_block_sizes = get_sizes(T2_new_rng_blocks);
      T_out_rng_block->insert(T_out_rng_block->end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      
      T2_unc_block_size = get_block_size( T2_new_rng_blocks, 0, T2_new_rng_blocks->size()); 
      
      auto T2_data_org = CTP2_data->get_block(*T2_org_rng_blocks);
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
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>> Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));}
  return new_ids;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<int>>(); 
  for( auto elem : *rngvec ) 
     size_vec->push_back(elem->size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>> Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<Index>> Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( auto elem : *Idvec ) 
     size_vec->push_back(elem.size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
//already a function for this.... so replace
////////////////////////////////////////////////////////////////////////////////////////
size_t Equation_Computer::Equation_Computer::get_block_size(shared_ptr<vector<Index>> Idvec, int startpos, int endpos) {
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>> Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map->at(rng)));

  return ranges_Bagel;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> Equation_Computer::Equation_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::inverse_reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
    newvec->at(neworder->at(ii)) =  origvec->at(ii);

  return newvec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
     newvec->at(ii) = origvec->at(neworder->at(ii));

  return newvec;
}


////////////////////////////////////////////////////
//template class CtrTensorPart<std::vector<double>>;
//template class CtrTensorPart<Tensor_<double>>;
///////////////////////////////////////////////////
#endif
