#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_tools.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Equation_Computer::Equation_Computer::Equation_Computer(std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<double>> eqn_info_in,
                                                        std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map_in,
                                                        std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  eqn_info =  eqn_info_in;
  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  cc_ = ref->ciwfn()->civectors();
  det_ = ref->ciwfn()->civectors()->det();

  GammaMap = eqn_info->GammaMap;
  CTP_map = eqn_info_in->CTP_map;
  CTP_data_map = CTP_data_map_in;
  
  range_conversion_map = range_conversion_map_in;

}  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_block_Tensor(string Tname){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_block_Tensor" << Tname << endl;
  
   if(  CTP_data_map->find(Tname) != CTP_data_map->end())
     return CTP_data_map->at(Tname);

   shared_ptr<vector<string>>     unc_ranges = CTP_map->at(Tname)->id_ranges;  
   cout << "unc_ranges = [ " ; cout.flush();  for  ( auto rng :  *unc_ranges )cout << rng << " " ; cout << " ] " << endl;
   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   cout << "Bagel_id_ranges->size() = "<< Bagel_id_ranges->size() << endl;
   cout << "range_lengths = " ;
   auto  range_lengths  = make_shared<vector<int>>(0); 
   for (auto idrng : *Bagel_id_ranges ){
      range_lengths->push_back(idrng.range().size()-1); cout << idrng.range().size() << " " ;
   }

   shared_ptr<Tensor_<double>> fulltens = CTP_data_map->at(Tname.substr(0,1));
   shared_ptr<Tensor_<double>> block_tensor = make_shared<Tensor_<double>>(*Bagel_id_ranges);
   block_tensor->allocate();

   cout << "allocated tensor" << endl;
   auto block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   auto mins = make_shared<vector<int>>(unc_ranges->size(),0);  
   do {
     
     cout << "block_pos = " ;cout.flush();  for (auto elem : *block_pos) { cout <<  elem <<  " "  ; } cout << endl;

     vector<Index> T_id_blocks(Bagel_id_ranges->size());
     for( int ii = 0 ;  ii != T_id_blocks.size(); ii++){
       T_id_blocks[ii] =  Bagel_id_ranges->at(ii).range(block_pos->at(ii));
     }    
     cout << "T_id_blocks sizes : " ; for (Index id : T_id_blocks){ cout << id.size() << " " ;} cout << endl;
  
     unique_ptr<double[]> T_block_data = fulltens->get_block(T_id_blocks);
     block_tensor->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));
 
   return block_tensor;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_on_same_tensor( pair<int,int> ctr_todo, std::string Tname) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::contract_on_same_tensor" <<endl; 

   cout << "Tname = " << Tname << endl;

   shared_ptr<CtrTensorPart<double>> CTP_old = CTP_map->at(Tname);
   shared_ptr<Tensor_<double>> CTP_data_old = find_or_get_CTP_data(Tname);

   // get original uncontracted ranges and positions of Ctrs relative to the current tensor
   pair<int,int> rel_ctr_todo;
   vector<int> unc_pos_new(0);
   vector<IndexRange> unc_ranges_old = CTP_data_old->indexrange(); 
   vector<IndexRange> unc_ranges_new(unc_ranges_old.size()-2);  

   vector<IndexRange>::iterator urn_iter = unc_ranges_new.begin();
   for ( int ii = 0 ; ii != CTP_old->unc_pos->size() ; ii++ ){
     if ( (CTP_old->unc_pos->at(ii) != ctr_todo.first) && (CTP_old->unc_pos->at(ii) != ctr_todo.second) ) {
       unc_pos_new.push_back(ii);
       *urn_iter++ = unc_ranges_old[ii];
     } else if (CTP_old->unc_pos->at(ii) == ctr_todo.first) {
       rel_ctr_todo.first  = ii;
     } else {
       rel_ctr_todo.second = ii;
     }
   }

   cout << "unc_pos_new = [ " ; cout.flush();  for  ( auto rng :  unc_pos_new ) {cout << rng << " " ;} cout << " ] " << endl;
 
   shared_ptr<Tensor_<double>> CTP_data_new = make_shared<Tensor_<double>>(unc_ranges_new);
   CTP_data_new->allocate();

   cout << "allocated tensor" << endl;
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges_old.size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges_old.size(),0);  
   shared_ptr<vector<int>> maxs = make_shared<vector<int>>(unc_ranges_old.size(),0);  
   for ( int ii = 0 ; ii != unc_ranges_old.size() ; ii++ ) 
      maxs->at(ii) = unc_ranges_old[ii].range().size()-1;

   cout << "maxs = [ " ; cout.flush();  for  ( auto sz :  *maxs) {cout << sz << " " ;} cout << " ] " << endl;
  
   int num_ctr_blocks = unc_ranges_old[rel_ctr_todo.first].size();

   for (int ii = 0 ; ii != num_ctr_blocks ; ii++){ 

     mins->at(rel_ctr_todo.first) = ii  ;
     maxs->at(rel_ctr_todo.first) = ii  ;
     mins->at(rel_ctr_todo.second) = ii ; 
     maxs->at(rel_ctr_todo.second) = ii ; 

     do {
       
       cout << "block_pos = " ;cout.flush();  for (auto elem : *block_pos) { cout <<  elem <<  " "  ; } cout << endl;
       vector<Index> CTP_id_blocks = get_rng_blocks(block_pos, unc_ranges_old); 
       cout << "CTP_id_blocks sizes : " ; for (Index id : CTP_id_blocks){ cout << id.size() << " " ;} cout << endl;
       unique_ptr<double[]> T_block_data = CTP_data_old->get_block(CTP_id_blocks);
     
       cout <<" fix lapack problem" << endl;
       // daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
       // daxpy_(&T1_ubs_int, &one_d, T1_data_new.get(), &one, T_out_data.get()+(kk*T1_unc_block_size), 1.0, &one);
       // void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
     } while (fvec_cycle_skipper(block_pos, maxs, mins ));

   
  } 
  return CTP_data_new;
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
cout << "Equation_Computer::contract_on_different_tensor" <<endl; 

  shared_ptr<CtrTensorPart<double>> CTP1 = CTP_map->at(T1name);  cout << "got info for " <<  T1name << endl;
  shared_ptr<CtrTensorPart<double>> CTP2 = CTP_map->at(T2name);  cout << "got info for " << T2name << endl; 

  // shift contraction indexes.
  pair<int,int> ctr_todo_rel = relativize_ctr_positions( ctr_todo, CTP1, CTP2 );

  cout << "orig_ctr = ( " << ctr_todo.first << ","<<ctr_todo.second << ")" << endl;
  cout << "rel  ctr = ( " << ctr_todo_rel.first << ","<<ctr_todo_rel.second << ")" << endl;

  shared_ptr<Tensor_<double>> CTP1_data = find_or_get_CTP_data(T1name);
  shared_ptr<vector<int>> T1_new_order  = put_ctr_at_back( CTP1->unc_pos, ctr_todo_rel.first);
  shared_ptr<vector<shared_ptr<const IndexRange>>> T1_org_rngs = Get_Bagel_const_IndexRanges(CTP1->full_id_ranges, CTP1->unc_pos) ;
  cout << "X1" << endl; 
  shared_ptr<vector<shared_ptr<const IndexRange>>> T1_new_rngs = reorder_vector(T1_new_order, T1_org_rngs);
  cout << "X2" << endl; 
  shared_ptr<vector<int>> maxs1 = get_num_index_blocks_vec(T1_new_rngs) ;

  cout << "T1_new_order = [" ; for (auto elem : *T1_new_order) { cout << elem << " " ; } cout << "]"<< endl;
  shared_ptr<Tensor_<double>> CTP2_data = find_or_get_CTP_data(T2name);
  shared_ptr<vector<int>> T2_new_order  = put_ctr_at_front( CTP2->unc_pos, ctr_todo_rel.second);
  cout << "X3" << endl; 
  shared_ptr<vector<shared_ptr<const IndexRange>>> T2_org_rngs = Get_Bagel_const_IndexRanges(CTP2->full_id_ranges, CTP2->unc_pos) ;
  cout << "X4" << endl; 
  shared_ptr<vector<shared_ptr<const IndexRange>>> T2_new_rngs = reorder_vector(T2_new_order, T2_org_rngs);
  cout << "X5" << endl; 
  shared_ptr<vector<int>> maxs2 = get_num_index_blocks_vec(T2_new_rngs) ;
  cout << "X6" << endl; 

  auto Tout_unc_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T2_new_rngs->begin()+1, T2_new_rngs->end());

  shared_ptr<Tensor_<double>> T_out;  

  {
  auto Tout_unc_rngs_raw = make_shared<vector<IndexRange>>(0); 
  for (auto id_rng : *Tout_unc_rngs)
    Tout_unc_rngs_raw->push_back(*id_rng);
  
  T_out = make_shared<Tensor_<double>>(*Tout_unc_rngs_raw); 
  }
  T_out->allocate();

  //loops over all index blocks of T1 and  T2; final index of T1 is same as first index of T2 due to contraction
  auto T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    cout << "T1_rng_block_pos = [ " ; cout.flush(); for (int pos : *T1_rng_block_pos ) {cout << pos << " " ; } cout << "] " <<  endl;
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector(T1_new_order, T1_new_rng_blocks); 
    
    size_t ctr_block_size = T1_new_rng_blocks->back().size(); 
    size_t T1_unc_block_size = get_block_size( T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1); 
    size_t T1_block_size =  get_block_size(T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 

    std::unique_ptr<double[]> T1_data_new(new double[ T1_block_size  ]);
    {
      std::unique_ptr<double[]> T1_data_org = CTP1_data->get_block(*T1_org_rng_blocks);
      T1_data_new = reorder_tensor_data_X(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }

    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    T2_rng_block_pos->front() = T1_rng_block_pos->back();

    do { 
      cout << "T2_rng_block_pos = [ " ; cout.flush(); for (int pos : *T2_rng_block_pos ) {cout << pos << " " ; } cout << "] " <<  endl;
      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, T2_new_rngs); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end() );
       
      std::unique_ptr<double[]> T2_data_new(  new double[ T2_unc_block_size*ctr_block_size]);
      {
        shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector(T2_new_order, T2_new_rng_blocks); 
        std::unique_ptr<double[]> T2_data_org = CTP2_data->get_block(*T2_org_rng_blocks);
        T2_data_new = reorder_tensor_data_X(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
      }
       
      std::unique_ptr<double[]> T_out_data(new double[T1_unc_block_size*T2_unc_block_size]);
     
      //should not use transpose; instead build T2_new_order backwards... 
      //dgemm_("N", "T", T1_unc_block_size, ctr_block_size,  T1_unc_block_size, 1.0, T1_data_new.get(), T2_unc_block_size,
      //        T2_data_new.get(), ctr_block_size, 1.0, T_out_data.get(), T1_unc_block_size);
   
      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      T_out->put_block( T_out_data, T_out_rng_block );

      //remove last index; contracted index is cycled in T1 loop
    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, 0 )) ;
  } while (fvec_cycle_test(T1_rng_block_pos, maxs1 ));
  return T_out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<int,int>
Equation_Computer::Equation_Computer::relativize_ctr_positions(pair <int,int> ctr_todo, 
                                                               shared_ptr<CtrTensorPart<double>>  CTP1,
                                                               shared_ptr<CtrTensorPart<double>>  CTP2){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Equation_Computer::Equation_Computer::find_or_get_CTP_data" << endl;
   pair<int,int> rel_ctr;

   int T1_orig_size = CTP1->full_id_ranges->size(); 
   int T2_orig_size = CTP2->full_id_ranges->size(); 

   if (ctr_todo.first >= T1_orig_size ){ 
     rel_ctr = make_pair(ctr_todo.second, ctr_todo.first-T2_orig_size);
   } else {
     rel_ctr = make_pair(ctr_todo.first, ctr_todo.second-T1_orig_size);
   }

  for ( int ii = 0 ; ii != CTP1->unc_pos->size() ; ii++ )
     if ( CTP1->unc_pos->at(ii) == rel_ctr.first){
       rel_ctr.first = ii;
       break;
     }

  for ( int ii = 0 ; ii != CTP2->unc_pos->size() ; ii++ )
     if (CTP2->unc_pos->at(ii) == rel_ctr.second){ 
       rel_ctr.second = ii;
       break;
     }

 return rel_ctr;


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::find_or_get_CTP_data(string CTP_name){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Equation_Computer::Equation_Computer::find_or_get_CTP_data" << endl;

  shared_ptr<Tensor_<double>> CTP_data;
  auto CTP_data_loc =  CTP_data_map->find(CTP_name); 
  if ( CTP_data_loc == CTP_data_map->end() ){
    CTP_data = get_block_Tensor(CTP_name); 
  } else {
    CTP_data = CTP_data_loc->second;
  }
  cout << "got data for " << CTP_name << endl; 
  
  return CTP_data;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector  to have ctr pos at back
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::put_ctr_at_back(shared_ptr<vector<int>> orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "put_ctr_at_back" << endl;
  vector<int> new_pos(orig_pos->size());
  
  vector<int>::iterator new_pos_it = new_pos.begin();
  for (int ii = 0; ii !=orig_pos->size(); ii++)
    if (orig_pos->at(ii) != ctr_pos ){
      *new_pos_it++ = orig_pos->at(ii);
    }
  *new_pos_it = ctr_pos;

  return make_shared<vector<int>>(new_pos);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::put_ctr_at_front(shared_ptr<vector<int>> orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "put_ctr_at_front  : " << "rel_ctr_pos = " << ctr_pos << "  " ; cout.flush();
 cout << "orig_pos = [ " ; for (int ii : *orig_pos ) {cout << ii << " " ;} cout << "]"<<endl;
  vector<int> new_pos(orig_pos->size());
  new_pos[0] = ctr_pos;
  vector<int>::reverse_iterator new_pos_it = new_pos.rbegin();
  for (int ii =orig_pos->size()-1; ii !=0 ;  ii--)
    if (orig_pos->at(ii) != ctr_pos ){
      *new_pos_it++ = orig_pos->at(ii);
    }
  return make_shared<vector<int>>(new_pos);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<Tensor_<double>>>>
Equation_Computer::Equation_Computer::get_gammas(int MM , int NN, string gamma_name){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  //Gets list of necessary gammas which are needed by truncating the gamma_range vector
  //e.g. [a,a,a,a,a,a] --> [a,a,a,a] --> [a,a]
  shared_ptr<vector<string>> gamma_ranges_str  = GammaMap->at(gamma_name)->id_ranges;
  auto gamma_ranges =  vector<shared_ptr<vector<IndexRange>>>(gamma_ranges_str->size()/2);
  for (int ii = 0 ; ii !=gamma_ranges.size(); ii++ ){ 
    auto gamma_ranges_str_tmp = make_shared<vector<string>>(gamma_ranges_str->begin(), gamma_ranges_str->end()-ii*2);
    gamma_ranges[ii] = Get_Bagel_IndexRanges( gamma_ranges_str_tmp); 
    cout << "gamma_ranges["<<ii<<"]  = [ "; cout.flush();
    for (auto elem : *gamma_ranges_str_tmp){  cout << elem << " " ;} cout<<  "]"<< endl;
  } 
  
  //Gamma data in vector format 
  shared_ptr<vector<shared_ptr<VectorB>>> gamma_data_vec = compute_gammas( MM, NN ) ;
  auto gamma_tensors = make_shared<vector<shared_ptr<Tensor_<double>>>>(0);
  for ( int ii = gamma_ranges.size()-1; ii != -1;  ii-- ) {
 
     cout << "get_gammas range_lengths = ";
     auto  range_lengths  = make_shared<vector<int>>(0); 
     for (auto idrng : *(gamma_ranges[ii]) ){
       range_lengths->push_back(idrng.range().size()-1); cout << idrng.range().size() << " " ;
     }

     Tensor_<double> new_gamma_tensor(*(gamma_ranges[ii]));
     new_gamma_tensor.allocate();
     
     auto block_pos = make_shared<vector<int>>(gamma_ranges[ii]->size(),0);  
     auto mins = make_shared<vector<int>>(gamma_ranges[ii]->size(),0);  

     do {
       
       cout << endl << "get_gammas block_pos = " ;cout.flush();  for (auto elem : *block_pos) { cout <<  elem <<  " "  ; } cout << endl;
       vector<Index> gamma_id_blocks(gamma_ranges[ii]->size());
       for( int jj = 0 ;  jj != gamma_id_blocks.size(); jj++)
         gamma_id_blocks[jj] =  gamma_ranges[ii]->at(jj).range(block_pos->at(jj));

       cout << endl << "get_gammas gamma_id_blocks sizes : " ; for (Index id : gamma_id_blocks){ cout << id.size() << " " ; cout.flush();}
       size_t gamma_block_size;
       size_t gamma_block_pos;
       tie(gamma_block_size, gamma_block_pos) = get_block_info( gamma_ranges[ii],  block_pos) ;

       cout << " gamma_block_size = " << gamma_block_size << endl;
       cout << " gamma_block_pos = "  << gamma_block_pos << endl;
       cout <<" gamma_data_vec->at("<<ii<<")->size() = "  << gamma_data_vec->at(ii)->size() << endl; 

       cout << endl << gamma_data_vec->size()  << gamma_data_vec->size() << endl; 
       unique_ptr<double[]> gamma_data_block(new double[gamma_block_size])  ;
//       copy_n(gamma_data_block.get(), gamma_block_size, gamma_data_vec->at(ii)->data());
       
       cout << " got gamma data successfully " << endl;
    //   new_gamma_tensor.put_block( gamma_data_block, gamma_id_blocks);
       cout << " put gamma data successfully " << endl;
     
     } while (fvec_cycle(block_pos, range_lengths, mins ));
     shared_ptr<Tensor_<double>> new_gamma_tensor_ptr = make_shared<Tensor_<double>>(new_gamma_tensor);
     gamma_tensors->push_back(new_gamma_tensor_ptr);
  }
  cout << "out of loop" << endl;
 
  return gamma_tensors;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple< size_t, size_t >
Equation_Computer::Equation_Computer::get_block_info(shared_ptr<vector<IndexRange>> id_ranges, 
                                                     shared_ptr<vector<int>> block_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::get_block_info" << endl;
  vector<size_t> id_pos(block_pos->size());
  cout <<endl << "id_pos = ";
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){

    const size_t range_size         = id_ranges->at(ii).size();
    const size_t biggest_block_size = id_ranges->at(ii).range(0).size();
    const size_t num_blocks         = id_ranges->at(ii).range().size();
    const size_t remainder          = num_blocks * biggest_block_size - range_size;

    if (block_pos->at(ii) <= remainder  ){
       id_pos[ii] = num_blocks*block_pos->at(ii);//  + id_ranges->at(ii).range(block_pos->at(ii)).offset();

    } else if ( block_pos->at(ii) > remainder ) {
       id_pos[ii] = num_blocks*(range_size - remainder)+(num_blocks-1)*(remainder - block_pos->at(ii));// + id_ranges->at(ii).range(block_pos->at(ii)).offset(); 
    }; 
    cout << id_pos[ii] << " " ;
  }

  cout << endl << "range_sizes = " ;
  // getting size of ranges (seems to be correctly offset for node)
  vector<size_t> range_sizes(block_pos->size());
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){
    range_sizes[ii]  = id_ranges->at(ii).size();
    cout << range_sizes[ii] << " " ;
  }

  size_t data_block_size = 1;
  size_t data_block_pos  = 0;
  for (int ii = 0 ; ii != id_ranges->size()-1 ; ii++){
    data_block_pos  += id_pos[ii]*(pow(range_sizes[ii] , id_ranges->size()-ii));
    data_block_size *= id_ranges->at(ii).range(block_pos->at(ii)).size();
  }

  data_block_size *= id_ranges->back().range(block_pos->back()).size();
 
  return tie(data_block_size, data_block_pos);
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
unique_ptr<double[]>
Equation_Computer::Equation_Computer::get_block_of_data( double* data_ptr ,
                                                         shared_ptr<vector<IndexRange>> id_ranges, 
                                                         shared_ptr<vector<int>> block_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::get_block_of_data" << endl;

 // merge this to one loop, but keep until debugged, as this is likely to go wrong...
 // getting id position of id_block ; list of block sizes looks like [n,n,n,n,... , n-1, n-1 n-1] 
  vector<size_t> id_pos(block_pos->size());
  cout <<endl << "id_pos = ";
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){

    const size_t range_size         = id_ranges->at(ii).size();
    const size_t biggest_block_size = id_ranges->at(ii).range(0).size();
    const size_t num_blocks         = id_ranges->at(ii).range().size();
    const size_t remainder          = num_blocks * biggest_block_size - range_size;

    if (block_pos->at(ii) <= remainder  ){
       id_pos[ii] = num_blocks*block_pos->at(ii);//  + id_ranges->at(ii).range(block_pos->at(ii)).offset();

    } else if ( block_pos->at(ii) > remainder ) {
       id_pos[ii] = num_blocks*(range_size - remainder)+(num_blocks-1)*(remainder - block_pos->at(ii));// + id_ranges->at(ii).range(block_pos->at(ii)).offset(); 
    }; 
    cout << id_pos[ii] << " " ;
  }

  cout << endl << "range_sizes = " ;
  // getting size of ranges (seems to be correctly offset for node)
  vector<size_t> range_sizes(block_pos->size());
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){
    range_sizes[ii]  = id_ranges->at(ii).size();
    cout << range_sizes[ii] << " " ;
  }

  size_t data_block_size = 1;
  size_t data_block_pos  = 1;
  for (int ii = 0 ; ii != id_ranges->size()-1 ; ii++){
    data_block_pos  *= id_pos[ii]*range_sizes[ii];
    data_block_size *= id_ranges->at(ii).range(block_pos->at(ii)).size();
  }

  data_block_size *= id_ranges->back().range(block_pos->back()).size();
 
  cout << "data_block_size = " << data_block_size << endl;
  cout << "data_block_pos = "  << data_block_pos << endl;
 
  unique_ptr<double[]> data_block(new double[data_block_size])  ;

  copy_n(data_block.get(), data_block_size, data_ptr+data_block_pos);

  return data_block; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma1  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<VectorB>>> 
Equation_Computer::Equation_Computer::compute_gammas(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12 MM = " << MM << " NN = " << NN  << endl;

  if (det_->compress()) { // uncompressing determinants
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }
  shared_ptr<Civec> ccbra = make_shared<Civec>(*cc_->data(MM));
  shared_ptr<Civec> ccket = make_shared<Civec>(*cc_->data(NN));

  shared_ptr<RDM<1>> gamma1x;
  shared_ptr<RDM<2>> gamma2x;
  shared_ptr<RDM<3>> gamma3x;
  tie(gamma1x, gamma2x, gamma3x) = compute_gamma12_from_civec(ccbra, ccket);
  
  auto  gamma1 = make_shared<VectorB>(norb_*norb_);
  copy_n(gamma1x->data(), norb_*norb_, gamma1->data());

  auto  gamma2 = make_shared<VectorB>(norb_*norb_*norb_*norb_);
  copy_n(gamma2x->data(), norb_*norb_*norb_*norb_, gamma2->data());

  auto  gamma3 = make_shared<VectorB>(norb_*norb_*norb_*norb_*norb_*norb_);
  copy_n(gamma3x->data(), norb_*norb_*norb_*norb_*norb_*norb_, gamma3->data());
  cc_->set_det(det_); 

  auto gamma_vec = make_shared<vector<shared_ptr<VectorB>>>(vector<shared_ptr<VectorB>> { gamma1, gamma2, gamma3}) ;
                    
 return gamma_vec;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma1  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12 MM = " << MM << " NN = " << NN  << endl;

  if (det_->compress()) { // uncompressing determinants
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }
  shared_ptr<Civec> ccbra = make_shared<Civec>(*cc_->data(MM));
  shared_ptr<Civec> ccket = make_shared<Civec>(*cc_->data(NN));
 
  shared_ptr<RDM<1>> gamma1;
  shared_ptr<RDM<2>> gamma2;
  shared_ptr<RDM<3>> gamma3;
  tie(gamma1, gamma2, gamma3) = compute_gamma12_from_civec(ccbra, ccket);
 
  cc_->set_det(det_); 

  return tie(gamma1, gamma2, gamma3);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "compute_gamma12_from_civec" << endl;

  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);
  
  shared_ptr<Dvec> dket;
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb_*norb_);
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }
  
  return compute_gamma12_last_step(dbra, dket, cbra);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12_last_step" << endl;

  const int nri = cibra->asize()*cibra->lenb();
  const int n1 = norb_;  
  const int n2 = n1*norb_;  
 
  // gamma1 c^dagger <I|\hat{E}|0>
  // gamma2 \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto gamma1 = make_shared<RDM<1>>(norb_);
  auto gamma2 = make_shared<RDM<2>>(norb_);
  auto gamma3 = make_shared<RDM<3>>(norb_);

  //section to be made recursive for arbitrary orders of gamma
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, n2);
    for (int i = 0; i != n2; ++i)
      copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
 
    auto gamma1t = btas::group(*gamma1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, gamma1t, {1});

    auto dbra_data = dket_data;
    if (dbra != dket) {
      dbra_data = make_shared<Matrix>(nri, n2);
      for (int i = 0; i != n2; ++i)
        copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
    }

    const char   transa = 'N';
    const char   transb = 'T';
    const double alpha = 1.0;
    const double beta = 0.0; 
    const int n3 = n2*norb_;  
    const int n4 = n2*n2;
 
    //Very bad way of getting gamma3  [sum_{K}<I|i*j|K>.[sum_{L}<K|k*l|L>.<L|m*n|J>]]
    for ( int mn = 0; mn!=n2 ; mn++){
      Matrix gamma3_klmn(n2, n2);
      int m = mn/n1;
      int n = mn-(m*n1);
      auto dket_klmn = make_shared<Dvec>(dket->det(), n2);
      sigma_2a1(dket->data(mn), dket_klmn);
      sigma_2a2(dket->data(mn), dket_klmn);

      //make gamma3_klmn block
      auto dket_data_klmn = make_shared<Matrix>(nri, n2);
      for (int kl = 0; kl != n2; kl++)
        copy_n(dbra->data(kl)->data(), nri, dket_data_klmn->element_ptr(0, kl));
      btas::contract(1.0, *dbra_data, {1,0}, *dket_data_klmn, {1,2}, 0.0, gamma3_klmn, {0,2});
       
      //copy gamma3_klmn block into gamma3 and undo transpose of dbra in btas::contract
      unique_ptr<double[]> buf(new double[norb_*norb_]);
      for (int k = 0; k != n1; ++k) {
        for (int l = 0; l != n1; ++l) {
          copy_n(gamma3_klmn.data(), n2, buf.get());
          blas::transpose(buf.get(), n1, n1, gamma3->element_ptr(0,0,k,l,m,n));
        }
      }
    }

    auto gamma2t = group(group(*gamma2, 2,4), 0,2);
    btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, gamma2t, {0,2});
  }
 
  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&gamma2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, gamma2->element_ptr(0,0,k,i));
    }
  }
 
  return tie(gamma1, gamma2, gamma3);
}
//////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
//////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a1(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma) const {
//////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "sigma_2a1" << endl;
  const int lb = sigma->lenb();
  const int ij = sigma->ij();
  const double* const source_base = cvec->data();

  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = sigma->data(ip)->data();

    for (auto& iter : cvec->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      blas::ax_plus_y_n(sign, source_base + iter.target*lb, lb, target_array);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
// I'm sure there was a version which used transposition of the civector; this looks slow.
///////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a2(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma) const {
///////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "sigma_2a2" << endl;
  const int la = sigma->lena();
  const int ij = sigma->ij();

  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cvec->element_ptr(0, i);

    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = sigma->data(ip)->element_ptr(0, i);

      for (auto& iter : cvec->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////
vector<Index> Equation_Computer::Equation_Computer::get_rng_blocks_raw(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr< IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
  vector<Index> new_ids(forvec->size()); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids[ii] = old_ids->at(ii)->range(forvec->at(ii));
  }

  return new_ids;
}
////////////////////////////////////////////////////////////////////////////////////////
vector<Index>
Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> forvec, vector<IndexRange>& id_ranges) {
////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::get_rng_blocks " << endl;

  vector<Index> new_ids(id_ranges.size()); 
  for( int ii =0; ii != forvec->size(); ii++ ) 
     new_ids[ii] = id_ranges[ii].range(forvec->at(ii));

  
  return new_ids;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>> Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::get_rng_blocks constver" << endl;

  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));
  }
  return new_ids;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>> Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::get_rng_blocks " << endl;

  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));
  }
  return new_ids;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec, int skip_id) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<int>>(); 
  for( int ii =0 ; ii!=rngvec->size();ii++)
    if(ii!=skip_id) 
     size_vec->push_back(rngvec->at(ii)->size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>> Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<Index>> Idvec, int skip_id) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( int ii =0 ; ii!=Idvec->size();ii++)
    if(ii!=skip_id) 
     size_vec->push_back(Idvec->at(ii).size());
  return size_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
vector<int>
Equation_Computer::Equation_Computer::get_num_index_blocks_vec(vector<IndexRange>& rngvec) {
////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec.size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( IndexRange id_rng : rngvec ) 
     *iter++ = id_rng.range().size()-1;

  return num_id_blocks_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>>
Equation_Computer::Equation_Computer::get_num_index_blocks_vec(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec->size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( shared_ptr<const IndexRange> id_rng : *rngvec ) 
     *iter++ = id_rng->range().size()-1;

  return make_shared<vector<int>>(num_id_blocks_vec);
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
shared_ptr<vector<size_t>> Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<Index>> Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( auto elem : *Idvec ) 
     size_vec->push_back(elem.size());
  return size_vec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t
Equation_Computer::Equation_Computer::get_unc_block_size( vector<Index>& idvec, pair<int,int> ctr ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  
  size_t block_size = 1; 
  for( int ii = 0 ; ii != idvec.size() ; ii++) 
    if ( (ii != ctr.first) && (ii !=ctr.second))
      block_size *= idvec[ii].size();
  
  return  block_size;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Equation_Computer::Equation_Computer::get_block_size( vector<Index>::iterator startpos,
                                                             vector<Index>::iterator endpos   ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( vector<Index>::iterator id_it = startpos ; id_it!=endpos; id_it++ ){ 
    block_size *= id_it->size();
  }
  return  block_size;
}

////////////////////////////////////////////////////////////////////////////////////////
//already a function for this.... so replace
////////////////////////////////////////////////////////////////////////////////////////
size_t Equation_Computer::Equation_Computer::get_block_size(shared_ptr<vector<Index>> Idvec, int startpos, int endpos) {
////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( int ii = startpos ; ii!=(endpos+1); ii++ ) 
    block_size *= Idvec->at(ii).size();
  return  block_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
 Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 1arg" << endl;

  auto ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map->at(rng)));

  return ranges_Bagel;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
 Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str, shared_ptr<vector<int>> unc_pos){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 2arg" << endl;

  vector<shared_ptr<const IndexRange>>  ranges_Bagel(unc_pos->size());
  for ( int ii = 0 ; ii != unc_pos->size() ; ii++) 
    ranges_Bagel[ii]=(make_shared<const IndexRange>(*range_conversion_map->at(ranges_str->at(unc_pos->at(ii)))));

  return make_shared<vector<shared_ptr<const IndexRange>>>(ranges_Bagel);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> Equation_Computer::Equation_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_IndexRanges 1arg" << endl;

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::inverse_reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "inverse_reorder_vector" << endl;

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
    newvec->at(neworder->at(ii)) =  origvec->at(ii);

  return newvec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "reorder_vector" << endl;

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
     newvec->at(ii) = origvec->at(neworder->at(ii));

  return newvec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//New version 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Equation_Computer::Equation_Computer::reorder_tensor_data_X(const DataType* orig_data, shared_ptr<vector<int>>  new_order_vec,
                                                            shared_ptr<vector<Index>> orig_index_blocks ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "reorder_tensor_data_X" << endl;

  shared_ptr<vector<size_t>> rlen = get_sizes(orig_index_blocks);
  shared_ptr<vector<size_t>> new_order_st = make_shared<vector<size_t>>(new_order_vec->size());   
  size_t block_size = get_block_size(orig_index_blocks->begin(), orig_index_blocks->end());
  array<int,4> sort_options = {0,1,1,1};

  unique_ptr<DataType[]> reordered_data(new DataType[block_size]);

   for ( int ii = 0 ; ii != new_order_vec->size(); ii++) 
     new_order_st->at(ii) = new_order_vec->at(ii);

  Tensor_Sorter::Tensor_Sorter<double> TS ;

  size_t num_ids =  orig_index_blocks->size();
  if ( num_ids == 2) { 
    TS.sort_indices_2( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 3 ) {
    TS.sort_indices_3( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 4 ) {
    TS.sort_indices_4( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 5 ) {
    TS.sort_indices_5( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 6 ) {
    TS.sort_indices_6( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 7 ) {
    TS.sort_indices_7( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  }

  return reordered_data;
}

#endif

 // int n1 = norb_;  
 // int n2 = norb_*norb_;  
 // int n3 = n2*norb_;  
 // int n4 = n3*norb_;  
 // int n5 = n4*norb_;  
 //                       
//  cout << endl<< "gamma3 " << endl;
//  for (int i = 0; i != norb_; ++i) 
//    for (int j = 0; j != norb_; ++j) 
//      for (int k = 0; k != norb_; ++k) 
//        for (int l = 0; l != norb_; ++l) 
//          for (int m = 0; m != norb_; ++m) {
//            cout << endl << i << " " << j << " "<<  k << " " << l << " " << m << "  ";
//            for (int n = 0; n != norb_; ++n){
//              cout << *(gamma3->data()+(i*n5 + j*n4 + k*n3 + l*n2 + m*n1 +n )) << " ";
//             }
//          }

 // cout <<endl<<"gamma2 " << endl;
 // for (int k = 0; k != norb_; ++k) 
 //   for (int l = 0; l != norb_; ++l) 
 //     for (int m = 0; m != norb_; ++m) {
 //       cout << endl << k << " " << l << " " << m << "  ";
 //       for (int n = 0; n != norb_; ++n){ 
 //         cout << *(gamma2->data()+ (k*n3 + l*n2 + m*n1 +n))  << " " ;
 ////       }                        
 ////     }
 ////
 //// cout << endl << "gamma1 " << endl;
 // for (int m = 0; m != norb_; ++m) {
 //   cout << endl << m << " " ;
 //   for (int n = 0; n != norb_; ++n){ 
 //     cout << *(gamma1->data()+ ( m*n1 +n )) << " "; 
 //   }
 // }
 //  cout << " maxs1 = [ ";   for (int max1 : *maxs1 ) { cout << max1 << " " ; } cout << "]" << endl;

 // for ( int jj = 0 ; jj != T1_new_rngs->size() ; jj++ ){ 
 //    cout << "T1_id_block_sizes->at("<<jj<<") = [" ; cout.flush();
 //   for ( int ii = 0 ; ii != T1_new_rngs->at(jj)->range().size(); ii++ ) {
 //     cout << T1_new_rngs->at(jj)->range(ii).size() << " " ;  cout.flush();
 //   }
 //   cout << "]" << endl;
 // }


