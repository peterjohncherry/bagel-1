#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/tensor_arithmetic.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;
using namespace bagel::SMITH::Tensor_Arithmetic_Utils; 
using namespace WickUtils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_on_same_tensor( shared_ptr<Tensor_<DataType>> Tens_in,  pair<int,int> ctr_todo) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::contract_on_same_tensor" << endl;

   // get original uncontracted ranges and positions of Ctrs relative to the current tensor
   int ctr1 = ctr_todo.first;
   int ctr2 = ctr_todo.second;
   vector<IndexRange> unc_ranges_old = Tens_in->indexrange(); 
   vector<IndexRange> unc_ranges_new(unc_ranges_old.size()-2);  
   vector<int> unc_pos_new(unc_ranges_old.size()-2);  
         
   vector<IndexRange>::iterator urn_iter = unc_ranges_new.begin();
   vector<int>::iterator upn_iter = unc_pos_new.begin();
   for ( int ii = 0 ; ii != unc_ranges_old.size() ; ii++ )
     if ( (ii != ctr_todo.first) && (ii != ctr_todo.second) ){ 
       *urn_iter++ = unc_ranges_old[ii];
       *upn_iter++ = ii;
     }

   shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(unc_ranges_new);
   Tens_out->allocate();
   int num_ctr_blocks = unc_ranges_old[ctr_todo.first].range().size();

   //loops over index blocks where ctr1 = ctr2 
   for (int ii = 0 ; ii != num_ctr_blocks ; ii++){ 
     shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> maxs = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     for ( int jj = 0 ; jj != unc_ranges_old.size() ; jj++ ) 
        maxs->at(jj) = unc_ranges_old[jj].range().size()-1;
 
     mins->at(ctr_todo.first)  = ii;
     mins->at(ctr_todo.second) = ii; 
     maxs->at(ctr_todo.first)  = ii;
     maxs->at(ctr_todo.second) = ii;
     block_pos->at(ctr_todo.first)  = ii;
     block_pos->at(ctr_todo.second) = ii;
     
     const int ione = 1; 
     const DataType done = 1.0; 
     do {
       
       shared_ptr<vector<pair<size_t, size_t>>> bob = get_block_start( make_shared<vector<IndexRange>>(unc_ranges_old), block_pos ) ;
       vector<Index> Tens_id_blocks_old = *(get_rng_blocks(block_pos, unc_ranges_old)); 
       vector<Index> Tens_id_blocks_new(Tens_id_blocks_old.size()-2);
       for (int kk = 0 ; kk != unc_pos_new.size(); kk++)       
         Tens_id_blocks_new[kk] = Tens_id_blocks_old[unc_pos_new[kk]];
     
       vector<int> range_sizes = get_sizes(Tens_id_blocks_old);
       int total_size = accumulate( range_sizes.begin(), range_sizes.end(), 1, std::multiplies<int>() );

       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(range_sizes);
       int inner_stride = Tens_strides->at(ctr1) < Tens_strides->at(ctr2) ? Tens_strides->at(ctr1) : Tens_strides->at(ctr2);

       shared_ptr<vector<int>> maxs2 = make_shared<vector<int>>(range_sizes.size(),0);
       shared_ptr<vector<int>> mins2 = make_shared<vector<int>>(maxs->size(), 0); 
       shared_ptr<vector<int>> fvec2 = make_shared<vector<int>>(*mins); 

       //within index block, loop through data copying chunks where ctr1 = ctr2 
       //Has odd striding; probably very inefficient when ctr1 or ctr2 are the leading index.
       for (int jj = ctr1+1 ; jj != maxs2->size(); jj++) 
         maxs2->at(jj) = range_sizes[jj]-1; // note, copying inner block as vector, so want max = min = 0 for those indexes;
       
       int ctr1_rlen = range_sizes[ctr1];          
       int tmp = total_size/(ctr1_rlen*ctr1_rlen); 
       unique_ptr<DataType[]>      Tens_data_block_old = Tens_in->get_block(Tens_id_blocks_old);
       std::unique_ptr<DataType[]> Tens_data_block_new(new DataType[tmp]);
       std::fill_n(Tens_data_block_new.get(), tmp, 0.0);
       shared_ptr<vector<int>>   CTens_strides   =  get_CTens_strides( range_sizes, ctr1, ctr2 );
       
       for ( int jj = 0 ; jj != ctr1_rlen ; jj++ ) {

         fill(fvec2->begin(), fvec2->end(), 0);
         maxs2->at(ctr1) = jj; 
         maxs2->at(ctr2) = jj; 
         mins2->at(ctr1) = jj;  
         mins2->at(ctr2) = jj; 
         fvec2->at(ctr1) = jj; 
         fvec2->at(ctr2) = jj; 

         do { 
           int pos   = inner_product( fvec2->begin(), fvec2->end(), CTens_strides->begin(), 0); 
           int inpos = inner_product( fvec2->begin(), fvec2->end(), Tens_strides->begin(),  0); 

           blas::ax_plus_y_n(1.0, Tens_data_block_old.get()+inpos, inner_stride, Tens_data_block_new.get()+pos);

         } while (fvec_cycle_skipper_f2b(fvec2, maxs2, mins2));
       }
 
       Tens_out->add_block( Tens_data_block_new, Tens_id_blocks_new );
     } while (fvec_cycle_skipper(block_pos, maxs, mins ));
  } 

  return Tens_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Will contract over all the indexes given in the conctracted indexes list. 
//Note these are indexes relative to the input tensor, not the original Tensor from which the CtrTensorPart was obtained.
//In principal this could replace the contract_same_tensor routine, as there's a lot of duplicated code
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_on_same_tensor( std::shared_ptr<Tensor_<DataType>> Tens_in, vector<int>& contracted_index_positions ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::get_trace_tensor" << endl;

   if ( contracted_index_positions.size()-  Tens_in->indexrange().size() )  {
     cout << "contracted_index_positions         = [ "; cout.flush();  for ( int pos : contracted_index_positions ) {cout << pos << " " ; } cout << " ] " << endl;
     cout << "contracted_index_positions.size()  = "<<  contracted_index_positions.size() << endl;
     cout << "Tens_in->indexrange().size()       = " << Tens_in->indexrange().size()  <<  endl;
     cout << "trying to contract more indexes than there are in the tensor!" << endl;
     assert(false);
   }

   int num_ctrs = contracted_index_positions.size();
   // get original uncontracted ranges and positions of Ctrs relative to the current tensor
   vector<IndexRange> unc_ranges_old = Tens_in->indexrange(); 
   vector<IndexRange> unc_ranges_new(unc_ranges_old.size() - num_ctrs);  
   vector<int> unc_pos_new(unc_ranges_old.size() - num_ctrs);  

   vector<bool> unc_get(Tens_in->indexrange().size(), true);
   for( int ii = 0; ii !=contracted_index_positions.size(); ii++ )
      unc_get[contracted_index_positions[ii]] = false;

   for ( int ii = 0 ; ii != unc_get.size(); ii++  )
     if (unc_get[ii]) {
       unc_ranges_new[ii] = unc_ranges_old[ii];
       unc_pos_new[ii] = ii;
     }

   shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(unc_ranges_new);
   Tens_out->allocate();
   Tens_out->zero();

   // loops over index blocks with same ctr
   int num_ctr_blocks = unc_ranges_old[contracted_index_positions[0]].range().size();
   for (int ii = 0 ; ii != num_ctr_blocks ; ii++){ 
     shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> maxs = make_shared<vector<int>>(unc_ranges_old.size());  
     for ( int jj = 0 ; jj != unc_ranges_old.size() ; jj++ ) 
        maxs->at(jj) = unc_ranges_old[jj].range().size()-1;
 
     for ( int ctr_idx_pos : contracted_index_positions ) {
        mins->at(ctr_idx_pos)       = ii;
        maxs->at(ctr_idx_pos)       = ii;
        block_pos->at(ctr_idx_pos)  = ii;
     }
     
     do {
       
       vector<Index> Tens_id_blocks_old = *(get_rng_blocks(block_pos, unc_ranges_old)); 
       vector<Index> Tens_id_blocks_new(Tens_id_blocks_old.size()- num_ctrs );
       for (int kk = 0 ; kk != unc_pos_new.size(); kk++)       
         Tens_id_blocks_new[kk] = Tens_id_blocks_old[unc_pos_new[kk]];
    
 
       vector<int> range_sizes = get_sizes(Tens_id_blocks_old);
       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(range_sizes);
   
       // To save time,  the contracted parts are copied in blocks i.e. (i, j, k, k, : , : ) 
       // Inner stride is the size of the data block to the right of the contracted indexes
       int inner_stride = Tens_strides->at(0); 
       int innermost_ctr = contracted_index_positions.back() ;
       for ( int qq = 0 ; qq!=contracted_index_positions.size() ; qq++ ){
         if ( (inner_stride) > Tens_strides->at(contracted_index_positions[qq])){
           inner_stride = Tens_strides->at(contracted_index_positions[qq]);
           innermost_ctr = contracted_index_positions[qq];
         }
       }

       cout << "inner_stride  = " << inner_stride << endl;
       cout << "innermost_ctr = " << innermost_ctr << endl;

       shared_ptr<vector<int>> mins2 = make_shared<vector<int>>(maxs->size(), 0); 
       shared_ptr<vector<int>> fvec2 = make_shared<vector<int>>(*mins); 
       shared_ptr<vector<int>> maxs2 = make_shared<vector<int>>(range_sizes.size(),0);
       // note: copying inner block as vector, so want max = min = 0 for those indexes  
       // looks backwards, but is countered by use of f2b fvec routine
       for (int jj = innermost_ctr ; jj != maxs2->size(); jj++) 
         maxs2->at(jj) = range_sizes[jj]-1; 

       shared_ptr<vector<int>> CTens_strides = get_CTens_strides( range_sizes, contracted_index_positions );

       int ctr_block_len = range_sizes[contracted_index_positions.front()];          
       int out_block_size = 1;
       for ( int pos : contracted_index_positions) 
          out_block_size *= range_sizes[pos];

       unique_ptr<DataType[]>      Tens_data_block_old = Tens_in->get_block(Tens_id_blocks_old);
       std::unique_ptr<DataType[]> Tens_data_block_new(new DataType[out_block_size]);
       std::fill_n(Tens_data_block_new.get(), out_block_size, 0.0);
       
       for ( int jj = 0 ; jj != ctr_block_len ; jj++ ) {

         fill(fvec2->begin(), fvec2->end(), 0);
         for ( int ctr_idx_pos : contracted_index_positions ) {
           maxs2->at(ctr_idx_pos) = jj; 
           mins2->at(ctr_idx_pos) = jj;  
           fvec2->at(ctr_idx_pos) = jj; 
         }

         //within index block, loop through data copying chunks where contracted_ctrs are equal
         do { 
           int pos   = inner_product( fvec2->begin(), fvec2->end(), CTens_strides->begin(), 0); 
           int inpos = inner_product( fvec2->begin(), fvec2->end(), Tens_strides->begin(),  0); 
           cout << "pos   = "<< pos <<   endl;
           cout << "inpos = "<< inpos <<   endl;
            
           blas::ax_plus_y_n(1.0, Tens_data_block_old.get()+inpos, inner_stride, Tens_data_block_new.get()+pos);
  
      
         } while (fvec_cycle_skipper_f2b(fvec2, maxs2, mins2));
       }
 
       Tens_out->add_block( Tens_data_block_new, Tens_id_blocks_new );
     } while (fvec_cycle_skipper(block_pos, maxs, mins ));
  } 

  return Tens_out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Code is much clearer if we have a seperate routine for this rather than pepper the generic version with ifs
// Could be accomplished with dot_product member of Tensor....
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_vectors( shared_ptr<Tensor_<DataType>> Tens1_in,
                                                                           shared_ptr<Tensor_<DataType>> Tens2_in ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::contract_vectors" <<endl; 

  vector<IndexRange> T1_org_rngs = Tens1_in->indexrange();
  vector<IndexRange> T2_org_rngs = Tens2_in->indexrange();

  assert(Tens1_in->size_alloc() == Tens2_in->size_alloc() );
  assert(T1_org_rngs.size() == T2_org_rngs.size() );

  DataType dot_product_value = 0.0 ;
  
  for (int ii = 0 ; ii != T1_org_rngs[0].range().size(); ii++ ){ 
     
      vector<Index> T1_block = { T1_org_rngs[0].range(ii) };
      vector<Index> T2_block = { T2_org_rngs[0].range(ii) }; 

      unique_ptr<DataType[]> T1_data_block = Tens1_in->get_block(T1_block); 
      unique_ptr<DataType[]> T2_data_block = Tens2_in->get_block(T2_block); 

      cout << "T1_data, T2_data  =  ";
      for ( int qq = 0 ; qq != T1_block[0].size() ; qq++ ) {
        cout <<  *(T1_data_block.get() + qq) << "    " <<  *(T2_data_block.get() + qq) <<  endl;
      }  
    
      dot_product_value += ddot_( T1_block[0].size(), T1_data_block.get(), 1, T2_data_block.get(), 1 ); 

      cout << endl << "dot_product_value = " << dot_product_value <<  endl;
  }
  
  return dot_product_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contracts a tensor with a vector
// Tens1 is the tensor ( order >1 ) 
// Tens2 is the vector 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_tensor_with_vector( shared_ptr<Tensor_<DataType>> Tens1_in,
                                                                             shared_ptr<Tensor_<DataType>> Tens2_in,
                                                                             pair<int,int> ctr_todo){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::contract_tensor_with_vector" <<endl; 

  assert(ctr_todo.second == 0 ) ;

  vector<IndexRange> T1_org_rngs = Tens1_in->indexrange();
  vector<IndexRange> T2_rng = Tens2_in->indexrange();

  shared_ptr<vector<int>> T1_org_order= make_shared<vector<int>>(T1_org_rngs.size());
  iota(T1_org_order->begin(), T1_org_order->end(), 0);

  //Fortran clomun-major ordering, swap indexes here, not later... 
  shared_ptr<vector<int>>        T1_new_order = put_ctr_at_back( T1_org_order, ctr_todo.first);
  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs) ;

  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin(), T1_new_rngs->end()-1);

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    cout << "T1 block pos =  [ " ;    for (int block_num : *block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( block_pos, *T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector( T1_new_order, T1_new_rng_blocks); 
    vector<Index> T_out_rng_blocks(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
    
    int ctr_block_size    = T1_new_rng_blocks->back().size(); cout << "ctr_block_size = " << ctr_block_size << endl; 
    int T1_block_size     = get_block_size( T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 
    int T_out_block_size = T1_block_size/ctr_block_size; 

    std::unique_ptr<DataType[]> T1_data_new;
    {
      std::unique_ptr<DataType[]> T1_data_org = Tens1_in->get_block(*T1_org_rng_blocks);
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }

    unique_ptr<DataType[]> T2_data = Tens2_in->get_block(vector<Index> { T2_rng[0].range(block_pos->back()) } ); 
 
    std::unique_ptr<DataType[]> T_out_data(new DataType[T_out_block_size]);
    std::fill_n(T_out_data.get(), T_out_block_size, 0.0);
    
    DataType dblone =  1.0;
    int int_one = 1; 
    

    dgemv_( "N", T_out_block_size, ctr_block_size, dblone, T1_data_new.get(), T_out_block_size, T2_data.get(), 1, 
             dblone, T_out_data.get(), int_one );  

    Tens_out->add_block( T_out_data, T_out_rng_blocks );

  } while (fvec_cycle_test(block_pos, maxs1 ));

  
  return Tens_out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Contracts tensors T1 and T2 over two specified indexes
//T1_org_rg T2_org_rg are the original index ranges for the tensors (not necessarily normal ordered).
//T2_new_rg T1_new_rg are the new ranges, with the contracted index at the end, and the rest in normal ordering.
//T1_new_order and T2_new_order are the new order of indexes, and are used for rearranging the tensor data.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_different_tensors( shared_ptr<Tensor_<DataType>> Tens1_in,
                                                                            shared_ptr<Tensor_<DataType>> Tens2_in,
                                                                            pair<int,int> ctr_todo){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::contract_on_different_tensor" <<endl; 

  vector<IndexRange> T1_org_rngs = Tens1_in->indexrange();
  vector<IndexRange> T2_org_rngs = Tens2_in->indexrange();

  if ( ( T1_org_rngs.size() == 1 ) && (T2_org_rngs.size() == 1) ) {
    cout << "vector dot product" << endl;
  } else if ( ( T1_org_rngs.size() == 1 ) && (T2_org_rngs.size() != 1) ) {
    cout << "tensor dot vector " << endl;
  } else if ( ( T1_org_rngs.size() == 1 ) && (T2_org_rngs.size() == 1) ) {
    cout << "vector dot tensor " << endl;
  }    

  shared_ptr<vector<int>> T1_org_order= make_shared<vector<int>>(T1_org_rngs.size());
  iota(T1_org_order->begin(), T1_org_order->end(), 0);

  shared_ptr<vector<int>> T2_org_order= make_shared<vector<int>>(T2_org_rngs.size());
  iota(T2_org_order->begin(), T2_org_order->end(), 0);

  shared_ptr<vector<int>>        T1_new_order = put_ctr_at_back( T1_org_order, ctr_todo.first);
  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs) ;

  shared_ptr<vector<int>>        T2_new_order = put_ctr_at_front( T2_org_order, ctr_todo.second);
  shared_ptr<vector<IndexRange>> T2_new_rngs  = reorder_vector(T2_new_order, T2_org_rngs);
  shared_ptr<vector<int>>        maxs2        = get_num_index_blocks_vec(T2_new_rngs) ;

  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs.insert(Tout_unc_rngs.end(), T2_new_rngs->begin()+1, T2_new_rngs->end());

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    cout << "T1 block pos =  [ " ;    for (int block_num : *T1_rng_block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, *T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector( T1_new_order, T1_new_rng_blocks); 
    
    size_t ctr_block_size    = T1_new_rng_blocks->back().size(); 
    size_t T1_unc_block_size = get_block_size( T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1); 
    size_t T1_block_size     = get_block_size( T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 

    std::unique_ptr<DataType[]> T1_data_new;
    {
      std::unique_ptr<DataType[]> T1_data_org = Tens1_in->get_block(*T1_org_rng_blocks);
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }
   
    cout <<" T1_data_block = [ " ; cout.flush();
    for ( int qq = 0 ; qq!=T1_unc_block_size; qq++ ) 
      cout << *(T1_data_new.get()) << " " ;
    cout << "] " << endl;
  
    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    T2_rng_block_pos->front() = T1_rng_block_pos->back();
  
    do { 

      cout << "T2 block pos =  [ " ; for (int block_num : *T2_rng_block_pos ){cout << block_num << " ";cout.flush();} cout << " ] " << endl;
      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, *T2_new_rngs); 
      shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector( T2_new_order, T2_new_rng_blocks); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      std::unique_ptr<DataType[]> T2_data_new;
      {
        std::unique_ptr<DataType[]> T2_data_org = Tens2_in->get_block(*T2_org_rng_blocks);
        T2_data_new = reorder_tensor_data(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
      }
      
      cout <<" T2_data_block = [ "; cout.flush(); 
      for ( int qq = 0 ; qq!=T2_unc_block_size; qq++ ) 
        cout << *(T2_data_new.get()+qq) << " " ;
      cout << "] " << endl;
 
      std::unique_ptr<DataType[]> T_out_data(new DataType[T1_unc_block_size*T2_unc_block_size]);
      std::fill_n(T_out_data.get(), T1_unc_block_size*T2_unc_block_size, 0.0);

      //should not use transpose; instead build T2_new_order backwards... 
      dgemm_( "N", "N", T1_unc_block_size, T2_unc_block_size, ctr_block_size, 1.0, T1_data_new, T1_unc_block_size,
              T2_data_new, ctr_block_size, 0.0, T_out_data, T1_unc_block_size );

      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      Tens_out->add_block( T_out_data, T_out_rng_block );

      cout <<" T_out_data = [ " ; cout.flush();
      for ( int qq = 0 ; qq!=T1_unc_block_size*T2_unc_block_size; qq++ ) 
        cout << *(T_out_data.get()+qq) << " " ;
      cout << "] " << endl;

      //remove last index; contracted index is cycled in T1 loop
    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, 0 )) ;

  } while (fvec_cycle_test(T1_rng_block_pos, maxs1 ));

  
  return Tens_out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sets all elements of input tensor Tens to the value specified by elem_val 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems(shared_ptr<Tensor_<DataType>> Tens, DataType elem_val  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::set_tensor_elems  " << endl;
  
   vector<IndexRange> id_ranges = Tens->indexrange();

   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges ) ; 

//   cout << " range_lengths = [ " ; cout.flush(); for ( int pos : *range_lengths ) { cout << pos << " "  ; cout.flush(); } cout << "] " << endl;

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
     // cout << " block_pos = [ " ; cout.flush(); for ( int pos : *block_pos ) { cout << pos << " "  ; cout.flush(); } cout << "] " << endl;
     vector<Index> id_blocks =  *(get_rng_blocks( block_pos, id_ranges ));
 
     if ( Tens->exists(id_blocks) ) { 
       unique_ptr<DataType[]> block_data = Tens->get_block(id_blocks);
       std::fill_n(block_data.get(), Tens->get_size(id_blocks), elem_val);
       Tens->put_block(block_data, id_blocks);
     }
   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor(shared_ptr<Tensor_<DataType>> Tens_in,  shared_ptr<vector<int>> new_order){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::reorder_block_Tensor "; cout.flush();

   shared_ptr<vector<IndexRange>> T_id_ranges = make_shared<vector<IndexRange>>(Tens_in->indexrange());
   shared_ptr<vector<int>> range_lengths = make_shared<vector<int>>(0); 
   for (auto idrng : *T_id_ranges )
     range_lengths->push_back(idrng.range().size()-1); 
   
   shared_ptr<vector<IndexRange>> reordered_ranges    = reorder_vector(new_order, T_id_ranges ) ;
   shared_ptr<Tensor_<DataType>> reordered_block_tensor = make_shared<Tensor_<DataType>>(*reordered_ranges);
   reordered_block_tensor->allocate();
   reordered_block_tensor->zero();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins      = make_shared<vector<int>>(T_id_ranges->size(),0);  
   do {
     cout << " block pos =  [ " ; for(int block_num : *block_pos ){ cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
     
     shared_ptr<vector<Index>> orig_id_blocks      = get_rng_blocks( block_pos, *T_id_ranges); 
     shared_ptr<vector<Index>> reordered_id_blocks = reorder_vector(new_order, orig_id_blocks ) ;

     unique_ptr<DataType[]> reordered_data_block; 
     {
     unique_ptr<DataType[]> orig_data_block = Tens_in->get_block(*orig_id_blocks);
     reordered_data_block = reorder_tensor_data( orig_data_block.get(), new_order, orig_id_blocks ) ;
     }
      
     reordered_block_tensor->put_block(reordered_data_block, *reordered_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));

   return reordered_block_tensor;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_tensor_data( const DataType* orig_data, shared_ptr<vector<int>>  new_order_vec,
                                                           shared_ptr<vector<Index>> orig_index_blocks ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "reorder_tensor_data" << endl;

  shared_ptr<vector<size_t>> rlen = get_sizes(orig_index_blocks);
  shared_ptr<vector<size_t>> new_order_st = make_shared<vector<size_t>>(new_order_vec->size());   
  size_t block_size = get_block_size(orig_index_blocks->begin(), orig_index_blocks->end());
  array<int,4> sort_options = {0,1,1,1};

  unique_ptr<DataType[]> reordered_data(new DataType[block_size]);

   for ( int ii = 0 ; ii != new_order_vec->size(); ii++) 
     new_order_st->at(ii) = new_order_vec->at(ii);

  Tensor_Sorter::Tensor_Sorter<DataType> TS ;

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with ranges specified by unc_ranges, where all values are equal to XX  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_Tensor(shared_ptr<vector<IndexRange>> T_id_ranges, DataType XX ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::get_uniform_Tensor" << endl;

   shared_ptr<vector<int>>  range_lengths  = make_shared<vector<int>>(0); 
   for ( IndexRange idrng : *T_id_ranges )
      range_lengths->push_back(idrng.range().size()-1); 

   shared_ptr<Tensor_<DataType>> block_tensor = make_shared<Tensor_<DataType>>(*T_id_ranges);
   block_tensor->allocate();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
   do {

     vector<Index> T_id_blocks(T_id_ranges->size());
     for( int ii = 0 ;  ii != T_id_blocks.size(); ii++)
       T_id_blocks[ii] =  T_id_ranges->at(ii).range(block_pos->at(ii));
     
     int out_size = 1;
     for ( Index id : T_id_blocks)
        out_size*= id.size();

     unique_ptr<DataType[]> T_block_data( new DataType[out_size] );
   
     DataType* dit = T_block_data.get() ;
     for ( int qq =0 ; qq != out_size; qq++ ) 
        T_block_data[qq] = XX; 

     block_tensor->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));
   return block_tensor;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_block_of_data( DataType* data_ptr ,
                                                         shared_ptr<vector<IndexRange>> id_ranges, 
                                                         shared_ptr<vector<int>> block_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::get_block_of_data" << endl;

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

  // getting size of ranges (seems to be correctly offset for node)
  vector<size_t> range_sizes(block_pos->size());
  for (int ii = 0 ; ii != id_ranges->size() ; ii++)
    range_sizes[ii]  = id_ranges->at(ii).size();

  size_t data_block_size = 1;
  size_t data_block_pos  = 1;
  for (int ii = 0 ; ii != id_ranges->size()-1 ; ii++){
    data_block_pos  *= id_pos[ii]*range_sizes[ii];
    data_block_size *= id_ranges->at(ii).range(block_pos->at(ii)).size();
  }

  data_block_size *= id_ranges->back().range(block_pos->back()).size();
 
  cout << "data_block_size = " << data_block_size << endl;
  cout << "data_block_pos = "  << data_block_pos << endl;
 
  unique_ptr<DataType[]> data_block(new DataType[data_block_size])  ;

  copy_n(data_block.get(), data_block_size, data_ptr+data_block_pos);

  return data_block; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Tensor_Arithmetic::Tensor_Arithmetic<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
 
#endif
