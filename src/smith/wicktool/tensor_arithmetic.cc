#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/tensor_arithmetic.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;
using namespace bagel::SMITH::Tensor_Arithmetic_Utils; 
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Specialized routine for summing over the whole tensor, should not be needed by handy for now
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType Tensor_Arithmetic::Tensor_Arithmetic<DataType>::sum_tensor_elems( shared_ptr<Tensor_<DataType>> Tens_in) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithemetic_Utils::sum_tensor_elems" << endl;

  shared_ptr<vector<IndexRange>> id_ranges = make_shared<vector<IndexRange>>( Tens_in->indexrange() );
  
  shared_ptr<vector<int>> range_maxs  =  get_range_lengths( id_ranges ) ;
  shared_ptr<vector<int>> block_pos   =  make_shared<vector<int>>(range_maxs->size(),0);  
  shared_ptr<vector<int>> mins        =  make_shared<vector<int>>(range_maxs->size(),0);  

  double sum_of_elems = (DataType)0.0;

  do { 

     vector<Index> id_blocks = *(get_rng_blocks( block_pos, *id_ranges ));
     unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
     Tens_in->get_size(id_blocks);
     double* bob = block.get();
     sum_of_elems += *bob++;

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return sum_of_elems;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Bad routine using reordering because I just want something which works
// All orders of ctrs_pos will give the same answer, but some will require fewer transposes. 
// Write (small) routine to optimize ordering of ctrs_pos and call at start of routine.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_on_same_tensor( shared_ptr<Tensor_<DataType>> Tens_in,  vector<int>& ctrs_pos) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::contract_on_same_tensor" << endl;

  vector<IndexRange> id_ranges_in = Tens_in->indexrange();
  int num_ids  = id_ranges_in.size();

  // Put contracted indexes at back so they change slowly; Fortran ordering in index blocks!
  shared_ptr<vector<int>> new_order = make_shared<vector<int>>(num_ids);
  iota(new_order->begin(), new_order->end(), 0);
  put_ctrs_at_back( *new_order, ctrs_pos );

  vector<int> unc_pos( new_order->begin(), new_order->end() - ctrs_pos.size() );
  vector<IndexRange> unc_idxrng( unc_pos.size() );
  for (int qq = 0; qq != unc_pos.size(); qq++ ) 
    unc_idxrng[qq] = id_ranges_in[unc_pos[qq]];

  shared_ptr<vector<int>> maxs      = get_range_lengths( id_ranges_in ) ;
  shared_ptr<vector<int>> mins      = make_shared<vector<int>>(maxs->size(),0);  
  shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(maxs->size(),0);  

  //set max = min so contracted blocks are skipped in fvec_cycle
  for ( vector<int>::iterator iter = ctrs_pos.begin(); iter != ctrs_pos.end(); iter++ ) {
     maxs->at(*iter) = 0;
     mins->at(*iter) = 0;
  } 

  int num_ctr_blocks = id_ranges_in[ctrs_pos.front()].range().size();

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(unc_idxrng);
  Tens_out->allocate();
  Tens_out->zero();

  // Outer forvec loop skips ctr ids due to min[i]=max[i].  Inner loop cycles ctr ids.
  do {
   
    shared_ptr<vector<Index>> id_blocks_in = get_rng_blocks( block_pos, id_ranges_in ); 

    int unc_block_size = 1;
    vector<Index> unc_id_blocks(unc_pos.size());
    for (int qq = 0 ; qq != unc_pos.size(); qq++ ) {
      unc_id_blocks[qq] = id_blocks_in->at(unc_pos[qq]);
      unc_block_size   *= unc_id_blocks[qq].size();
    }

    unique_ptr<DataType[]> contracted_block(new DataType[unc_block_size]);
    fill_n(contracted_block.get(), unc_block_size, 0.0);
    //loop over ctr blocks
    for ( int ii = 0 ; ii != num_ctr_blocks; ii++) { 
      
      for ( int pos : ctrs_pos )  
         block_pos->at(pos) = ii;

      //redefine block position
      id_blocks_in = get_rng_blocks(block_pos, id_ranges_in);

      vector<int> ctr_block_strides(ctrs_pos.size(),1);
      ctr_block_strides.front() = unc_block_size;
      int ctr_total_stride = unc_block_size;
      for (int qq = 1; qq != ctrs_pos.size() ; qq++ ) {
        ctr_block_strides[qq] = (ctr_block_strides[qq-1] * id_blocks_in->at(ctrs_pos[qq-1]).size());
        ctr_total_stride += ctr_block_strides[qq]; 
      }

      {
      unique_ptr<DataType[]> orig_block = Tens_in->get_block(*id_blocks_in);
      unique_ptr<DataType[]> reordered_block = reorder_tensor_data( orig_block.get(), new_order, id_blocks_in );

      //looping over the id positions within the block
      for ( int ctr_id = 0 ; ctr_id != id_blocks_in->at(ctrs_pos[0]).size(); ctr_id++ )
        blas::ax_plus_y_n(1.0, reordered_block.get() + (ctr_total_stride*ctr_id), unc_block_size, contracted_block.get() );
      
      }
       
    }
    
    Tens_out->put_block(contracted_block, unc_id_blocks );
   
  } while( fvec_cycle_skipper(block_pos, maxs, mins));

  return Tens_out;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Lazy, so can substitute this in to test easily, should remove later
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_on_same_tensor( shared_ptr<Tensor_<DataType>> Tens_in,  pair<int,int> ctrs_pair) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::contract_on_same_tensor , pair " << endl;

   vector<int> ctrs_pos = { ctrs_pair.first, ctrs_pair.second };

   return contract_on_same_tensor( Tens_in, ctrs_pos) ;
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

  //Fortran column-major ordering, swap indexes here, not later... 
  shared_ptr<vector<int>>        T1_new_order = put_ctr_at_front( T1_org_order, ctr_todo.first);
  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs) ;

  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin()+1, T1_new_rngs->end());

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    cout << "T1 block pos =  [ " ;    for (int block_num : *block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( block_pos, *T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector( T1_new_order, T1_new_rng_blocks); 
    vector<Index> T_out_rng_blocks(T1_new_rng_blocks->begin()+1, T1_new_rng_blocks->end());
    
    int ctr_block_size    = T1_new_rng_blocks->front().size(); cout << " ctr_block_size = " << ctr_block_size << endl; 
    int T1_block_size     = get_block_size( T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 
    int T_out_block_size  = T1_block_size/ctr_block_size;      cout << " T_out_block_size = " << T_out_block_size << endl; 

    std::unique_ptr<DataType[]> T1_data_new;
    {
      std::unique_ptr<DataType[]> T1_data_org = Tens1_in->get_block(*T1_org_rng_blocks);
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }

    unique_ptr<DataType[]> T2_data = Tens2_in->get_block(vector<Index> { T2_rng[0].range(block_pos->front()) } ); 
 
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
cout << "Tensor_Arithmetic::contract_on_different_tensor_column_major" <<endl; 

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
 
  // front for column major ordering
  shared_ptr<vector<int>>        T1_new_order = put_ctr_at_front( T1_org_order, ctr_todo.first);
  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs) ;

  shared_ptr<vector<int>>        T2_new_order = put_ctr_at_back( T2_org_order, ctr_todo.second);
  shared_ptr<vector<IndexRange>> T2_new_rngs  = reorder_vector(T2_new_order, T2_org_rngs);
  shared_ptr<vector<int>>        maxs2        = get_num_index_blocks_vec(T2_new_rngs) ;

   cout << "established ordering and got ranges" <<endl;
 
  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin()+1, T1_new_rngs->end());
  Tout_unc_rngs.insert(Tout_unc_rngs.end(), T2_new_rngs->begin(), T2_new_rngs->end()-1);

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    
    print_vector( *T1_rng_block_pos , "T1_block_pos" ); cout << endl;

    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, *T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector( T1_new_order, T1_new_rng_blocks); 
    
    size_t ctr_block_size    = T1_new_rng_blocks->front().size(); 
    size_t T1_unc_block_size = get_block_size( T1_new_rng_blocks->begin()+1, T1_new_rng_blocks->end()); 
    size_t T1_block_size     = get_block_size( T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 

    std::unique_ptr<DataType[]> T1_data_new;
    {
      std::unique_ptr<DataType[]> T1_data_org = Tens1_in->get_block(*T1_org_rng_blocks);
      cout << "got T1 block, reordering to  "; print_vector(*T1_new_order , "T1_new_order" ); cout << endl;
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
      cout << "reordered T1 block" << endl;
    }
   
    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    T2_rng_block_pos->back() = T1_rng_block_pos->front();
  
    do { 

      print_vector( *T2_rng_block_pos , "T2_block_pos" ); cout << endl;

      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, *T2_new_rngs); 
      shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector( T2_new_order, T2_new_rng_blocks); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin(), T2_new_rng_blocks->end()-1);



      std::unique_ptr<DataType[]> T2_data_new;
      {
        cout << "getting T2 block" << endl;
        std::unique_ptr<DataType[]> T2_data_org = Tens2_in->get_block(*T2_org_rng_blocks);
        cout << "got T2 block, reordering to  "; print_vector(*T2_new_order , "T2_new_order" ); cout << endl;
        T2_data_new = reorder_tensor_data(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
        cout << "reordered T2 block" << endl;
      }
      
      std::unique_ptr<DataType[]> T_out_data(new DataType[T1_unc_block_size*T2_unc_block_size]);
      std::fill_n(T_out_data.get(), T1_unc_block_size*T2_unc_block_size, 0.0);

      dgemm_( "N", "N", T1_unc_block_size, T2_unc_block_size, ctr_block_size, 1.0, T1_data_new, T1_unc_block_size,
              T2_data_new, ctr_block_size, 0.0, T_out_data, T1_unc_block_size );

      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin()+1, T1_new_rng_blocks->end());
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin(), T2_new_rng_blocks->end()-1);
      cout << "putting block" << endl;
      Tens_out->add_block( T_out_data, T_out_rng_block );
      cout << "put block" << endl;

    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, 0 )) ;//remove last index; contracted index is cycled in T1 loop

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
   shared_ptr<vector<int>> range_lengths = get_range_lengths( T_id_ranges ); 
   
   shared_ptr<vector<IndexRange>> reordered_ranges       = reorder_vector(new_order, T_id_ranges ) ;
   shared_ptr<Tensor_<DataType>>  reordered_block_tensor = make_shared<Tensor_<DataType>>(*reordered_ranges);
   reordered_block_tensor->allocate();
   reordered_block_tensor->zero();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins      = make_shared<vector<int>>(T_id_ranges->size(),0);  
   do {
     cout << " block pos =  [ " ; for(int block_num : *block_pos ){ cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
     
     shared_ptr<vector<Index>> orig_id_blocks      = get_rng_blocks( block_pos, *T_id_ranges ); 

     unique_ptr<DataType[]> reordered_data_block; 
     {
     unique_ptr<DataType[]> orig_data_block = Tens_in->get_block( *orig_id_blocks );
     reordered_data_block = reorder_tensor_data( orig_data_block.get(), new_order, orig_id_blocks );
     }
      
     shared_ptr<vector<Index>> reordered_id_blocks = reorder_vector( new_order, orig_id_blocks );
     reordered_block_tensor->put_block( reordered_data_block, *reordered_id_blocks );

   } while ( fvec_cycle_skipper( block_pos, range_lengths, mins ) );

   cout << "reordered_block_tensor" << endl;   
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
  fill_n(reordered_data.get(), block_size, 0.0 );

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
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_Tensor(shared_ptr<vector<IndexRange>> T_id_ranges, DataType XX ){
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with element with distinct elements, all blcoks have similarly generated elems though
//C ordering (row major)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_Tensor_row_major(shared_ptr<vector<IndexRange>> T_id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::get_test_Tensor" << endl;

   shared_ptr<Tensor_<DataType>> Tens = make_shared<Tensor_<DataType>>(*T_id_ranges);
   Tens->allocate();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( T_id_ranges ); 

   double put_after_decimal_point  = pow(10 , block_pos->size());
   vector<double> power_10( block_pos->size() );
   for (int ii = 0  ; ii != power_10.size() ; ii++ ) 
     power_10[power_10.size()-1-ii] = pow(10, (double)ii );///put_after_decimal_point ;

   do {

     vector<Index> T_id_blocks = *(get_rng_blocks( block_pos, *T_id_ranges )); 
     size_t out_size = Tens->get_size(T_id_blocks); 

     unique_ptr<DataType[]> T_block_data( new DataType[out_size] );
  
     shared_ptr<vector<int>> id_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
     shared_ptr<vector<int>> id_mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
     shared_ptr<vector<int>> id_maxs = make_shared<vector<int>>(T_id_ranges->size());
     for (int xx = 0 ; xx !=id_maxs->size(); xx++ ) 
       id_maxs->at(xx)  =   T_id_blocks[xx].size()-1;

     int  qq =0;
     //TODO add in offsets at some point
     do {
        T_block_data[qq++] = inner_product (id_pos->begin(), id_pos->end(), power_10.begin(), 0 );
     } while (fvec_cycle( id_pos, id_maxs, id_mins ));

     Tens->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));
   return Tens;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with element with distinct elements, all blcoks have similarly generated elems though
//C ordering (row major)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_Tensor_column_major(shared_ptr<vector<IndexRange>> T_id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::get_test_Tensor" << endl;

   shared_ptr<Tensor_<DataType>> Tens = make_shared<Tensor_<DataType>>(*T_id_ranges);
   Tens->allocate();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( T_id_ranges ); 

   double put_after_decimal_point  = pow(10 , block_pos->size());
   vector<double> power_10( block_pos->size() );
   for (int ii = 0  ; ii != power_10.size() ; ii++ ) 
     power_10[ii] = pow(10, (double)ii );///put_after_decimal_point ;

   do {

     vector<Index> T_id_blocks = *(get_rng_blocks( block_pos, *T_id_ranges )); 
     size_t out_size = Tens->get_size(T_id_blocks); 

     unique_ptr<DataType[]> T_block_data( new DataType[out_size] );
  
     shared_ptr<vector<int>> id_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
     shared_ptr<vector<int>> id_mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
     shared_ptr<vector<int>> id_maxs = make_shared<vector<int>>(T_id_ranges->size());
     for (int xx = 0 ; xx !=id_maxs->size(); xx++ ) 
       id_maxs->at(xx)  =   T_id_blocks[xx].size()-1;

     int  qq =0;
     //TODO add in offsets at some point
     do {
        T_block_data[qq++] = inner_product (id_pos->begin(), id_pos->end(), power_10.begin(), 0 );
     } while (fvec_cycle( id_pos, id_maxs, id_mins ));

     Tens->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));
   return Tens;
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
