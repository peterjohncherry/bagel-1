#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::Tensor_Sorter;
using namespace bagel::Tensor_Arithmetic_Utils; 
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

  DataType sum_of_elems = (DataType)0.0;

  do { 

     vector<Index> id_blocks = *(get_rng_blocks( block_pos, *id_ranges ));
     unique_ptr<DataType[]> block = Tens_in->get_block( id_blocks );
     Tens_in->get_size(id_blocks);
     DataType* bob = block.get();
     sum_of_elems += *bob++;

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return sum_of_elems;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Full takes the trace of the tensor, i.e. sets all indexes equal. Use different routine as reordering of indexes is not necessary here.
// Returns a number.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::trace_tensor__number_return( shared_ptr<Tensor_<DataType>> Tens_in ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::trace_tensor__number_return" << endl;
  assert ( Tens_in->indexrange().size() > 0 );

#ifndef NDEBUG
  for ( int ii =1 ; ii != Tens_in->indexrange().size(); ii++ ) 
    assert( Tens_in->indexrange()[ii].size() == Tens_in->indexrange()[ii-1].size() ) ; 
# endif 

  vector<IndexRange> id_ranges = Tens_in->indexrange();

  DataType Tens_trace = 0.0;

  for ( int ctr_block_pos = 0 ; ctr_block_pos != id_ranges[0].range().size(); ctr_block_pos++ ) {
  
    vector<Index> id_blocks( id_ranges.size() );
    for ( int qq = 0; qq != id_ranges.size(); qq++ ) 
      id_blocks[qq] = id_ranges[qq].range(ctr_block_pos);
      
    int ctr_total_stride = 1;
    unique_ptr<DataType[]> data_block = Tens_in->get_block( id_blocks ) ; 
    for ( int ctr_id = 0 ; ctr_id != id_blocks[0].size(); ctr_id++ )
      Tens_trace += *(data_block.get() + (ctr_total_stride*ctr_id));
   
  } 

  return Tens_trace;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Full takes the trace of the tensor, i.e. sets all indexes equal. Use different routine as reordering of indexes is not necessary here.
// Returns Tensor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::trace_tensor__tensor_return( shared_ptr<Tensor_<DataType>> Tens_in ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::trace_tensor__tensor_return" << endl;
  assert ( Tens_in->indexrange().size() > 1 );
 
#ifndef NDEBUG
  for ( int ii =1 ; ii != Tens_in->indexrange().size(); ii++ ) 
    assert( Tens_in->indexrange()[ii].size() == Tens_in->indexrange()[ii-1].size() ) ; 
# endif 

  vector<IndexRange> id_ranges = Tens_in->indexrange();
  vector<IndexRange> unit_range = { IndexRange( 1, 1, 0, 1 ) };

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(unit_range);
  Tens_out->allocate();
  DataType Tens_trace = 0.0; 

  for ( int ctr_block_pos = 0 ; ctr_block_pos != id_ranges[0].range().size(); ctr_block_pos++ ) {
   
    vector<Index> id_blocks( id_ranges.size() );
    for ( int qq = 0; qq != id_ranges.size(); qq++ ) 
      id_blocks[qq] = id_ranges[qq].range(ctr_block_pos);
      

    int ctr_total_stride = 1;
    unique_ptr<DataType[]> data_block = Tens_in->get_block( id_blocks ) ; 
    for ( int ctr_id = 0 ; ctr_id != id_blocks[0].size(); ctr_id++ ) {
      Tens_trace += *(data_block.get() + (ctr_total_stride*ctr_id));
    }
 
  } 

  unique_ptr<DataType[]> data_block( new DataType[1] ); 
  data_block[0] = Tens_trace;

  vector<Index> unit_index_block = { unit_range[0].range(0) }; 
  Tens_out->put_block( data_block, unit_index_block );

  cout << " out of Tensor trace contract" << endl;
  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sum over indexes on a tensor, is used for excitation operators. Note the output tensor may have an odd number of indexes,
//as a rule, this may mess up some of the other routines....
//Currently using reordering, as other wise there's going to be a huge amount of hopping about (weird stride). Don't know what is best.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::sum_over_idxs( shared_ptr<Tensor_<DataType>> Tens_in, vector<int>& summed_idxs_pos) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::contract_on_same_tensor" << endl;
  assert ( Tens_in->indexrange().size() > 0 );

  vector<IndexRange> id_ranges = Tens_in->indexrange();
  shared_ptr<Tensor_<DataType>> Tens_out ;

  int num_ids  = id_ranges.size();
  int num_ctrs = summed_idxs_pos.size();
  int num_uncs = num_ids - num_ctrs;
  
  // Put contracted indexes at back so they change slowly; Fortran ordering in index blocks!
  shared_ptr<vector<int>> new_order = make_shared<vector<int>>(num_ids);
  iota(new_order->begin(), new_order->end(), 0);
  put_ctrs_at_back( *new_order, summed_idxs_pos );


  vector<int> unc_pos( new_order->begin(), new_order->end() - num_ctrs );
  vector<IndexRange> unc_idxrng( unc_pos.size() );
  for ( int qq = 0; qq != unc_pos.size(); qq++ )
    unc_idxrng[qq] = id_ranges[unc_pos[qq]];
  
  shared_ptr<vector<int>> unc_maxs      = get_range_lengths( id_ranges );
  shared_ptr<vector<int>> unc_mins      = make_shared<vector<int>>(unc_maxs->size(),0);
  shared_ptr<vector<int>> unc_block_pos = make_shared<vector<int>>(unc_maxs->size(),0);
  
  shared_ptr<vector<int>> ctr_maxs      = make_shared<vector<int>>(*unc_maxs);
  shared_ptr<vector<int>> ctr_mins      = make_shared<vector<int>>(unc_maxs->size(),0);
  //set max = min so contracted blocks are skipped in fvec_cycle
  for ( int qq = num_ctrs-1; qq!= num_ids; qq++ ) {
    unc_maxs->at(qq) = 0;
    unc_mins->at(qq) = 0;
  }
  
  Tens_out = make_shared<Tensor_<DataType>>(unc_idxrng);
  Tens_out->allocate();
  Tens_out->zero();
  
  // Outer forvec loop skips ctr ids due to min[i]=max[i].  Inner loop cycles ctr ids.
  do {
   
    shared_ptr<vector<Index>> id_blocks = get_rng_blocks( unc_block_pos, id_ranges );
  
    int unc_block_size = 1;
    vector<Index> unc_id_blocks(unc_pos.size());
    for (int qq = 0 ; qq != unc_pos.size(); qq++ ) {
      unc_id_blocks[qq] = id_blocks->at(unc_pos[qq]);
      unc_block_size   *= unc_id_blocks[qq].size();
    }
 
    int ctr_block_size = 1;
    for (int qq = 0 ; qq != summed_idxs_pos.size(); qq++ ) 
      ctr_block_size   *= unc_id_blocks[qq].size();


    unique_ptr<DataType[]> contracted_block(new DataType[unc_block_size]);
    fill_n(contracted_block.get(), unc_block_size, 0.0);
    //set max = min so unc ids are skipped due to min = max
    for ( int qq = 0 ; qq != unc_pos.size(); qq++ ) {
      ctr_maxs->at(qq) = unc_block_pos->at(qq);
      ctr_mins->at(qq) = unc_block_pos->at(qq);
    }

    shared_ptr<vector<int>> ctr_block_pos = make_shared<vector<int>>(unc_maxs->size(),0);
    //loop over ctr blocks
    do {

      //redefine block position
      id_blocks = get_rng_blocks(ctr_block_pos, id_ranges);

      {

      unique_ptr<DataType[]> orig_block = Tens_in->get_block(*id_blocks);
      unique_ptr<DataType[]> reordered_block = reorder_tensor_data( orig_block.get(), new_order, id_blocks );
      //looping over the id positions within the block
      for ( int flattened_ctr_id = 0 ; flattened_ctr_id != ctr_block_size; flattened_ctr_id++ )
        blas::ax_plus_y_n(1.0, reordered_block.get() + (unc_block_size*flattened_ctr_id), unc_block_size, contracted_block.get() );
      }
       
    } while( fvec_cycle_skipper(ctr_block_pos, ctr_maxs, ctr_mins ));

    Tens_out->put_block(contracted_block, unc_id_blocks );

  } while( fvec_cycle_skipper(unc_block_pos, unc_maxs, unc_mins));

  return Tens_out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Bad routine using reordering because I just want something which works
// TODO Write (small) routine to optimize ordering of ctrs_pos and call at start of routine.
//      -> All orders of ctrs_pos will give the same answer, but some will require fewer transposes. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_on_same_tensor( shared_ptr<Tensor_<DataType>> Tens_in,  vector<int>& ctrs_pos) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::contract_on_same_tensor" << endl;
  assert ( Tens_in->indexrange().size() > 0 );

#ifndef NDEBUG
  cout << "id_ranges_sizes = [ " << Tens_in->indexrange()[ctrs_pos[0]].size() << " " ;
  for ( int ii =1 ; ii != ctrs_pos.size(); ii++ ) {
    cout << Tens_in->indexrange()[ctrs_pos[ii]].size() << " " ;
    assert( Tens_in->indexrange()[ctrs_pos[ii]].size() == Tens_in->indexrange()[ctrs_pos[ii-1]].size() ) ; 
  } 
  cout << " ] " << endl;
#endif 

  vector<IndexRange> id_ranges_in = Tens_in->indexrange();
  shared_ptr<Tensor_<DataType>> Tens_out ;

  if ( ctrs_pos.size() == Tens_in->indexrange().size() ) {

    DataType Tensor_trace = trace_tensor__number_return( Tens_in ) ;
    cout << "Tensor_trace = " << Tensor_trace << endl;

    Tens_out = trace_tensor__tensor_return( Tens_in ) ;

  } else {

    int num_ids  = id_ranges_in.size();
  
    // Put contracted indexes at back so they change slowly; Fortran ordering in index blocks!
    shared_ptr<vector<int>> new_order = make_shared<vector<int>>(num_ids);
    iota(new_order->begin(), new_order->end(), 0);
    put_ctrs_at_back( *new_order, ctrs_pos );
    
    vector<int> unc_pos( new_order->begin(), new_order->end() - ctrs_pos.size() );
    vector<IndexRange> unc_idxrng( unc_pos.size() );
    for ( int qq = 0; qq != unc_pos.size(); qq++ ) 
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
    
    Tens_out = make_shared<Tensor_<DataType>>(unc_idxrng);
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
  }

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

      dot_product_value += dot_product( T1_block[0].size(), T1_data_block.get(), T2_data_block.get() ); 

      cout << endl << "dot_product_value = " << dot_product_value <<  endl;
  }
  
  return dot_product_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contracts a tensor with a vector
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_tensor_with_vector( shared_ptr<Tensor_<DataType>> TensIn,
                                                                             shared_ptr<Tensor_<DataType>> VecIn,
                                                                             int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::contract_tensor_with_vector" <<endl; 

  vector<IndexRange> TensIn_org_rngs = TensIn->indexrange();
  vector<IndexRange> VecIn_rng = VecIn->indexrange();

  shared_ptr<vector<int>> TensIn_org_order= make_shared<vector<int>>(TensIn_org_rngs.size());
  iota(TensIn_org_order->begin(), TensIn_org_order->end(), 0);

  //Fortran column-major ordering, swap indexes here, not later... 
  shared_ptr<vector<int>>        TensIn_new_order = put_ctr_at_back( TensIn_org_order, ctr_pos);
  shared_ptr<vector<IndexRange>> TensIn_new_rngs  = reorder_vector(TensIn_new_order, TensIn_org_rngs);

  vector<IndexRange> TensOut_rngs(TensIn_new_rngs->begin(), TensIn_new_rngs->end()-1);

  shared_ptr<Tensor_<DataType>> TensOut = make_shared<Tensor_<DataType>>(TensOut_rngs);  
  TensOut->allocate();
  TensOut->zero();

  shared_ptr<vector<int>> maxs = get_num_index_blocks_vec(TensIn_new_rngs) ;
  shared_ptr<vector<int>> mins = make_shared<vector<int>>(maxs->size(),0);

  shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(TensIn_new_order->size(),0);
  int num_ctr_blocks = maxs->back()+1; 

  //This loop looks silly; trying to avoid excessive reallocation, deallocation and zeroing of TensOut_data
  //However, when parallelized, will presumably need to do all that in innermost loop anyway....
  do { 

     shared_ptr<vector<Index>> TensIn_new_rng_blocks = get_rng_blocks( block_pos, *TensIn_new_rngs); 

     shared_ptr<vector<int>>   TensOut_block_pos = make_shared<vector<int>>(block_pos->begin(), block_pos->end()-1);
     shared_ptr<vector<Index>> TensOut_rng_blocks = get_rng_blocks(TensOut_block_pos, TensOut_rngs);

     int TensOut_block_size  =  get_block_size( TensOut_rng_blocks->begin(), TensOut_rng_blocks->end() );

     std::unique_ptr<DataType[]> TensOut_data(new DataType[TensOut_block_size]);
     std::fill_n(TensOut_data.get(), TensOut_block_size, 0.0);

     for ( int ii = 0 ; ii != num_ctr_blocks; ii++) {

       maxs->front() = ii;
       mins->front() = ii;
       block_pos->front() = ii;

       shared_ptr<vector<Index>> TensIn_org_rng_blocks = inverse_reorder_vector( TensIn_new_order, TensIn_new_rng_blocks); 
       
       int TensIn_block_size = get_block_size( TensIn_org_rng_blocks->begin(), TensIn_org_rng_blocks->end()); 
       int ctr_block_size    = TensIn_new_rng_blocks->back().size();  
 
       cout << "-----------Tensor block in old order----------- " << endl;
       std::unique_ptr<DataType[]> TensIn_data_reord;
       {
       std::unique_ptr<DataType[]> TensIn_data_org = TensIn->get_block(*TensIn_org_rng_blocks);
       TensIn_data_reord = reorder_tensor_data(TensIn_data_org.get(), TensIn_new_order, TensIn_org_rng_blocks);
       }
       
       unique_ptr<DataType[]> VecIn_data = VecIn->get_block(vector<Index> { VecIn_rng[0].range(block_pos->front()) } ); 
       
       DataType dblone =  1.0;
       int int_one = 1; 
       
       gemv( 'N', TensOut_block_size, ctr_block_size, TensIn_data_reord.get(), VecIn_data.get(), TensOut_data.get(), dblone, dblone ); 
       
       TensIn_new_rng_blocks = get_rng_blocks( block_pos, *TensIn_new_rngs); 

     }

     TensOut->put_block( TensOut_data, *TensOut_rng_blocks );

  } while (fvec_cycle_skipper(block_pos, maxs, mins ));
  
  return TensOut;
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

  shared_ptr<vector<int>> T1_org_order= make_shared<vector<int>>(T1_org_rngs.size());
  iota(T1_org_order->begin(), T1_org_order->end(), 0);

  shared_ptr<vector<int>> T2_org_order= make_shared<vector<int>>(T2_org_rngs.size());
  iota(T2_org_order->begin(), T2_org_order->end(), 0);
 
  //note column major ordering
  shared_ptr<vector<int>>        T1_new_order = put_ctr_at_back( T1_org_order, ctr_todo.first);
  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs);
  shared_ptr<vector<int>>        mins1        = make_shared<vector<int>>( maxs1->size(), 0 );  

  shared_ptr<vector<int>>        T2_new_order = put_ctr_at_front( T2_org_order, ctr_todo.second);
  shared_ptr<vector<IndexRange>> T2_new_rngs  = reorder_vector(T2_new_order, T2_org_rngs);
  shared_ptr<vector<int>>        maxs2        = get_num_index_blocks_vec(T2_new_rngs) ;
  shared_ptr<vector<int>>        mins2        = make_shared<vector<int>>( maxs2->size(), 0 );  

  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs.insert(Tout_unc_rngs.end(), T2_new_rngs->begin()+1, T2_new_rngs->end());

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    
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
    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    T2_rng_block_pos->front() = T1_rng_block_pos->back();
    mins2->front() = maxs1->back();
    maxs2->front() = maxs1->back();
  
    do { 

      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, *T2_new_rngs); 
      shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector( T2_new_order, T2_new_rng_blocks); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());

      std::unique_ptr<DataType[]> T2_data_new;   
      {
        std::unique_ptr<DataType[]> T2_data_org = Tens2_in->get_block(*T2_org_rng_blocks);  
        T2_data_new = reorder_tensor_data(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
      }
      
      std::unique_ptr<DataType[]> T_out_data(new DataType[T1_unc_block_size*T2_unc_block_size]);
      std::fill_n(T_out_data.get(), T1_unc_block_size*T2_unc_block_size, 0.0);

      gemm( 'N', 'N', T1_unc_block_size, T2_unc_block_size, ctr_block_size, T1_data_new.get(), T2_data_new.get(), T_out_data.get(), (DataType)1.0, (DataType)0.0 ) ;

      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      Tens_out->add_block( T_out_data, T_out_rng_block );

    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, mins2 ));

  } while (fvec_cycle_skipper(T1_rng_block_pos, maxs1, mins1 ));

  
  return Tens_out;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::contract_different_tensors( shared_ptr<Tensor_<DataType>> Tens1_in,
                                                                            shared_ptr<Tensor_<DataType>> Tens2_in,
                                                                            pair< vector<int>, vector<int> >& ctrs_todo){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::contract_different_tensors_general" <<endl; 
  assert( ctrs_todo.first.size() ==  ctrs_todo.second.size() );

  vector<IndexRange> T1_org_rngs = Tens1_in->indexrange();
  vector<IndexRange> T2_org_rngs = Tens2_in->indexrange();

  shared_ptr<vector<int>> T1_org_order= make_shared<vector<int>>(T1_org_rngs.size());
  iota(T1_org_order->begin(), T1_org_order->end(), 0);

  //note column major ordering
  shared_ptr<vector<int>> T1_new_order = make_shared<vector<int>>( *T1_org_order );
  put_ctrs_at_back( *T1_new_order, ctrs_todo.first);

  shared_ptr<vector<IndexRange>> T1_new_rngs  = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>>        maxs1        = get_num_index_blocks_vec(T1_new_rngs) ;
  shared_ptr<vector<int>>        mins1        = make_shared<vector<int>>(maxs1->size(), 0 );  

  shared_ptr<vector<int>> T2_org_order= make_shared<vector<int>>(T2_org_rngs.size());
  iota( T2_org_order->begin(), T2_org_order->end(), 0);
 
  shared_ptr<vector<int>> T2_new_order = make_shared<vector<int>>( *T2_org_order );
  put_reversed_ctrs_at_front( *T2_new_order, ctrs_todo.second);

  shared_ptr<vector<IndexRange>> T2_new_rngs  = reorder_vector(T2_new_order, T2_org_rngs);
  shared_ptr<vector<int>>        maxs2        = get_num_index_blocks_vec(T2_new_rngs) ;
  shared_ptr<vector<int>>        mins2        = make_shared<vector<int>>(maxs2->size(), 0 );  

  int num_ctrs = ctrs_todo.first.size();

  vector<IndexRange> Tout_unc_rngs(T1_new_rngs->begin(), T1_new_rngs->end()-num_ctrs);
  Tout_unc_rngs.insert(Tout_unc_rngs.end(), T2_new_rngs->begin()+num_ctrs, T2_new_rngs->end());

  shared_ptr<Tensor_<DataType>> Tens_out = make_shared<Tensor_<DataType>>(Tout_unc_rngs);  
  Tens_out->allocate();
  Tens_out->zero();

  //loops over all index blocks of T1 and T2; final index of T1 is same as first index of T2 due to contraction
  shared_ptr<vector<int>> T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, *T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector( T1_new_order, T1_new_rng_blocks); 
    
    size_t ctr_block_size    = 1;
    for ( vector<Index>::reverse_iterator ctr_id_it = T1_new_rng_blocks->rbegin() ; ctr_id_it != T1_new_rng_blocks->rbegin()+num_ctrs; ctr_id_it++) 
       ctr_block_size *= ctr_id_it->size();      

    size_t T1_unc_block_size = get_block_size( T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-num_ctrs); 
    size_t T1_block_size     = get_block_size( T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 

    std::unique_ptr<DataType[]> T1_data_new;
    {
      std::unique_ptr<DataType[]> T1_data_org = Tens1_in->get_block(*T1_org_rng_blocks); 
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }
   
    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    //remember end is one past the end whilst begin really is the beginning
    for ( int qq = 0; qq != num_ctrs+1; qq++ ){
      *(mins2->begin()+qq) = *(maxs1->end()-(qq+1));
      *(maxs2->begin()+qq) = *(maxs1->end()-(qq+1));
      *(T2_rng_block_pos->begin()+qq) = *(T1_rng_block_pos->end()-(qq+1));
    }

    do { 

      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, *T2_new_rngs); 
      shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector( T2_new_order, T2_new_rng_blocks); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin()+num_ctrs, T2_new_rng_blocks->end());

      std::unique_ptr<DataType[]> T2_data_new;   
      {
        std::unique_ptr<DataType[]> T2_data_org = Tens2_in->get_block(*T2_org_rng_blocks);  
        T2_data_new = reorder_tensor_data(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
      }
      
      std::unique_ptr<DataType[]> T_out_data(new DataType[T1_unc_block_size*T2_unc_block_size]);
      std::fill_n(T_out_data.get(), T1_unc_block_size*T2_unc_block_size, 0.0);

      gemm( 'N', 'N', T1_unc_block_size, T2_unc_block_size, ctr_block_size, T1_data_new.get(), T2_data_new.get(), T_out_data.get(), (DataType)1.0, (DataType)0.0  );

      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-num_ctrs);
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin()+num_ctrs, T2_new_rng_blocks->end());
      Tens_out->add_block( T_out_data, T_out_rng_block );

    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, mins2 ));

  } while (fvec_cycle_skipper(T1_rng_block_pos, maxs1, mins1 ));

  
  return Tens_out;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sets all elements of input tensor Tens to the value specified by elem_val 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems(shared_ptr<Tensor_<DataType>> Tens, DataType elem_val  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   cout << "Tensor_Arithmetic::set_tensor_elems all " << endl;
  
   vector<IndexRange> id_ranges = Tens->indexrange();

   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges ) ; 
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
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
//Sets all elements of input tensor Tens to the value specified by elem_val 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( shared_ptr<Tensor_<DataType>> Tens, vector<IndexRange>& id_ranges,
                                                                       DataType elem_val  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   cout << "Tensor_Arithmetic::set_tensor_elems range_block_specific  " << endl;
  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges ) ; 
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
     vector<Index> id_blocks =  *(get_rng_blocks( block_pos, id_ranges ));
     if ( Tens->exists( id_blocks ) ) { 
       unique_ptr<DataType[]> block_data = Tens->get_block(id_blocks);
       std::fill_n( block_data.get(), Tens->get_size(id_blocks), elem_val );
       Tens->put_block(block_data, id_blocks);
     }
   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copies the data of Tens_sub into the appropriate block of Tens_main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( shared_ptr<Tensor_<DataType>> Tens_sub, shared_ptr<Tensor_<DataType>> Tens_main ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::put_sub_tensor range_block_specific  " << endl;
  
   vector<IndexRange> id_ranges = Tens_sub->indexrange();
   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges ) ; 
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
     vector<Index> id_blocks =  *(get_rng_blocks( block_pos, id_ranges ));
     assert( Tens_main->exists( id_blocks ) );
     assert( Tens_sub->exists( id_blocks ) );
     unique_ptr<DataType[]> data_block = Tens_sub->get_block(id_blocks);
     Tens_main->put_block(data_block, id_blocks);
   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copies block of index ranges from Tens1 into the appropriate block of Tens2
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_tensor_range_block( shared_ptr<Tensor_<DataType>> Tens1, shared_ptr<Tensor_<DataType>> Tens2,
                                                                             vector<IndexRange>& id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::put_sub_tensor range_block_specific  " << endl;
  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges ) ; 
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
     vector<Index> id_blocks =  *(get_rng_blocks( block_pos, id_ranges ));
     assert( Tens1->exists( id_blocks ) );
     assert( Tens2->exists( id_blocks ) );
     unique_ptr<DataType[]> data_block = Tens1->get_block(id_blocks);
     Tens2->put_block(data_block, id_blocks);
   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copies block of index ranges from Tens1, reorder, and then put into the appropriate block of Tens2 
// TODO clear up the mess with the arguments; pick either shared_ptrs or references and stick to it.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_reordered_range_block( shared_ptr<Tensor_<DataType>> T1, vector<IndexRange>& id_ranges_T1,
                                                                                shared_ptr<Tensor_<DataType>> T2, vector<IndexRange>& id_ranges_T2,
                                                                                shared_ptr<vector<int>> new_order ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic::put_reordered_range_block range_block_specific  " << endl;
  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( id_ranges_T1 ) ; 
   shared_ptr<vector<int>> block_pos_T1 = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {

     shared_ptr<vector<Index>> id_blocks_T1 = get_rng_blocks( block_pos_T1, id_ranges_T1 );

     shared_ptr<vector<int>> block_pos_T2 = reorder_vector( *new_order, *block_pos_T1 );
     shared_ptr<vector<Index>> id_blocks_T2 = get_rng_blocks( block_pos_T2, id_ranges_T2 );

     assert( T1->exists( *id_blocks_T1 ) );
     assert( T2->exists( *id_blocks_T2 ) );

     unique_ptr<DataType[]> reordered_data_block;
     {
       unique_ptr<DataType[]> data_block = T1->get_block(*id_blocks_T1);
       reordered_data_block =  reorder_tensor_data( data_block.get(), new_order, id_blocks_T1 );
     }
     T2->put_block(reordered_data_block, *id_blocks_T2);

   } while (fvec_cycle_skipper(block_pos_T1, range_lengths, mins ));

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
//TODO Write a modified version of this which erases the old block immediately after the relevant reordered block is created.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( shared_ptr<Tensor_<DataType>> Tens_in, shared_ptr<vector<int>> new_order ){
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
 
     shared_ptr<vector<Index>> orig_id_blocks      = get_rng_blocks( block_pos, *T_id_ranges );

     if ( Tens_in->exists(*orig_id_blocks) ){
       unique_ptr<DataType[]> reordered_data_block;
       {
       unique_ptr<DataType[]> orig_data_block = Tens_in->get_block( *orig_id_blocks );
       reordered_data_block = reorder_tensor_data( orig_data_block.get(), new_order, orig_id_blocks );
       }
       shared_ptr<vector<Index>> reordered_id_blocks = reorder_vector( new_order, orig_id_blocks );
       reordered_block_tensor->put_block( reordered_data_block, *reordered_id_blocks );
     }

   } while ( fvec_cycle_skipper( block_pos, range_lengths, mins ) );

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
  } else if ( num_ids == 8 ) {
    TS.sort_indices_8( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
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
   cout << "Tensor_Arithmetic::get_test_Tensor_row_major" << endl;

   shared_ptr<Tensor_<DataType>> Tens = make_shared<Tensor_<DataType>>(*T_id_ranges);
   Tens->allocate();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> range_lengths = get_range_lengths( T_id_ranges ); 

   DataType put_after_decimal_point  = pow(10 , block_pos->size());
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
   cout << "Tensor_Arithmetic::get_test_Tensor_column_major" << endl;

   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets(*T_id_ranges) ;

   shared_ptr<Tensor_<DataType>> Tens = make_shared<Tensor_<DataType>>(*T_id_ranges);
   Tens->allocate();

   shared_ptr<vector<int>> block_pos  = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> mins       = make_shared<vector<int>>(T_id_ranges->size(),0);  
   shared_ptr<vector<int>> range_maxs = get_range_lengths( T_id_ranges ); 

   DataType put_after_decimal_point  = pow(10 , block_pos->size());
   vector<double> power_10( block_pos->size() );
   for (int ii = 0  ; ii != power_10.size() ; ii++ ) 
     power_10[ii] = pow(10, (double)(power_10.size()-1-ii) );

    // This could be simplified a lot and is slow, however, it is only for testing.
    // Having seperate id_pos and rel_id_pos seperate is the logically simplest way of doing things.
   do {

     vector<Index> T_id_blocks = *(get_rng_blocks( block_pos, *T_id_ranges )); 
     size_t out_size = Tens->get_size(T_id_blocks); 
     unique_ptr<DataType[]> T_block_data( new DataType[out_size] );

     vector<int> id_blocks_sizes(T_id_blocks.size());
     for( int ii = 0 ;  ii != T_id_blocks.size(); ii++)
       id_blocks_sizes[ii] = T_id_blocks[ii].size();

     shared_ptr<vector<int>> id_strides = get_Tens_strides_column_major(id_blocks_sizes);
 
     shared_ptr<vector<int>> id_pos = make_shared<vector<int>>(T_id_ranges->size(),0);  
     for ( int qq = 0 ; qq != block_pos->size(); qq++)
        id_pos->at(qq) =  block_offsets->at(qq).at(block_pos->at(qq));
    
     shared_ptr<vector<int>> id_mins = make_shared<vector<int>>(*id_pos);  
     shared_ptr<vector<int>> id_maxs = make_shared<vector<int>>(T_id_ranges->size());
     for (int xx = 0 ; xx !=id_maxs->size(); xx++ ) 
       id_maxs->at(xx)  =  id_mins->at(xx)+id_blocks_sizes[xx]-1;
        
     shared_ptr<vector<int>> rel_id_pos  = make_shared<vector<int>>(T_id_ranges->size(), 0);  
     shared_ptr<vector<int>> rel_id_mins = make_shared<vector<int>>(T_id_ranges->size(), 0);  
     shared_ptr<vector<int>> rel_id_maxs = make_shared<vector<int>>(T_id_ranges->size());
     for (int xx = 0 ; xx !=rel_id_maxs->size(); xx++ ) 
       rel_id_maxs->at(xx)  =  id_blocks_sizes[xx]-1;

     do {
        T_block_data[inner_product( rel_id_pos->begin(), rel_id_pos->end(), id_strides->begin(), 0 )] = inner_product (id_pos->begin(), id_pos->end(), power_10.begin(), 0 );
        fvec_cycle_skipper(rel_id_pos, rel_id_maxs, rel_id_mins);
     } while (fvec_cycle_skipper( id_pos, id_maxs, id_mins ));

     Tens->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_maxs, mins ));

   return Tens;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_block_of_data( DataType* data_ptr, shared_ptr<vector<IndexRange>> id_ranges, 
                                                                   shared_ptr<vector<int>> block_pos) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::get_block_of_data" << endl;

 // merge this to one loop, but keep until debugged, as this is likely to go wrong...
 // getting id position of id_block ; list of block sizes looks like [n,n,n,n,... , n-1, n-1 n-1] 
 // TODO on closer inspection, the list of block sizes _doesn't_ look like that.
  vector<size_t> id_pos(block_pos->size());
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
 
  unique_ptr<DataType[]> data_block(new DataType[data_block_size])  ;

  copy_n(data_block.get(), data_block_size, data_ptr+data_block_pos);

  return data_block; 
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cannot find Lapack routine for Kronecker products so doing this.  There is probably a better way.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::direct_tensor_product( shared_ptr<Tensor_<DataType>> Tens1,
                                                                       shared_ptr<Tensor_<DataType>> Tens2  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::direct_tensor_product" <<endl; 
  
  shared_ptr<Tensor_<DataType>> Tens_out;

  if ( Tens1->size_alloc() != 1 && Tens2->size_alloc() != 1 ) {

    //note column major ordering
    vector<IndexRange>      T1_rngs = Tens1->indexrange();
    shared_ptr<vector<int>> T1_maxs = get_num_index_blocks_vec( make_shared<vector<IndexRange>>(T1_rngs) );
    shared_ptr<vector<int>> T1_mins = make_shared<vector<int>>(T1_maxs->size(), 0 );  
    shared_ptr<vector<int>> T1_block_pos = make_shared<vector<int>>(T1_maxs->size(), 0 );
    
    vector<IndexRange>      T2_rngs = Tens2->indexrange();
    shared_ptr<vector<int>> T2_maxs = get_num_index_blocks_vec( make_shared<vector<IndexRange>>(T2_rngs) );
    shared_ptr<vector<int>> T2_mins = make_shared<vector<int>>(T2_maxs->size(), 0 );  
    shared_ptr<vector<int>> T2_block_pos = make_shared<vector<int>>(T2_maxs->size(), 0 );
    
    vector<IndexRange> Tout_rngs(T1_rngs.size() + T2_rngs.size());
    copy( T1_rngs.begin(), T1_rngs.end(), Tout_rngs.begin() );
    copy( T2_rngs.begin(), T2_rngs.end(), Tout_rngs.begin()+T1_rngs.size() );
    
    Tens_out = make_shared<Tensor_<DataType>>(Tout_rngs);  
    Tens_out->allocate();
    
    do { 
      
      shared_ptr<vector<Index>> T1_id_blocks = get_rng_blocks( T1_block_pos, T1_rngs); 
      size_t T1_block_size = Tens1->get_size( *T1_id_blocks ); 
    
      std::unique_ptr<DataType[]> T1_data = Tens1->get_block(*T1_id_blocks);  
    
      do { 
    
        shared_ptr<vector<Index>> T2_id_blocks = get_rng_blocks( T2_block_pos, T2_rngs ); 
        size_t T2_block_size = Tens2->get_size( *T2_id_blocks ); 
    
        std::unique_ptr<DataType[]> T2_data = Tens2->get_block( *T2_id_blocks ); 
   
        std::unique_ptr<DataType[]> Tout_data(new DataType[T1_block_size*T2_block_size]);
        
        vector<Index> Tout_id_blocks( T1_id_blocks->size() + T2_id_blocks->size() );
        copy( T1_id_blocks->begin(), T1_id_blocks->end(), Tout_id_blocks.begin() );
        copy( T2_id_blocks->begin(), T2_id_blocks->end(), Tout_id_blocks.begin()+T1_id_blocks->size() );
    
        DataType* T2_data_ptr = T2_data.get();
        DataType* Tout_data_ptr = Tout_data.get();
    
        for ( int qq = 0; qq != T2_block_size ; qq++ ){
          copy_n( T1_data.get(), T1_block_size, Tout_data_ptr );
          scaler( T1_block_size, *T2_data_ptr, Tout_data_ptr ); 
          T2_data_ptr++;
          Tout_data_ptr += T1_block_size;
        }
     
        Tens_out->put_block( Tout_data, Tout_id_blocks );
    
      } while(fvec_cycle_skipper(T2_block_pos, T2_maxs, T2_mins ));
    
    } while (fvec_cycle_skipper(T1_block_pos, T1_maxs, T1_mins ));

  } else  if ( Tens1->size_alloc() == 1 && Tens2->size_alloc() != 1 ) { //silly way of doing things, a stopgap solution

    Tens_out = Tens2->copy(); 
    vector<Index> id_block(1,Tens1->indexrange()[0].range(0));
    unique_ptr<DataType[]> factor_ptr = Tens1->get_block( id_block ); 
    Tens_out->scale(factor_ptr[0]);

  } else if ( (Tens1->size_alloc() != 1 && Tens2->size_alloc() == 1) || (Tens1->size_alloc() == 1 && Tens2->size_alloc() == 1 )  ) {

    Tens_out = Tens1->copy(); 
    vector<Index> id_block(1,Tens2->indexrange()[0].range(0));
    unique_ptr<DataType[]> factor_ptr = Tens2->get_block( id_block ); 
    Tens_out->scale(factor_ptr[0]);

  }   


  
  return Tens_out;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gets a single element of the tensor, inefficient. Try to use other version which takes block_offsets as arg if getting
// more than one element
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType
Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_tensor_element( shared_ptr<Tensor_<DataType>> Tens, vector<int>& id_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::get_tensor_element" ; cout.flush(); print_vector(id_pos, "id_pos"); cout << endl; 

   vector<IndexRange> id_ranges = Tens->indexrange();
   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( id_ranges ) ;
 
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(id_pos.size());  
   vector<int> rel_id_pos(id_pos.size());  
   for ( int ii = 0; ii != id_pos.size(); ii++ ) 
     for ( int jj = 0 ; jj != block_offsets->at(ii).size() ; jj++) 
       if ( (id_pos.at(ii) < block_offsets->at(ii)[jj]) || (jj == block_offsets->at(ii).size()-1) ){ 
         block_pos->at(ii) = jj-1 == -1 ? 0 : jj-1;
         rel_id_pos[ii] == id_pos[ii]-block_offsets->at(ii)[jj];
       } 

   shared_ptr<vector<Index>> id_blocks = get_rng_blocks( block_pos, id_ranges );

   shared_ptr<vector<int>> id_blocks_sizes= make_shared<vector<int>>(id_blocks->size());
   for( int ii = 0 ;  ii != id_blocks->size(); ii++) id_blocks_sizes->at(ii) = id_blocks->at(ii).size();
   shared_ptr<vector<int>> Tens_strides = get_Tens_strides_column_major(*id_blocks_sizes);

   int shift = inner_product(rel_id_pos.begin(), rel_id_pos.end(), Tens_strides->begin(), 0 );

   unique_ptr<DataType[]> T_data_block = Tens->get_block(*id_blocks);
      
   return *(T_data_block.get()+shift);
;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// alpha.[op1(A_ij).op2(B_jl)] + beta.[C_il] = C_il
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<double>::gemm( char op1, char op2, int size_i, int size_l, int size_j, 
                                                         double* A_data, double* B_data, double* C_data, double alpha, double beta ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic<double>::gemm_interface( shared_ptr<Tensor_<double>> Tens, vector<int>& id_pos)" << endl;
  dgemm_( &op1, &op2, size_i, size_l, size_j, alpha, A_data, size_i, B_data, size_j, beta, C_data, size_i);  
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// alpha.[op1(A_ij).op2(B_jl)] + beta.[C_il] = C_il
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<std::complex<double>>::gemm( char op1, char op2, int size_i, int size_l, int size_j, 
                                                                       std::complex<double>* A_data, std::complex<double>* B_data,
                                                                       std::complex<double>* C_data,
                                                                       std::complex<double> alpha, std::complex<double> beta ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic<std::complex<double>>::gemm_interface( shared_ptr<Tensor_<double>> Tens, vector<int>& id_pos)" << endl;
  zgemm_( &op1, &op2, size_i, size_l, size_j, alpha, A_data, size_i,  B_data, size_j, beta, C_data, size_i);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// alpha.[op1(A_ij)B_j] + beta.[C_il] = C_il
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<double>::gemv( char op1, int size_i, int size_j, double* A_data, double* B_data, double* C_data, double alpha, double beta ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic<double>::gemv( shared_ptr<Tensor_<double>> Tens, vector<int>& id_pos)" << endl;
  int stride = 1 ; 
  int int_one = 1 ; 
  dgemv_( &op1, size_i, size_j, alpha, A_data , size_i, B_data, stride, beta, C_data, int_one );  
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// alpha.[op1(A_ij)B_j] + beta.[C_il] = C_il
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<std::complex<double>>::gemv( char op1, int size_i, int size_j, std::complex<double>* A_data,
                                                                       std::complex<double>* B_data, std::complex<double>* C_data, 
                                                                       std::complex<double> alpha, std::complex<double> beta ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic<std::complex<double>>::gemv( shared_ptr<Tensor_<std::complex<double>>> Tens, vector<int>& id_pos)" << endl;
  int stride = 1 ; 
  int int_one = 1 ; 
  zgemv_( &op1, size_i, size_j, alpha, A_data, size_i, B_data, stride, beta, C_data, int_one ); 
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<double>::scaler( int T1_block_size, double T2_data_ptr, double* Tout_data_ptr){  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic::Tensor_Arithmetic<double>::scaler( int T1_block_size, double* T2_data_ptr, double* Tout_data_ptr) " << endl;
  dscal_( T1_block_size, T2_data_ptr, Tout_data_ptr, 1); 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void Tensor_Arithmetic::Tensor_Arithmetic<std::complex<double>>::scaler( int T1_block_size, std::complex<double> T2_data_ptr, std::complex<double>* Tout_data_ptr){  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::Tensor_Arithmetic<double>::scaler( int T1_block_size, complex<double>* T2_data_ptr, complex<double>* Tout_data_ptr) " << endl;
  zscal_( T1_block_size, T2_data_ptr, Tout_data_ptr, 1); 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
double
Tensor_Arithmetic::Tensor_Arithmetic<double>::dot_product( size_t vec_size, double* v1, double* v2){  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::Tensor_Arithmetic<double>::dot_product " << endl;
  return ddot_( vec_size, v1, 1, v2, 1 ); 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
std::complex<double>
Tensor_Arithmetic::Tensor_Arithmetic<std::complex<double>>::dot_product( size_t vec_size, std::complex<double>* v1, std::complex<double>* v2){  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Tensor_Arithmetic::Tensor_Arithmetic<double>::dot_product " << endl;
  return zdotc_( vec_size, v1, 1, v2, 1 ); 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Tensor_Arithmetic::Tensor_Arithmetic<double>;
template class Tensor_Arithmetic::Tensor_Arithmetic<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
