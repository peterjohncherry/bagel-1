#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/debugging_utils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;
using namespace Debugging_Utils;
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elems( unique_ptr<double[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems double" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ) };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   double* ptr = some_data.get();
   double total = 0.0;
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += *ptr ;

   return total; 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
complex<double> Tensor_Arithmetic_Utils::sum_elems( unique_ptr<complex<double>[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems complex<double>" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ) };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   complex<double>* ptr = some_data.get();
   complex<double> total = (complex<double>)(0.0);
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += *ptr ;

   return total; 
}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elem_norms( unique_ptr<double[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems double" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ) };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   return sum_elems( some_data, length ); 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elem_norms( unique_ptr<complex<double>[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems complex<double>" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ) };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   complex<double>* ptr = some_data.get();
   double total = 0.0;
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += abs(*ptr) ;

   return total; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<complex<double>>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<complex<double>>> Tens_in, vector<string>& range_names,
                                                                               shared_ptr<map< string, shared_ptr<IndexRange> > > range_conversion_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 3arg" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> id_ranges(range_names.size());
  for ( int ii = 0; ii != range_names.size() ; ii++ )
    id_ranges[ii] = *range_conversion_map->at(range_names[ii]);

  shared_ptr<Tensor_<complex<double>>> Tens_out = make_shared<Tensor_<complex<double>>>(id_ranges);
  Tens_out->allocate();

  vector<int> range_maxs  = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<complex<double>[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<complex<double>>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<complex<double>>> Tens_in, vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 2arg" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<complex<double>>> Tens_out = make_shared<Tensor_<complex<double>>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<complex<double>[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<double>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<double>> Tens_in, vector<string>& range_names,
                                                                      shared_ptr<map< string, shared_ptr<IndexRange> > > range_conversion_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 3arg" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> id_ranges(range_names.size());
  for ( int ii = 0; ii != range_names.size() ; ii++ )
    id_ranges[ii] = *range_conversion_map->at(range_names[ii]);

  shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(id_ranges);
  Tens_out->allocate();

  vector<int> range_maxs = get_range_lengths( id_ranges );
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<double>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<double>> Tens_in, vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 2arg" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs  = get_range_lengths( id_ranges );
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor( shared_ptr<Tensor_<double>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
   cout << "Tensor_Arithmetic_Utils::Print_Tensor " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
  
   if ( Bagel_id_ranges.size() == 1 ) { //lazy switch; should have generic routine 

     assert( Bagel_id_ranges[0].range(0).size() > 0 ); 

     if ( name == "ID" || ( ( Bagel_id_ranges[0].range().size() == 1 ) && ( Bagel_id_ranges[0].range(0).size() == 1 ) ) ) { 
       cout << "Bagel_id_ranges[0].range(0).size() == 1  " << endl; 
       unique_ptr<double[]> data = Tens->get_block( vector<Index> ( 1, Bagel_id_ranges[0].range(0) ) );
       cout << *data.get()  << endl;

     } else { 

       Print_Vector_Tensor_Format( Tens, name );

     }

   } else {

     shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( Bagel_id_ranges) ;
     
     vector<int> range_lengths  = get_range_lengths( Bagel_id_ranges );
     vector<int> block_pos(range_lengths.size(),0);  
     vector<int> mins(range_lengths.size(),0);  
      
     do {
     
      
       vector<Index> id_blocks = get_rng_blocks( block_pos, Bagel_id_ranges );
      
       if ( !(Tens->exists(id_blocks)) ) { 
     
         cout << "WARNING : TENSOR BLOCK AT " ; cout.flush() ; print_vector( block_pos ) ; cout << " IS NOT STORED FOR PRINTING!! " << endl;
         cout << "NOTE : THIS IS OK IF THERE IS SYMMETRY, BUT OTHERWISE A BIG PROBLEM!! " << endl;
          
       } else {
     
     
         shared_ptr<vector<int>> id_blocks_sizes= make_shared<vector<int>>(id_blocks.size());
         for( int ii = 0 ;  ii != id_blocks.size(); ii++)
           id_blocks_sizes->at(ii) = id_blocks[ii].size();
        
         int id_block_size = Tens->get_size(id_blocks);
       
         unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
         vector<int> Tens_strides = get_Tens_strides_column_major(*id_blocks_sizes);
         double* pos = T_data_block.get();
         int shift = 0 ; 
         if ( block_pos.size() >2 ) { 
           do {
     
             int id_pos_tmp = shift;
             vector<int> id_pos_rel(Tens_strides.size(), 0 );
             for ( int kk = Tens_strides.size()-1 ; kk != 1 ; kk--){
               id_pos_rel[kk] = id_pos_tmp/Tens_strides.at(kk);
               id_pos_tmp -= id_pos_rel[kk]*Tens_strides.at(kk); 
             }

             vector<int> id_pos(block_pos.size());
             for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
               id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii)) + id_pos_rel[ii];

             print_vector( id_pos_rel, "id_pos_rel" ) ; cout << endl;
             print_vector( id_pos, "id_pos" ) ; cout << endl;
     
             for (int jj = 0 ; jj != id_blocks_sizes->at(0) ; jj++ ){
               for (int kk = 0 ; kk != id_blocks_sizes->at(1) ; kk++ ){
                 cout << T_data_block[shift+kk*id_blocks_sizes->at(0)+jj] << " " ; 
               }
               cout <<endl;
             }
     
             shift +=  id_blocks_sizes->at(1)*id_blocks_sizes->at(0);
     
           } while (shift < id_block_size );
     
         } else {
         
           vector<int> id_pos_rel(Tens_strides.size(), 0 );

           vector<int> id_pos(block_pos.size());
           for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
             id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii)) + id_pos_rel[ii];
           
           print_vector( id_pos, "id_pos" ) ; cout << endl;
    
           for (int jj = 0 ; jj != id_blocks_sizes->at(0) ; jj++ ){
             for (int kk = 0 ; kk != id_blocks_sizes->at(1) ; kk++ ){
               cout << T_data_block[shift+kk*id_blocks_sizes->at(0)+jj] << " " ; 
             }
             cout <<endl;
           }
     
           shift += id_blocks_sizes->at(1)*id_blocks_sizes->at(0);
           cout << endl;
         }
       } 
     
     } while ( fvec_cycle_skipper_f2b( block_pos, range_lengths, mins ));
  
   }  

   return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor( shared_ptr<Tensor_<complex<double>>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error( " have not implemented print routines for complex tensors yet!! Aborting!! " ) ;
  return ;
}   

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format( shared_ptr<Tensor_<double>> VecIn, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format" <<endl; 
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  IndexRange VecIn_Idrng = VecIn->indexrange()[0];
  
  for (int ii = 0 ; ii != VecIn_Idrng.range().size(); ii++ ){ 

    vector<Index> id_block = { VecIn_Idrng.range(ii) } ;
    unique_ptr<double[]> VecIn_data_block = VecIn->get_block( id_block );
    double* data_ptr = VecIn_data_block.get();

    for (int jj = 0 ; jj != id_block[0].size() ; jj++ ){ 
      cout << *data_ptr++ << " " ; 
    }   

    cout.flush();
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format( shared_ptr<Tensor_<complex<double>>> VecIn, string name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error("have not implemented print routines for complex tensors yet!! Aborting !! " ); 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor_row_major( shared_ptr<Tensor_<double>> Tens, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
   cout << "Tensor_Arithmetic_Utils::Print_Tensor " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( Bagel_id_ranges) ;

   vector<int> range_lengths = get_range_lengths( Bagel_id_ranges ) ;
   vector<int> block_pos (range_lengths.size(),0);  
   vector<int> mins(range_lengths.size(),0);  
    
   do {

     vector<int> id_pos(block_pos.size());
     for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
       id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii));
    
     vector<Index> id_blocks = get_rng_blocks( block_pos, Bagel_id_ranges );
    
     if ( Tens->exists(id_blocks) ) { 

       vector<int> id_blocks_sizes(id_blocks.size());

       for( int ii = 0 ;  ii != id_blocks.size(); ii++)
         id_blocks_sizes[ii] = id_blocks[ii].size();
       
       int id_block_size = accumulate(id_blocks_sizes.begin(), id_blocks_sizes.end(), 1 , std::multiplies<int>());


       unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
       vector<int> Tens_strides = get_Tens_strides(id_blocks_sizes);

       for (int ii = 0 ; ii != id_blocks_sizes.size() ; ii++ ) 
         id_blocks_sizes[ii]--;       

       vector<int> rel_id_pos(id_blocks_sizes.size(),0);  
    
       double* pos = T_data_block.get();

       vector<int> rel_id_pos_old(id_blocks_sizes.size()-1, 0);

       do {

          if ( *(rel_id_pos.end()-2) != rel_id_pos_old.back() ){
            rel_id_pos_old.back() = *(rel_id_pos.end()-2);
            cout << endl; 
          }
          
          bool newmat = false;
          for ( int ii = 0; ii != rel_id_pos.size()-2 ; ii++  ){
            if ( rel_id_pos.at(ii) != rel_id_pos_old[ii] ) {
              newmat= true;
            }
          }
          if (newmat ){
            rel_id_pos_old = vector<int>(rel_id_pos.begin(), rel_id_pos.end()-1);
            cout << endl <<  "[ " ; cout.flush();
            for ( int kk = 0 ; kk != rel_id_pos.size() ; kk++ ) { 
              cout << rel_id_pos.at(kk) + block_offsets->at(kk).at(block_pos.at(kk)) << " " ; cout.flush();
            }
            cout << "] " << endl; 
          }   

          cout << *pos++ << " " ; cout.flush();        
 
       } while (fvec_cycle_skipper(rel_id_pos, id_blocks_sizes, mins ) ) ;
 
     } else { 
       
       cout << "WARNING : TENSOR BLOCK AT " ; cout.flush() ; print_vector( block_pos ) ; cout << " IS NOT STORED FOR PRINTING!! " << endl;
       cout << "NOTE : THIS IS OK IF THERE IS SYMMETRY, BUT OTHERWISE A BIG PROBLEM!! " << endl;
        
     }

   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));
 
   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor_row_major( shared_ptr<Tensor_<complex<double>>> Tens, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error("have not implemented print routines for complex tensors yet!! Aborting !! " ); 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_Tens_strides( const vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::get_Tens_strides row major" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<int> Tens_strides(range_sizes.size());
  Tens_strides.front() = 1;
  for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
    Tens_strides.at(ii) = Tens_strides.at(ii-1) * range_sizes[ii-1];
   
  return Tens_strides;
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_Tens_strides_column_major( const vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::get_Tens_strides_column_major " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<int> Tens_strides(range_sizes.size());
  Tens_strides.front() = 1;
  if (range_sizes.size() > 1 ) {
    for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
      Tens_strides.at(ii) = Tens_strides.at(ii-1) * range_sizes[ii-1];
  }
  return Tens_strides;
 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at back
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::put_ctr_at_back( const vector<int>& orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctr_at_back" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> new_pos = orig_pos;
  if (ctr_pos == orig_pos.back()){
   return new_pos;
 
  } else {
   
   int ii =0;
   for (  ; ii != orig_pos.size(); ii++)
     if ( orig_pos.at(ii) == ctr_pos ) 
       break; 
 
   new_pos.erase(new_pos.begin()+ii);
   new_pos.push_back(ctr_pos); 
 
   if (new_pos.size() != orig_pos.size())
     throw runtime_error("something has gone wrong in index reshuffling");
 
   return new_pos;
 
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at front
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::put_ctr_at_front( const vector<int>& orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::put_ctr_at_front" <<  endl;  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (ctr_pos == orig_pos.front()){

    vector<int> new_pos = orig_pos; 
    return new_pos;
 
  } else { 
    vector<int> tmp_pos = orig_pos;
    int ii =0;
    for (  ; ii != orig_pos.size(); ii++)
      if ( tmp_pos[ii] == ctr_pos ) 
        break; 
    
    tmp_pos.erase(tmp_pos.begin()+ii);
    vector<int> new_pos(1, ctr_pos);
    new_pos.insert(new_pos.end(), tmp_pos.begin(), tmp_pos.end() );
   
    assert(new_pos.size() == orig_pos.size());

    return new_pos;
  }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Index>
Tensor_Arithmetic_Utils::get_rng_blocks( const vector<int>& block_pos, const vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "TensOp_Computer::get_rng_blocks " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<Index> id_blocks(block_pos.size());
  for( int ii =0; ii != block_pos.size(); ii++ ) 
     id_blocks[ii] = id_ranges[ii].range(block_pos[ii]);
  return id_blocks;
}
////////////////////////////////////////////////////////////////////////////////////////
vector<int>
Tensor_Arithmetic_Utils::get_num_index_blocks_vec(vector<IndexRange>& rngvec) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "TensOp_Computer::get_num_index_blocks" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec.size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( IndexRange& id_rng : rngvec ){
     *iter = id_rng.range().size()-1;
      ++iter;
  }

  return num_id_blocks_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>>
Tensor_Arithmetic_Utils::get_sizes(const vector<Index>& Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_sizes" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  shared_ptr<vector<size_t>> size_vec = make_shared<vector<size_t>>(Idvec.size()); 
  for( int ii = 0 ; ii != Idvec.size() ; ii++ ) 
     (*size_vec)[ii] = Idvec[ii].size();
  return size_vec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t
Tensor_Arithmetic_Utils::get_unc_block_size( vector<Index>& idvec, pair<int,int> ctr ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_unc_block_size" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  size_t block_size = 1; 
  for( int ii = 0 ; ii != idvec.size() ; ii++) 
    if ( (ii != ctr.first) && (ii !=ctr.second))
      block_size *= idvec[ii].size();
  
  return  block_size;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Tensor_Arithmetic_Utils::get_block_size( vector<Index>::iterator startpos, vector<Index>::iterator endpos   ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_block_size" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( vector<Index>::iterator id_it = startpos ; id_it!=endpos; id_it++ ){ 
    block_size *= id_it->size();
  }
  return  block_size;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_range_lengths(const vector<IndexRange>& indexranges ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_range_lengths" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> range_lengths(indexranges.size() ); 
  vector<IndexRange>::const_iterator ir_it = indexranges.begin();
  for (vector<int>::iterator rl_it = range_lengths.begin(); rl_it!= range_lengths.end(); rl_it++, ir_it++ )
    *rl_it = ir_it->range().size()-1; 

  return range_lengths;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Checks all the index ranges in idx_pos have the same length and number of index blocks
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::check_contracted_indexes( std::vector<IndexRange>&  idx_block, std::vector<int>& idx_pos ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::check_contracted_indexes " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  size_t block_len = idx_block[idx_pos.front()].size();
  for ( int ii = 1; ii != idx_pos.size() ; ii++ ) 
    if (block_len != idx_block[idx_pos[ii]].size()){
      cout << "ranges of contracted indexes " << idx_pos[ii] << " and " <<  idx_pos[ii-1]  << "have different lengths" << endl;
      assert(block_len == idx_block[idx_pos[ii]].size());
    }

  block_len = idx_block[idx_pos.front()].range().size();
  for ( int ii = 1; ii != idx_pos.size() ; ii++ ) 
    if (block_len != idx_block[idx_pos[ii]].range().size()){
      cout << "ranges of contracted indexes " << idx_pos[ii] << " and " <<  idx_pos[ii-1]  << "have same lengths but different number of index blocks";
      assert(block_len == idx_block[idx_pos[ii]].range().size());
     }

  return;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Outputs the interval [start_orb, end_orb) of the index blocks at block_pos
// Note the interval is relative to the input IndexRanges.
////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<size_t, size_t>>>
Tensor_Arithmetic_Utils::get_block_start( shared_ptr<vector<IndexRange>> id_ranges, 
                                          shared_ptr<vector<int>> block_pos         ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "TensOp_Arithmetic_Utils::get_block_start" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<pair<size_t,size_t>> block_start_end(block_pos->size());
  for (int ii = 0 ; ii != block_start_end.size() ; ii++){
    size_t block_start = 0;

    for (int jj = 0 ; jj != block_pos->at(ii); jj++) 
      block_start += id_ranges->at(ii).range(jj).size();

    block_start_end[ii] = make_pair(block_start, block_start+id_ranges->at(ii).range(block_pos->at(ii)).size());

  }

  return make_shared<vector<pair<size_t,size_t>>>(block_start_end);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<vector<int>>> Tensor_Arithmetic_Utils::get_block_offsets( const vector<IndexRange>&  ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
 cout << "Tensor_Arithmetic_Utils::get_block_offsets" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<vector<vector<int>>> block_offsets = make_shared<vector<vector<int>>>(ranges.size());
   
  for ( int  ii = 0 ; ii != ranges.size() ; ii++ ) { 
    vector<int> block_offset(ranges[ii].range().size()); 
    block_offset[0] = 0;
    for ( int  jj = 1 ; jj != ranges[ii].range().size() ; jj++ ) { 
      block_offset[jj] = block_offset[jj-1]+ranges[ii].range(jj).size(); 
    }    
    block_offsets->at(ii) =  block_offset;
  }
  
  return block_offsets ; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front, reverse vector version. Note, this will change the input vector orig_pos. Do not change ctrs_pos. 
// The reversal so that when the indices are flattened matrix multiplication can be used.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_reversed_ctrs_at_front( vector<int>& ids_pos, vector<int> ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_reversed_ctrs_at_front" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }

  reverse(ctrs_pos.begin(), ctrs_pos.end());
  ids_pos.insert(ids_pos.begin(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_front( vector<int>& ids_pos, vector<int>& ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_front" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }
  ids_pos.insert(ids_pos.begin(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_front( vector<int>& ids_pos, pair<int,int>& ctr_pair) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_front pair" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   vector<int> ctrs_pos = {ctr_pair.first, ctr_pair.second};

   put_ctrs_at_front(ids_pos, ctrs_pos );
 
   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at back,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_back( vector<int>& ids_pos, vector<int>& ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_back" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }

  ids_pos.insert(ids_pos.end(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
