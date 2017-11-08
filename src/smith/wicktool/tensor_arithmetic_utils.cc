#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/tensor_arithmetic_utils.h>
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::Print_Tensor( shared_ptr<Tensor_<double>> Tens , string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic_Utils::Print_Tensor " << endl;
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( Bagel_id_ranges) ;

   shared_ptr<vector<int>> range_lengths  =  get_range_lengths( Bagel_id_ranges ) ;
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
    
   cout << "Tdata_block : " << endl; 
   do {

     vector<int> id_pos(block_pos->size());
     for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
       id_pos[ii] = block_offsets->at(ii).at(block_pos->at(ii));
    
     vector<Index> id_blocks = *(get_rng_blocks( block_pos, Bagel_id_ranges ));
    
     if ( Tens->exists(id_blocks) ) { 
       vector<int> id_blocks_sizes(id_blocks.size());
       for( int ii = 0 ;  ii != id_blocks.size(); ii++)
         id_blocks_sizes[ii] = id_blocks[ii].size();

       int id_block_size = accumulate(id_blocks_sizes.begin(), id_blocks_sizes.end(), 1 , std::multiplies<int>());
       
       unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(id_blocks_sizes);
       
       cout <<  "[ " ; for ( int elem : id_pos ) { cout << elem << " " ; } cout  << "] " << endl;
       for ( int jj = 0 ; jj != id_block_size ; jj++) {
         cout <<  T_data_block[jj] << " "; cout.flush();
         for ( int kk = 1; kk != Tens_strides->size(); kk++ ) {
           if ( (((jj+1) % Tens_strides->at(kk)) == 0)  && (jj > 0) ){ 
             cout << endl;
             if ( kk > 1 ) {
               for ( int  ll = kk ; ll !=  Tens_strides->size() ; ll++){
                 if ( ( jj+1 % Tens_strides->at(ll) == 0) ) 
                   id_pos[ll]++;
               }
               cout <<  "[ " ; for ( int elem : id_pos ) { cout << elem << " " ; } cout  << "] " << endl;
             }
           }
         }
       }
     }

   } while (fvec_cycle(block_pos, range_lengths, mins ));
 
   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::Print_Tensor_test( shared_ptr<Tensor_<double>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Tensor_Arithmetic_Utils::Print_Tensor_test " << endl;
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( Bagel_id_ranges) ;

   shared_ptr<vector<int>> range_lengths  =  get_range_lengths( Bagel_id_ranges ) ;
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
    
   cout << "Tdata_block : " << endl; 
   do {

     vector<int> id_pos(block_pos->size());
     for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
       id_pos[ii] = block_offsets->at(ii).at(block_pos->at(ii));
    
     print_vector( *block_pos , "block_pos" ); cout <<  endl;
     print_vector( id_pos , "id_pos"); cout<< endl;
     
     vector<Index> id_blocks = *(get_rng_blocks( block_pos, Bagel_id_ranges ));
    
     if ( Tens->exists(id_blocks) ) { 

       vector<int> id_blocks_sizes(id_blocks.size());

       for( int ii = 0 ;  ii != id_blocks.size(); ii++)
         id_blocks_sizes[ii] = id_blocks[ii].size();
       
       int id_block_size = accumulate(id_blocks_sizes.begin(), id_blocks_sizes.end(), 1 , std::multiplies<int>());


       unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(id_blocks_sizes);

       for (int ii = 0 ; ii != id_blocks_sizes.size() ; ii++ ) 
         id_blocks_sizes[ii]--;       

       shared_ptr<vector<int>> block_lengths  = make_shared<vector<int>>( id_blocks_sizes ) ;
       shared_ptr<vector<int>> rel_id_pos     = make_shared<vector<int>>( id_blocks_sizes.size(),0);  
       shared_ptr<vector<int>> rel_id_mins    = make_shared<vector<int>>( id_blocks_sizes.size(),0);  
    
       double* pos = T_data_block.get();

       vector<int> rel_id_pos_old(id_blocks_sizes.size()-1, 0);

       do {


          if ( *(rel_id_pos->end()-2) != rel_id_pos_old.back() ){
            print_vector(rel_id_pos_old , "rel_id_pos_old" ); cout << "   ";  print_vector(*rel_id_pos , "rel_id_pos" );  cout <<  *(rel_id_pos->end()-2)  << "    "  << rel_id_pos_old.back() ;

            rel_id_pos_old.back() = *(rel_id_pos->end()-2);
            cout << endl; 
          }
          
          bool newmat = false;
          for ( int ii = 0; ii != rel_id_pos->size()-2 ; ii++  ){
            if ( rel_id_pos->at(ii) != rel_id_pos_old[ii] ) {
              newmat= true;
            }
          }
          if (newmat ){
            rel_id_pos_old = vector<int>(rel_id_pos->begin(), rel_id_pos->end()-1);
            cout << endl <<  "[ " ; cout.flush();
            for ( int kk = 0 ; kk != rel_id_pos->size() ; kk++ ) { 
              cout << rel_id_pos->at(kk) + block_offsets->at(kk).at(block_pos->at(kk)) << " " ; cout.flush();
            }
            cout << "] " << endl; 
          }   

          cout << *pos++ << " " ; cout.flush();        
 
       } while (fvec_cycle_skipper(rel_id_pos, block_lengths, mins ) ) ;
 
     } else { 
       
       cout << "WARNING : TENSOR BLOCK AT " ; cout.flush() ; print_vector( *block_pos ) ; cout << " IS NOT STORED FOR PRINTING!! " << endl;
        
     }

   } while (fvec_cycle(block_pos, range_lengths, mins ));
 
   return ;
}
               //for ( int  ll = kk ; ll !=  Tens_strides->size() ; ll++){
               //  if ( ( jj+1 % Tens_strides->at(ll) == 0) ) 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_Tens_strides(vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "get_Tens_strides " << endl;

 shared_ptr<vector<int>> Tens_strides= make_shared<vector<int>>(range_sizes.size());
 Tens_strides->front() = 1;
 for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
   Tens_strides->at(ii) = Tens_strides->at(ii-1) * range_sizes[ii-1];
  
 return Tens_strides;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_CTens_strides( vector<int>& range_sizes, int ctr1 , int ctr2 ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "get_CTens_strides " <<  endl;

  shared_ptr<vector<int>>  CTens_strides = make_shared<vector<int>>(range_sizes.size(), 1); 
  for ( int ii = 0  ; ii!= range_sizes.size(); ii++ ) 
    for ( int jj = ii-1  ; jj!= -1; jj-- ) 
      if (jj!= ctr1 && jj!=ctr2) 
        CTens_strides->at(ii) *= range_sizes[jj];
  
  CTens_strides->at(ctr1) = 0;
  CTens_strides->at(ctr2) = 0;

  return CTens_strides;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the vector of strides of the tensor for which ctr_idxs_pos are contracted from the 
//range sizes of the tensor where they are not contracted.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_CTens_strides( vector<int>& range_sizes, vector<int>& ctr_idxs_pos ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "get_CTens_strides " <<  endl;

  shared_ptr<vector<int>>  CTens_strides = make_shared<vector<int>>(range_sizes.size(), 1); 
  vector<bool> is_unc(ctr_idxs_pos.size(), false);
  for ( int pos : ctr_idxs_pos)  
    is_unc[pos] = false;

  for ( int ii = 0  ; ii!= range_sizes.size(); ii++ ) 
    for ( int jj = ii-1  ; jj!= -1; jj-- ) 
      if (is_unc[jj]) 
        CTens_strides->at(ii) *= range_sizes[jj];
  
  for ( int ii : ctr_idxs_pos ) 
    CTens_strides->at(ii) = 0 ; 

  return CTens_strides;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_CTens_strides( shared_ptr<vector<int>> range_sizes, int ctr1 , int ctr2 ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "get_CTens_strides " <<  endl;

  shared_ptr<vector<int>>  CTens_strides = make_shared<vector<int>>(range_sizes->size(), 1); 
  for ( int ii = 0  ; ii!= range_sizes->size(); ii++ ) 
    for ( int jj = ii-1  ; jj!= -1; jj-- ) 
      if (jj!= ctr1 && jj!=ctr2) 
        CTens_strides->at(ii) *= range_sizes->at(jj);
  
  CTens_strides->at(ctr1) = 0;
  CTens_strides->at(ctr2) = 0;

  return CTens_strides;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at back
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::put_ctr_at_back(shared_ptr<vector<int>> orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "put_ctr_at_back" << endl;

  if (ctr_pos == orig_pos->back()){
 
   return make_shared<vector<int>>(*orig_pos);
 
  } else {
   
   vector<int> new_pos = *orig_pos;
   int ii =0;
   for (  ; ii != orig_pos->size(); ii++)
     if ( orig_pos->at(ii) == ctr_pos ) 
       break; 
 
   new_pos.erase(new_pos.begin()+ii);
   new_pos.push_back(ctr_pos); 
 
   if (new_pos.size() != orig_pos->size())
     throw runtime_error("something has gone wrong in index reshuffling");
 
   return make_shared<vector<int>>(new_pos);
 
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at front
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::put_ctr_at_front(shared_ptr<vector<int>> orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "put_ctr_at_front" <<  endl;  
  if (ctr_pos == orig_pos->front()){
 
    return make_shared<vector<int>>(*orig_pos);
 
  } else { 
    vector<int> tmp_pos = *orig_pos;
    int ii =0;
    for (  ; ii != orig_pos->size(); ii++)
      if ( tmp_pos[ii] == ctr_pos ) 
        break; 
    
    tmp_pos.erase(tmp_pos.begin()+ii);
    vector<int> new_pos(1,  ctr_pos);
    new_pos.insert(new_pos.end(), tmp_pos.begin(), tmp_pos.end() );
   
    assert(new_pos.size() == orig_pos->size());

    return make_shared<vector<int>>(new_pos);
  }
  

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>>
Tensor_Arithmetic_Utils::get_rng_blocks(shared_ptr<vector<int>> block_pos, vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Equation_Computer::get_rng_blocks " << endl;
  shared_ptr<vector<Index>> id_blocks = make_shared<vector<Index>>(block_pos->size());
  for( int ii =0; ii != block_pos->size(); ii++ ) 
     id_blocks->at(ii) = id_ranges[ii].range(block_pos->at(ii));
  return id_blocks;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>>
Tensor_Arithmetic_Utils::get_rng_blocks(shared_ptr<vector<int>> block_pos, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "Equation_Computer::get_rng_blocks constver" << endl;

  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != block_pos->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(block_pos->at(ii)));
  }
  return new_ids;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>>
Tensor_Arithmetic_Utils::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Equation_Computer::get_rng_blocks " << endl;

  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != forvec->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(forvec->at(ii)));
  }
  return new_ids;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>>
Tensor_Arithmetic_Utils::get_sizes(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec, int skip_id) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<int>>(); 
  for( int ii =0 ; ii!=rngvec->size();ii++)
    if(ii!=skip_id) 
     size_vec->push_back(rngvec->at(ii)->size());
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>> 
Tensor_Arithmetic_Utils::get_sizes(shared_ptr<vector<Index>> Idvec, int skip_id) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( int ii =0 ; ii!=Idvec->size();ii++)
    if(ii!=skip_id) 
     size_vec->push_back(Idvec->at(ii).size());
  return size_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
vector<int>
Tensor_Arithmetic_Utils::get_num_index_blocks_vec(vector<IndexRange>& rngvec) {
////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec.size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( IndexRange id_rng : rngvec ) 
     *iter++ = id_rng.range().size()-1;

  return num_id_blocks_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>>
Tensor_Arithmetic_Utils::get_num_index_blocks_vec(shared_ptr<vector<IndexRange>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec->size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( IndexRange id_rng : *rngvec ) 
     *iter++ = id_rng.range().size()-1;

  return make_shared<vector<int>>(num_id_blocks_vec);
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>>
Tensor_Arithmetic_Utils::get_num_index_blocks_vec(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec->size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( shared_ptr<const IndexRange> id_rng : *rngvec ) 
     *iter++ = id_rng->range().size()-1;

  return make_shared<vector<int>>(num_id_blocks_vec);
}


////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>>
Tensor_Arithmetic_Utils::get_sizes(shared_ptr<vector<shared_ptr<const IndexRange>>> rngvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<int>>(); 
  for( auto elem : *rngvec ) 
     size_vec->push_back(elem->size());
  return size_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
vector<int>
Tensor_Arithmetic_Utils::get_sizes(vector<Index>& Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  vector<int> size_vec(Idvec.size()); 
  for( int ii = 0 ; ii != Idvec.size() ; ii++ ) 
     size_vec[ii] = Idvec[ii].size();
  return size_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>>
Tensor_Arithmetic_Utils::get_sizes(shared_ptr<vector<Index>> Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  auto size_vec = make_shared<vector<size_t>>(); 
  for( auto elem : *Idvec ) 
     size_vec->push_back(elem.size());
  return size_vec;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t
Tensor_Arithmetic_Utils::get_unc_block_size( vector<Index>& idvec, pair<int,int> ctr ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  
  size_t block_size = 1; 
  for( int ii = 0 ; ii != idvec.size() ; ii++) 
    if ( (ii != ctr.first) && (ii !=ctr.second))
      block_size *= idvec[ii].size();
  
  return  block_size;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Tensor_Arithmetic_Utils::get_block_size( vector<Index>::iterator startpos,
                                                             vector<Index>::iterator endpos   ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( vector<Index>::iterator id_it = startpos ; id_it!=endpos; id_it++ ){ 
    block_size *= id_it->size();
  }
  return  block_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_range_lengths(shared_ptr<vector<IndexRange>> indexranges ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<int>> range_lengths  = make_shared<vector<int>>(indexranges->size() ); 
  for (int jj = 0 ; jj != indexranges->size() ; jj++ )                                                        
    range_lengths->at(jj) = indexranges->at(jj).range().size()-1; 

  return range_lengths;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Tensor_Arithmetic_Utils::get_range_lengths(vector<IndexRange>& indexranges ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<int>> range_lengths  = make_shared<vector<int>>(indexranges.size() ); 
  for (int jj = 0 ; jj != indexranges.size() ; jj++ )                                                        
    range_lengths->at(jj) = indexranges[jj].range().size()-1; 

  return range_lengths;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Checks all the index ranges in idx_pos have the same length and number of index blocks
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::check_contracted_indexes( std::vector<IndexRange>&  idx_block, std::vector<int>& idx_pos ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
//  cout << "Equation_Computer::get_block_info" << endl;

  vector<pair<size_t,size_t>> block_start_end(block_pos->size());
  for (int ii = 0 ; ii != block_start_end.size() ; ii++){
    size_t block_start = 0;

    for (int jj = 0 ; jj != block_pos->at(ii); jj++) 
      block_start += id_ranges->at(ii).range(jj).size();

    block_start_end[ii] = make_pair(block_start, block_start+id_ranges->at(ii).range(block_pos->at(ii)).size());

  }
//  cout << "block_start_end = [ " ; cout.flush(); for ( auto elem : block_start_end  ) {cout << "(" << elem.first <<  ","  << elem.second << ") " ; } cout << " ] " << endl;

  return make_shared<vector<pair<size_t,size_t>>>(block_start_end);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<vector<int>>> Tensor_Arithmetic_Utils::get_block_offsets(vector<IndexRange>&  ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic_Utils::get_block_offsets" << endl;
 
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<vector<int>>> Tensor_Arithmetic_Utils::get_block_offsets(shared_ptr<vector<IndexRange>>  ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Tensor_Arithmetic_Utils::get_block_offsets" << endl;
 
  shared_ptr<vector<vector<int>>> block_offsets = make_shared<vector<vector<int>>>(ranges->size());
   
  for ( int  ii = 0 ; ii != ranges->size() ; ii++ ) { 
    vector<int> block_offset(ranges->at(ii).range().size()); 
    block_offset[0] = 0;
    for ( int  jj = 1 ; jj != ranges->at(ii).range().size() ; jj++ ) { 
      block_offset[jj] = block_offset[jj-1]+ranges->at(ii).range(jj).size(); 
    }    
    block_offsets->at(ii) =  block_offset;
  }
  
  return block_offsets ; 

}

#endif
