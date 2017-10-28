#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_computer.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::Print_Tensor( shared_ptr<Tensor_<double>> Tens ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::Print_Tensor " << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();

   shared_ptr<vector<int>> range_lengths  = make_shared<vector<int>>(Tens->indexrange().size()); 
   for (int ii = 0 ; ii != range_lengths->size(); ii++)
      range_lengths->at(ii) =  Tens->indexrange()[ii].range().size()-1;
    
   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   vector<vector<int>> all_id_block_pos(0);
   vector<vector<int>> all_id_block_sizes(0);
    
   for (int ii = 0 ; ii != Tens->indexrange().size(); ii++) {

     vector<int> id_block_pos_part(Tens->indexrange()[ii].range().size());
     vector<int> id_block_sizes_part(Tens->indexrange()[ii].range().size());
     id_block_pos_part[0] = 0;
     id_block_sizes_part[0] = Tens->indexrange()[ii].range(0).size() ;

     for (int jj = 1 ; jj != Tens->indexrange()[ii].range().size(); jj++){ 
        id_block_pos_part[jj]   = id_block_pos_part[jj-1] + Tens->indexrange()[ii].range(jj-1).size();
        id_block_sizes_part[jj] = Tens->indexrange()[ii].range(jj).size();
     }

     all_id_block_pos.push_back(id_block_pos_part);
     all_id_block_sizes.push_back(id_block_sizes_part);
   } 

   cout << "Tdata_block : " << endl; 
   do {

     vector<int> block_id_start(block_pos->size());
     for ( int ii = 0 ; ii != block_id_start.size(); ii++)
       block_id_start[ii] = all_id_block_pos[ii].at(block_pos->at(ii));
    
     vector<Index> id_blocks( block_pos->size() );
     for( int ii = 0 ;  ii != id_blocks.size(); ii++)
       id_blocks[ii] = Tens->indexrange()[ii].range(block_pos->at(ii));
    
     if ( Tens->exists(id_blocks) ) { 
       vector<int> id_blocks_sizes(id_blocks.size());
       shared_ptr<vector<int>> id_blocks_maxs = make_shared<vector<int>>(id_blocks.size());
       for( int ii = 0 ;  ii != id_blocks.size(); ii++){
         id_blocks_sizes[ii] = id_blocks[ii].size();
         id_blocks_maxs->at(ii) = id_blocks[ii].size()-1;
       }
      
       shared_ptr<vector<int>> id_pos = make_shared<vector<int>>(range_lengths->size(),0);
       int id_block_size = accumulate(id_blocks_sizes.begin(), id_blocks_sizes.end(), 1 , std::multiplies<int>());
       
       unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(id_blocks_sizes);
       
       for ( int jj = 0 ; jj != id_block_size ; jj++) {
         cout <<  T_data_block[jj] << " "; cout.flush();
         for ( vector<int>::iterator Titer =  Tens_strides->begin()+1;  Titer != Tens_strides->end(); Titer++ ) {
           if ( (jj > 0) && ( ( (jj+1) % *Titer) == 0) ){ 
             cout << endl;
           }
         }
       }
     }

   } while (fvec_cycle(block_pos, range_lengths, mins ));
 
   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_Tens_strides(vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "get_Tens_strides " << endl;

 shared_ptr<vector<int>> Tens_strides= make_shared<vector<int>>(range_sizes.size());
 Tens_strides->front() = 1;
 for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
   Tens_strides->at(ii) = Tens_strides->at(ii-1) * range_sizes[ii-1];
  
 return Tens_strides;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_CTens_strides( vector<int>& range_sizes, int ctr1 , int ctr2 ) {
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
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_CTens_strides( vector<int>& range_sizes, vector<int>& ctr_idxs_pos ) {
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
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_CTens_strides( shared_ptr<vector<int>> range_sizes, int ctr1 , int ctr2 ) {
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
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::put_ctr_at_back(shared_ptr<vector<int>> orig_pos , int ctr_pos){
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
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::put_ctr_at_front(shared_ptr<vector<int>> orig_pos , int ctr_pos){
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
   
    if (new_pos.size() != orig_pos->size()){
      throw runtime_error("something has gone wrong in index reshuffling");
    } else { 
      return make_shared<vector<int>>(new_pos);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>>
Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> block_pos, vector<IndexRange>& id_ranges) {
////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Equation_Computer::get_rng_blocks " << endl;
  shared_ptr<vector<Index>> id_blocks = make_shared<vector<Index>>(block_pos->size());
  for( int ii =0; ii != block_pos->size(); ii++ ) 
     id_blocks->at(ii) = id_ranges[ii].range(block_pos->at(ii));
  return id_blocks;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>>
Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> block_pos, shared_ptr<vector<shared_ptr<const IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
//cout << "Equation_Computer::get_rng_blocks constver" << endl;

  auto new_ids = make_shared<vector<Index>>(); 
  for( int ii =0; ii != block_pos->size(); ii++ ) {
     new_ids->push_back(old_ids->at(ii)->range(block_pos->at(ii)));
  }
  return new_ids;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<Index>> Equation_Computer::Equation_Computer::get_rng_blocks(shared_ptr<vector<int>> forvec, shared_ptr<vector<shared_ptr<IndexRange>>> old_ids) {
////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Equation_Computer::get_rng_blocks " << endl;

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
vector<int> Equation_Computer::Equation_Computer::get_sizes(vector<Index>& Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
  vector<int> size_vec(Idvec.size()); 
  for( int ii = 0 ; ii != Idvec.size() ; ii++ ) 
     size_vec[ii] = Idvec[ii].size();
  return size_vec;
}

////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>>
Equation_Computer::Equation_Computer::get_sizes(shared_ptr<vector<Index>> Idvec) {
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> Equation_Computer::Equation_Computer::get_range_lengths(shared_ptr<vector<IndexRange>> indexranges ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<int>> range_lengths  = make_shared<vector<int>>(indexranges->size() ); 
  for (int jj = 0 ; jj != indexranges->size() ; jj++ )                                                        
    range_lengths->at(jj) = indexranges->at(jj).range().size()-1; 

  return range_lengths;
}

#ifndef NDEBUG
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Checks all the index ranges in idx_pos have the same length and number of index blocks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::check_contracted_indexes( std::vector<IndexRange>&  idx_block, std::vector<int>& idx_pos ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

#endif
#endif

 // int n1 = norb_;  
 // int n2 = norb_*norb_;  
 // int n3 = n2*norb_;  
 // int n4 = n3*norb_;  
 // int n5 = n4*norb_;  
 //                       
 // cout << endl<< "gamma3 " << endl;
 // for (int i = 0; i != norb_; ++i) 
 //   for (int j = 0; j != norb_; ++j) 
 //     for (int k = 0; k != norb_; ++k) 
 //       for (int l = 0; l != norb_; ++l) 
 //         for (int m = 0; m != norb_; ++m) {
 //           cout << endl << i << " " << j << " "<<  k << " " << l << " " << m << "  ";
 //           for (int n = 0; n != norb_; ++n){
 //             cout << *(gamma3->data()+(i*n5 + j*n4 + k*n3 + l*n2 + m*n1 +n )) << " ";
 //            }
 //         }
 //
 // cout <<endl<<"gamma2 " << endl;
 // for (int k = 0; k != norb_; ++k) 
 //   for (int l = 0; l != norb_; ++l) 
 //     for (int m = 0; m != norb_; ++m) {
 //       cout << endl << k << " " << l << " " << m << "  ";
 //       for (int n = 0; n != norb_; ++n){ 
 //         cout << *(gamma2->data()+ (k*n3 + l*n2 + m*n1 +n))  << " " ;
 //       }                        
 //     }
 //
 // cout << endl << "gamma1 " << endl;
 // for (int m = 0; m != norb_; ++m) {
 //   cout << endl << m << " " ;
 //   for (int n = 0; n != norb_; ++n){ 
 //     cout << *(gamma1->data()+ ( m*n1 +n )) << " "; 
 //   }
 // }
 //  cout << " maxs1 = [ ";   for (int max1 : *maxs1 ) { cout << max1 << " " ; } cout << "]" << endl;
