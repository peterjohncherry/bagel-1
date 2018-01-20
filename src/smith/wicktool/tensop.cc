#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/tensop.h>

using namespace std;
using namespace WickUtils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_General::TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<const range_block_info >>> all_ranges) : 
                                idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges), orig_factor_(factor), num_idxs_(idxs.size()),
                                unique_range_blocks_(*unique_range_blocks),  all_ranges_(*all_ranges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_General::TensOp_General "<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  unique_range_blocks_ptr_ = make_shared<const vector< shared_ptr< const vector<string>>>>(unique_range_blocks_);

  all_ranges_ptr_ =  make_shared<const map< const vector<string>, shared_ptr<const range_block_info >>>(all_ranges_);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MultiTensOp_General::MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                          std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                                          std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                          std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr< const range_block_info>>> all_ranges ) :
                                          TensOp_General::TensOp_General( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, unique_range_blocks, all_ranges ),
                                          cmlsizevec_(cmlsizevec) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp_General::MultiTensOp_General "<< endl;
           
  cmlsizevec_ptr_ =  make_shared<const vector<int>>(cmlsizevec_);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                  vector<bool>& aops, DataType orig_factor,
                                  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>), int, int > >& symmfuncs, 
                                  vector<bool(*)(shared_ptr<vector<string>>) >& constraints,
                                  string& Tsymm, shared_ptr<StatesInfo<DataType>> target_states ) :
                                  name_(name), spinfree_(true), symmfuncs_(symmfuncs), constraints_(constraints), Tsymm_(Tsymm),
                                  target_states_(target_states)  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp" <<   endl;
          
  CTP_map  = make_shared< map< string, shared_ptr<CtrTensorPart<DataType>> >>();
   
  std::vector<int> plus_ops;
  std::vector<int> kill_ops;
  for ( int ii =0 ; ii != aops.size() ; ii++ ) {
    if(aops[ii]) {
      plus_ops.push_back(ii);
    } else {
      kill_ops.push_back(ii);
    }
  }

  auto test_var = generate_ranges( idxs, idx_ranges, aops);

  pair<double,double> orig_factor_tmp =  make_pair(1.0, 1.0);

  Op_dense_ = make_shared<const TensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, orig_factor_tmp, get<1>(test_var), get<0>(test_var) );

  cout << "TensOp::TensOp" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO this is too slow; just build ranges using symmetry, don't check.
template<typename DataType>
tuple< std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<const range_block_info> >>,
       std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>> >
TensOp::TensOp<DataType>::generate_ranges( vector<string>& idxs,  vector<vector<string>>& idx_ranges, vector<bool>& aops ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::generate_ranges" <<   endl;
  //set up loop utils

  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
  int num_idxs = orig_idxs->size();

  vector< shared_ptr< const vector<string>> > unique_range_blocks;
  map< const vector<string>, shared_ptr<const range_block_info> > all_ranges;
  
  auto apply_symmetry = [ &all_ranges, &orig_idxs, &aops ]( const vector<string>& ranges_1,  const vector<string>& ranges_2 ) {
     const bool a_unique_range_block = true; //TODO put symmetry back in
     if (a_unique_range_block) {
       const pair<int,int> fac(1,1);
       const pair<double,double> fac_new(1.0,1.0);
       shared_ptr<const vector<string>> unique_ranges = make_shared<const vector<string>>(ranges_2); 
       shared_ptr<const vector<string>> trans_idxs = orig_idxs;
       bool survives = WickUtils::RangeCheck( ranges_2, aops );
       all_ranges.emplace(ranges_2, make_shared< const range_block_info >( a_unique_range_block, survives, fac_new, unique_ranges, unique_ranges, orig_idxs, trans_idxs ) );
     }
     return false;
  };

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( num_idxs, 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( num_idxs );
  int  num_cycles = 1;
  bool is_unique = true;
  const pair<int,int> fac(1,1);
  const pair<double,double> fac_new(1.0,1.0);

  for( int ii = 0; ii != idx_ranges.size(); ii++ ){
    maxs->at(ii) = ( idx_ranges[ii].size()-1 );
    num_cycles *= idx_ranges[ii].size();
  }

  //generate all possible ranges
  vector<vector<string>> possible_ranges(0);
  for (int ii = 0 ; ii != num_cycles ; ii++){
    vector<string> new_range(idx_ranges.size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range[jj] = idx_ranges[jj][fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges.push_back(new_range);
    fvec_cycle( fvec, maxs );
  } 

  // Initialize range maps and add in first range to start things off.
  const vector<string> init_range = possible_ranges[0];
  shared_ptr<const vector<string>> orig_range = make_shared<const vector<string>>(possible_ranges[0]);
  shared_ptr<const vector<string>> trans_range = orig_range;
  shared_ptr<const vector<string>> trans_idxs = orig_idxs;
   
  bool survives = WickUtils::RangeCheck( *orig_range, aops );
  all_ranges.emplace(init_range, make_shared< const range_block_info >( is_unique, survives, fac_new, orig_range, trans_range, orig_idxs, trans_idxs ) );
 
  unique_range_blocks = vector< shared_ptr< const vector<string>>>(1, trans_range );
  //Apply symmetry operations to remove unnecessary ranges 
  for (int ii = 1 ; ii!=possible_ranges.size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
      if( apply_symmetry(possible_ranges[kk], possible_ranges[ii]) ){ //note that this sets the elements of all_ranges_
        break;
      }
      if(kk == ii-1){
        unique_range_blocks.push_back(make_shared<const vector<string>>(possible_ranges[ii]));
        break;
      }
    } 
  }

  auto all_ranges_ptr = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<const range_block_info >>> (all_ranges);
  auto unique_ranges_ptr = make_shared<std::vector< std::shared_ptr< const std::vector<std::string>> >>(unique_range_blocks);                                                                                                                
  return tie ( all_ranges_ptr, unique_ranges_ptr );
              
  cout << "TensOp::generate_ranges()" <<   endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
bool TensOp::TensOp<DataType>::satisfies_constraints( vector<string>& ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  shared_ptr<vector<string>> ranges_ptr = make_shared<vector<string>>(ranges);
  for (auto OK : constraints_) 
    if( !OK(ranges_ptr))
      return false;
  return true;
}; 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and ranges; note this is done by TensOp not TensOp gen as the contractions needed will be 
//dependent on the sparsity associated with the relevant state
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
  cout << "TensOp get_ctrs_tens_ranges" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp get_ctrs_tens_ranges" << endl;
  cout << name_ << endl;

  cout << " X1" << endl;  
 //puts uncontracted ranges into map 
 shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector< pair<int,int>>>(0);
 for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
   shared_ptr<vector<pair<int,int>>>  ReIm_factors = make_shared< vector<pair<int,int>>>(1, rng_it->second->factors()); 
   shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
   shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );

   shared_ptr<CtrTensorPart<DataType>>  CTP       = make_shared< CtrTensorPart<DataType> >( full_idxs, full_ranges, noctrs, ReIm_factors ); 
   cout << "CTP->myname() = " << CTP->myname() << " is going into CTP_map" <<  endl;
   CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs
 }
  cout << " X2" << endl;  
 
  //puts_contracted ranges into map
  for ( int nctrs = 1 ; nctrs != (num_idxs()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>  ctr_lists = get_unique_pairs(plus_ops(), kill_ops(), nctrs);
    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
        bool valid =true;
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
            valid = false;
            break;
          }
        }

        if (valid) {
          shared_ptr<vector<pair<int,int>>> ReIm_factors = make_shared<vector<pair<int,int>>>(1, rng_it->second->factors()); 
          shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
          shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );

          shared_ptr<CtrTensorPart<DataType>> CTP = make_shared<CtrTensorPart<DataType>>( full_idxs, full_ranges, ctr_vec, ReIm_factors ); 
          cout << "CTP->myname() = " << CTP->myname() << "is going into CTP_map " << endl;
          CTP_map->emplace(CTP->myname(), CTP); 

        }
      }
    }
  }
  cout << "leaving TensOp::TensOp get_ctrs_tens_ranges" << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
MultiTensOp::MultiTensOp<DataType>::MultiTensOp( std::string name, bool spinfree,
                                                 std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>& orig_tensors,
                                                 shared_ptr<StatesInfo<DataType>> target_states ):
                                                 TensOp::TensOp<DataType>( name, spinfree, target_states ),
                                                 orig_tensors_(orig_tensors), num_tensors_(orig_tensors.size()) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::MultiTensOp" << endl;
    
  vector<string> idxs;
  vector<vector<string>> idx_ranges;
  vector<bool> aops;
 
  vector<int> cmlsizevec(num_tensors_);
  vector<int> plus_ops;
  vector<int> kill_ops;
 
  int num_idxs = 0;
  for ( int ii = 0;  ii != orig_tensors_.size() ; ii++ ) {
  
    cmlsizevec[ii]  = num_idxs;
    idx_ranges.insert(  idx_ranges.end(), orig_tensors_[ii]->idx_ranges()->begin(), orig_tensors_[ii]->idx_ranges()->end() );
    idxs.insert( idxs.end(),orig_tensors_[ii]->idxs()->begin(), orig_tensors_[ii]->idxs()->end() );
    aops.insert( aops.end(),orig_tensors_[ii]->aops()->begin(), orig_tensors_[ii]->aops()->end() );
  
    for ( int elem : *(orig_tensors_[ii]->plus_ops()) )
      plus_ops.push_back(elem + cmlsizevec[ii]);
  
    for ( int elem : *(orig_tensors_[ii]->kill_ops()) )
      kill_ops.push_back(elem + cmlsizevec[ii]);
  
    num_idxs += orig_tensors_[ii]->num_idxs();
  } 
  shared_ptr<const map < const vector<string> , shared_ptr<const range_block_info > > >  all_ranges_ptr = generate_ranges( idxs, aops, cmlsizevec );

  auto  unique_range_blocks = make_shared<std::vector< std::shared_ptr<const std::vector<std::string>>>>(0);
  for ( auto all_range_map_it : *all_ranges_ptr ) 
    unique_range_blocks->push_back( make_shared<const vector<string>>(all_range_map_it.first) ); 
  
  Op_dense_ = make_shared<const MultiTensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, make_pair(1.0,1.0), cmlsizevec, unique_range_blocks, all_ranges_ptr ); 
  CTP_map  = make_shared< map< string, shared_ptr<CtrTensorPart<DataType>> >>();
  CMTP_map = make_shared< map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>();
   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr< const map< const vector<string>, shared_ptr<const range_block_info> >>
MultiTensOp::MultiTensOp<DataType>::generate_ranges( vector<string>& idxs, vector<bool>& aops, vector<int>& cmlsizevec ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;

  int num_idxs = idxs.size();
  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
  vector< map< const  vector<string>, shared_ptr< const range_block_info >>::const_iterator> rng_maps(num_tensors_);  
  vector<int> posvec( num_tensors_, 0 ); 

  shared_ptr< const map <const vector<string>, shared_ptr<const range_block_info >>>  all_ranges_tmp_ptr;   
  cout << " num_tensors_ = " << num_tensors_ << endl;

  // Such an ugly loop that this should be fixed when switch ranges to numbers
  // The reason for it is that maps can't map->begin()+ii , so must iterate
  // This is why you had the ugly thing before; it was much faster than this stupid thing with all it's checking.
  // Maybe add some comments in future so you can understand things...
  if ( num_tensors_ > 1 ) { 
    map< const vector<string>, shared_ptr<const range_block_info> > all_ranges;

    shared_ptr<vector<int>> forvec = make_shared<vector<int>>(num_tensors_, 0 ); 
    shared_ptr<vector<int>> old_forvec = make_shared<vector<int>>(*forvec);

    shared_ptr<vector<int>> mins   = make_shared<vector<int>>(num_tensors_, 0 );  
    shared_ptr<vector<int>> maxs   = make_shared<vector<int>>(num_tensors_, 0 );  
    for ( int ii = 0 ; ii != orig_tensors_.size() ; ii++ ){ 
      maxs->at(ii) = orig_tensors_[ii]->all_ranges()->size()-1;
      rng_maps[ii] = orig_tensors_[ii]->all_ranges()->begin() ; 
    }

    do {

      for ( int ii = 0 ; ii != forvec->size() ; ii++ ) {
        if ( old_forvec->at(ii) != forvec->at(ii) ) {
          if ( forvec->at(ii) == 0 ) {
            rng_maps[ii] = orig_tensors_[ii]->all_ranges()->begin() ;
          } else {
            rng_maps[ii]++;
          }
        }
      }
      old_forvec = make_shared<vector<int>>(*forvec);

      int  Re_factor = 1;
      int  Im_factor = 1;
      
      vector<string> merged_oranges(num_idxs);
      vector<string> merged_uranges(num_idxs);
      vector<string> merged_uqidxs(num_idxs);
      vector<string> merged_org_idxs(num_idxs);
       
      auto split_block = make_shared< vector <shared_ptr<const range_block_info >>>(0);
      for (int jj = 0 ; jj != num_tensors_ ; jj++){
       
        shared_ptr<const range_block_info> block = rng_maps[jj]->second;

        //merged ranges
        copy( rng_maps[jj]->first.begin(), rng_maps[jj]->first.end(), merged_oranges.begin() + cmlsizevec[jj] );             
        copy( block->unique_block()->begin(), block->unique_block()->end(), merged_uranges.begin() + cmlsizevec[jj] ); 
        copy( block->transformed_idxs()->begin(), block->transformed_idxs()->end(), merged_uqidxs.begin() + cmlsizevec[jj] );     
      
        int Re_factor_buff = Re_factor;  
        int Im_factor_buff = Im_factor; //TODO only need one buff, but I'm very tired so could be wrong. fix later.
        Re_factor = Re_factor_buff * block->Re_factor() - Im_factor_buff * block->Im_factor() ;
        Im_factor = Re_factor_buff * block->Im_factor() + Im_factor_buff * block->Re_factor();
      
      }

      //TODO Must obtain from constraint functions 
      bool survives = WickUtils::RangeCheck( merged_oranges, aops ) ;
      pair<double,double> combined_factor_new = make_pair(Re_factor, Im_factor);
      bool merged_unique = ( merged_oranges == merged_uranges ) ?  false : true ;
      all_ranges.emplace(merged_oranges, make_shared<const range_block_info>( merged_unique, survives, combined_factor_new,
                                                                              make_shared<const vector<string>>(merged_oranges), make_shared<const vector<string>>(merged_uranges),
                                                                              orig_idxs, make_shared<const vector<string>>(merged_uqidxs) ));

      split_ranges_.emplace( merged_oranges, split_block );
 
    } while( fvec_cycle_skipper( forvec, maxs, mins ) );
  
    all_ranges_tmp_ptr =  make_shared< const map < const vector<string> , shared_ptr<const range_block_info > > >( all_ranges );

  } else { 

    all_ranges_tmp_ptr =  orig_tensors_[0]->all_ranges();
   
    for ( auto  elem : *all_ranges_tmp_ptr ) {
      auto tmp = make_shared<vector<shared_ptr<const range_block_info >>>(1, elem.second);
      split_ranges_.emplace( elem.first , tmp ); 
    }
   
  }
  cout << "leaving MultiTensOp::generate_ranges()" << endl;

  return all_ranges_tmp_ptr ;
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and range, looks totally inefficient, but generating everything and then checking in sequence
//seems to work out faster; it's either check everything lots of times, or regenerate it lots of times. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp get_ctrs_tens_ranges" << endl;
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp::get_ctrs_tens_ranges " <<  endl;

  //puts uncontracted ranges into map 
  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);
 
  //silly, should just test {act,act,....} against contraints, and check if act in each ranges... 
  for (auto rng_it = Op_dense_->all_ranges()->begin(); rng_it != Op_dense_->all_ranges()->end(); rng_it++) {
    bool check = true;
    for (int xx = 0; xx !=rng_it->first.size() ; xx++ ) {
      if (rng_it->first[xx][0] != 'a') {
        check=false;
        break;
      }
    }
    if(!check){
      continue;
    } else {
      
      enter_into_CMTP_map(*noctrs, rng_it->second->factors(), rng_it->first );
    }
  }

  //puts_contractions, with specified ranges, into the map
  for ( int nctrs = 1 ; nctrs != (Op_dense_->num_idxs()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> ctr_lists = get_unique_pairs( Op_dense_->plus_ops(), Op_dense_->kill_ops(), nctrs );

    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = Op_dense_->all_ranges()->begin(); rng_it != Op_dense_->all_ranges()->end(); rng_it++) {
  
        bool valid = true;
        //checks ranges for contractions match
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
            valid = false;
            break;
          }
        }
        if (valid)
          enter_into_CMTP_map(*ctr_vec, rng_it->second->factors(), rng_it->first );
        
      }
    }
  }
  return;
}  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::enter_into_CMTP_map(pint_vec ctr_pos_list, pair<int,int> ReIm_factors, const vector<string>& id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::enter_into_CMTP_map" << endl;
#endif 
////////////////////////////////////////////////////////////////////////////////////////

//cout << "MultiTensOp::enter_into_CMTP_map" << endl;
  shared_ptr<vector<shared_ptr<CtrTensorPart<DataType>>>> CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart<DataType>> >> (num_tensors_); 
  vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
  vector<vector<pair<int,int>>> sameT_ctrs_pos( num_tensors_,  pint_vec(0));
  shared_ptr<vector<pair<int,int>>> no_ctrs =  make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> ReIm_factor_vec = make_shared<vector<pair<int,int>>>(1, ReIm_factors ) ; 
  //seperate contractions into those on the same tensor, and those between different tensors 
  // TODO tidy this up, e.g., use lambda to get cross pos
  for ( pair<int,int> ctr_pos : ctr_pos_list ) {
 
    pair<int,int> ctr1;
    pair<int,int> ctr2;

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.first -= orig_tensors_[ii]->num_idxs()) < 0 ){
         ctr1 = make_pair(ii, ctr_pos.first + orig_tensors_[ii]->num_idxs());
        break;
      }

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.second -= orig_tensors_[ii]->num_idxs()) < 0 ){
         ctr2 = make_pair(ii, ctr_pos.second + orig_tensors_[ii]->num_idxs());
        break;
      } 

    if (ctr1.first == ctr2.first) {
      sameT_ctrs_pos.at(ctr1.first).push_back(make_pair(ctr1.second, ctr2.second));  
    } else {
      diffT_ctrs_pos.push_back(make_pair(ctr1,ctr2));
    }
  } 

  //get_ranges for individual tensors
  for (int ii = 0 ; ii !=num_tensors_; ii++ ){ 
    shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(id_ranges.begin()+Op_dense_->cmlsizevec(ii), id_ranges.begin()+cmlsizevec(ii)+orig_tensors_[ii]->num_idxs());
    shared_ptr<vector<string>>  TS_idxs = make_shared<vector<string>>( *(orig_tensors_[ii]->idxs()) );
   
    // TODO Fix constructor so it will not add if CTP is sparse block
    if( sameT_ctrs_pos[ii].size() != 0 ) {
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, make_shared<vector<pair<int,int>>>(sameT_ctrs_pos[ii]), ReIm_factor_vec ) ; 
    } else { 
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, no_ctrs, ReIm_factor_vec ); 
    }
    CTP_map->emplace(CTP_vec->at(ii)->name, CTP_vec->at(ii)); //TODO silly hack way of dealing with states find a better way 
  }
  

  //TODO silly hack way of dealing with states find a better way
  shared_ptr<CtrMultiTensorPart<DataType>> CMTP = make_shared<CtrMultiTensorPart<DataType> >(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos) ); 
  CMTP_map->emplace(CMTP->myname(), CMTP);
  CTP_map->emplace(CMTP->myname(), CMTP);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Only build the CTP for the unique ranges; these are the only ones needed. If gamma generator
// asks for any other ones, then you have a bug.
///////////////////////////.///////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DType>
void TensOp::TensOp<DType>::get_ctp_idxs_ranges( shared_ptr<vector<pair<int,int>>> ctrs_pos, shared_ptr<const range_block_info> block_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<bool> get_unc(block_info->orig_idxs()->size(), true);
  for (int ii =0; ii != ctrs_pos->size() ; ii++){
    get_unc[ctrs_pos->at(ii).first] = false;
    get_unc[ctrs_pos->at(ii).second] = false;
  }

  int num_unc_ids =  get_unc.size() - ctrs_pos->size()*2;

  vector<string> unc_id_ranges( num_unc_ids );
  vector<string> unc_idxs( num_unc_ids );
  vector<int> unc_pos( num_unc_ids );
  map<int,int> unc_rel_pos;

  int jj = 0;
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      unc_id_ranges[jj] = block_info->orig_block()->at(ii);
      unc_idxs[jj]      = block_info->orig_idxs()->at(ii);
      unc_pos[jj]       = ii;
      unc_rel_pos.emplace(ii, jj);
      jj++;
    }
  } 
 
  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp::TensOp<double>;
template class MultiTensOp::MultiTensOp<double>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
