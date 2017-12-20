#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/tensop.h>

using namespace std;
using namespace WickUtils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_General::TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                std::shared_ptr<std::map< const std::vector<std::string>, 
                                std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>> all_ranges):  
                                idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges), orig_factor_(factor), num_idxs_(idxs.size()),
                                unique_range_blocks_(*unique_range_blocks), all_ranges_(*all_ranges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_General::TensOp_General "<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  unique_range_blocks_ptr_ = make_shared<const vector< shared_ptr< const vector<string>>>>(unique_range_blocks_);

  all_ranges_ptr_ =
  make_shared<const map< const vector<string>, tuple<bool, shared_ptr<const vector<string>>, shared_ptr< const vector<string>>, pair<int,int>>>>(all_ranges_);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MultiTensOp_General::MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                          std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                                          std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                          std::shared_ptr<std::map< const std::vector<std::string>, 
                                          std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>> all_ranges, 
                                          std::shared_ptr<std::map< const std::vector<std::string>,
                                          std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >>> split_ranges_map ) :
                         TensOp_General::TensOp_General( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, unique_range_blocks, all_ranges ),
                         cmlsizevec_(cmlsizevec),  split_ranges_map_(*split_ranges_map){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp_General::MultiTensOp_General "<< endl;
           
  cmlsizevec_ptr_ =  make_shared<const vector<int>>(cmlsizevec_);

  split_ranges_map_ptr_ = split_ranges_map;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                                 vector<bool>& aops, DataType orig_factor,
                                                 vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>), int, int > >& symmfuncs, 
                                                 vector<bool(*)(shared_ptr<vector<string>>) >& constraints,
                                                 string& Tsymm ) :
                                                 name_(name), spinfree_(true), symmfuncs_(symmfuncs), constraints_(constraints), Tsymm_(Tsymm)  {
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


  tuple< std::shared_ptr<std::map< const std::vector<std::string>, std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>,
         std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>> > test_var = generate_ranges(idxs.size(), idxs, idx_ranges);

  pair<double, double> orig_factor_tmp =  make_pair(1.0, 1.0);

  Op_dense_ = make_shared<const TensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, orig_factor_tmp, get<1>(test_var), get<0>(test_var) );

  cout << "TensOp::TensOp" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
tuple< std::shared_ptr<std::map< const std::vector<std::string>, std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>,
       std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>> >
TensOp::TensOp<DataType>::generate_ranges(int num_idxs, vector<string>& idxs,  vector<vector<string>>& idx_ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::generate_ranges" <<   endl;
  //set up loop utils

  shared_ptr<vector<string>> idxs_ptr = make_shared<vector<string>>(idxs);
  std::vector< std::shared_ptr< const std::vector<std::string>> > unique_range_blocks;
  std::map< const std::vector<std::string>, 
            std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > > all_ranges;

  
  auto apply_symmetry = [ &all_ranges, &idxs_ptr ]( const vector<string>& ranges_1,  const vector<string>& ranges_2 ) {
     const bool a_unique_range_block = true; //TODO put symmetry back in
     if (a_unique_range_block) {
       const pair<int,int> fac(1,1);
       shared_ptr<const vector<string>> ranges_2_ptr = make_shared<const vector<string>>(ranges_2);
       all_ranges.emplace(ranges_2, tie( a_unique_range_block,  ranges_2_ptr, idxs_ptr, fac ) );
     }
     return false;
  };

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( num_idxs, 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( num_idxs );
  int  num_cycles = 1;
  bool tbool = true;
  const pair<int,int> fac(1,1);

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
  shared_ptr<const vector<string>> tmp_range = make_shared<const vector<string>>(possible_ranges[0]);
  all_ranges.emplace( init_range, tie( tbool, tmp_range, idxs_ptr, fac ) );
 
  unique_range_blocks = vector< shared_ptr< const vector<string>>>(1, tmp_range );
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

  auto all_ranges_ptr = make_shared<std::map< const std::vector<std::string>, std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > >> (all_ranges);
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
  cout << "TensOp get_ctrs_tens_ranges" << endl;
 
    //puts uncontracted ranges into map 
 shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector< pair<int,int>>>(0);
 for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
   shared_ptr<vector<pair<int,int>>>  ReIm_factors = make_shared< vector<pair<int,int>>>(1, get<3>(rng_it->second)); 
   shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
   shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );
   shared_ptr<CtrTensorPart<DataType>>  CTP       = make_shared< CtrTensorPart<DataType> >( full_idxs, full_ranges, noctrs, ReIm_factors ); 
   cout << "CTP->myname() = " << CTP->myname() << " is going into CTP_map" <<  endl;
   CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs
 }
 
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
          shared_ptr<vector<pair<int,int>>> ReIm_factors = make_shared<vector<pair<int,int>>>(1, get<3>(rng_it->second)); 
          shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
          shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );
          shared_ptr<CtrTensorPart<DataType>> CTP = make_shared<CtrTensorPart<DataType>>( full_idxs, full_ranges, ctr_vec, ReIm_factors ); 
          cout << "CTP->myname() = " << CTP->myname() << "is going into CTP_map " << endl;
          CTP_map->emplace(CTP->myname(), CTP); //maybe should be added in with ctr_idxs.
        }
      }
    }
  }
  cout << "TensOp_Prep::get_ctrs_tens_ranges out" << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
MultiTensOp::MultiTensOp<DataType>::MultiTensOp( std::string name, bool spinfree,
                                                 std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>& orig_tensors ):
                                                 TensOp::TensOp<DataType>( name, spinfree ),
                                                 orig_tensors_(orig_tensors), num_tensors_(orig_tensors.size()) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int num_idxs = 0;
    
  vector<string> idxs;
  vector<vector<string>> idx_ranges;
  vector<bool> aops;
 
  vector<int> cmlsizevec(num_tensors_);
  vector<int> plus_ops;
  vector<int> kill_ops;
 
  int xx = 0;
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
  auto bob = generate_ranges(num_idxs, cmlsizevec);

  auto  unique_range_blocks = make_shared<std::vector< std::shared_ptr<const std::vector<std::string>>>>(0) ;
//  for ( auto elem : *(get<2>bob) ) 
//     unique_range_blocks->push_back( make_shared<std::vector<std::string>>(elem.first) ); 

  Op_dense_ = make_shared<const MultiTensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, make_pair(1.0,1.0), cmlsizevec,  unique_range_blocks, get<1>(bob), get<2>(bob) ); 
  CTP_map  = make_shared< map< string, shared_ptr<CtrTensorPart<DataType>> >>();
  CMTP_map = make_shared< map< string, shared_ptr<CtrMultiTensorPart<DataType>> >>();
   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
tuple< shared_ptr<map< const vector<string>, tuple<bool,                            shared_ptr<const vector<string>>,                     shared_ptr<const vector<string>>, shared_ptr<vector<pair<int,int>>> > >>,
       shared_ptr<map< const vector<string>, tuple<bool,                            shared_ptr<const vector<string>>,                     shared_ptr<const vector<string>>, pair<int,int>>  >>, 
       shared_ptr<map< const vector<string>, tuple< shared_ptr<const vector<bool>>, shared_ptr<vector<shared_ptr<const vector<string>>>>, shared_ptr<vector<pair<int,int>>> >>>>
MultiTensOp::MultiTensOp<DataType>::generate_ranges(int num_idxs_, vector<int>& cmlsizevec   ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;

  vector< map< const vector<string>, tuple<bool, shared_ptr<const vector<string>>, shared_ptr<const vector<string>>, pair<int,int>> >::const_iterator> rng_maps;

  vector<int> posvec( num_tensors_, 0 ); 

  auto combined_ranges  = make_shared<map< const vector<string>, tuple< bool,                           shared_ptr<const vector<string>>,                     shared_ptr<const vector<string>>, shared_ptr<vector<pair<int,int>>> >>>();
  auto all_ranges       = make_shared<map< const vector<string>, tuple< bool,                           shared_ptr<const vector<string>>,                     shared_ptr<const vector<string>>, pair<int,int>> >>();
  auto split_ranges_map = make_shared<map< const vector<string>, tuple< shared_ptr<const vector<bool>>, shared_ptr<vector<shared_ptr<const vector<string>>>>, shared_ptr<vector<pair<int,int>>> >>>();

  if ( num_tensors_ > 1 ) { 
 
    for (shared_ptr<TensOp::TensOp<DataType>> Ten : orig_tensors_ )
      rng_maps.push_back(Ten->all_ranges()->begin());
    
    bool ham = true; // Generate all possible combinations of different ranges on different tensors in multitensop
                     // TODO recheck this; I remember it took a while, and there was a reason for doing this without forvecs, but 
                     //      cannot for the life of me work out or remember why... The reason involved duplication, I think. 
    do{
      for (auto ii = rng_maps.size()-1; ii != 0; ii--){
        
        if ( (ii == posvec.size()-1) && (posvec[ii] == orig_tensors_[ii]->all_ranges()->size()-1 && (posvec[ii-1] == orig_tensors_[ii-1]->all_ranges()->size()-1  )  ) ) {
          ham =false;  
          break;
        } else if ( posvec[ii-1] != orig_tensors_[ii-1]->all_ranges()->size()-1 ) {
  
          rng_maps[ii-1]++;
          posvec[ii-1]++;
  
        } else if (posvec[ii-1] == orig_tensors_[ii-1]->all_ranges()->size()-1){
  
          rng_maps[ii-1] = orig_tensors_[ii-1]->all_ranges()->begin();
          posvec[ii-1] = 0;
  
          rng_maps[ii]++;
          posvec[ii]++;
  
          continue;
        }     
  
         shared_ptr<vector<shared_ptr<const vector<string>>>> unique_ranges = make_shared<vector<shared_ptr<const vector<string>>>>(num_tensors_);
         shared_ptr<vector<shared_ptr<const vector<string>>>> unique_idxs   = make_shared<vector<shared_ptr<const vector<string>>>>(num_tensors_);
         shared_ptr<pint_vec>     factors  = make_shared<pint_vec>(num_tensors_);
         shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(num_tensors_);
    
         int  Re_factor = 1;
         int  Im_factor = 1;

         vector<string>             merged_oranges (num_idxs_);
         shared_ptr<vector<string>> merged_uranges = make_shared<vector<string>>(num_idxs_);
         shared_ptr<vector<string>> merged_uqidxs  = make_shared<vector<string>>(num_idxs_);


         for (int jj = 0 ; jj != num_tensors_ ; jj++){
           isunique->at(jj)      = get<0>(rng_maps[jj]->second);
           unique_ranges->at(jj) = get<1>(rng_maps[jj]->second);  
           unique_idxs->at(jj)   = get<2>(rng_maps[jj]->second);  
        
           copy( rng_maps[jj]->first.begin(), rng_maps[jj]->first.end(), merged_oranges.begin() + cmlsizevec[jj] );             
           copy( unique_ranges->at(jj)->begin(), unique_ranges->at(jj)->end(), merged_uranges->begin() + cmlsizevec[jj] ); 
           copy( unique_idxs->at(jj)->begin(), unique_idxs->at(jj)->end(), merged_uqidxs->begin() + cmlsizevec[jj] );     

           factors->at(jj)       = get<3>(rng_maps[jj]->second);
           int Re_factor_buff = Re_factor;  
           int Im_factor_buff = Im_factor; //TODO only need one buff, but I'm very tired so could be wrong. fix later.
           Re_factor = ( Re_factor_buff * (get<3>(rng_maps[jj]->second)).first ) - ( Im_factor_buff * (get<3>(rng_maps[jj]->second)).second );
           Im_factor = ( Im_factor_buff * (get<3>(rng_maps[jj]->second)).first ) + ( Re_factor_buff * (get<3>(rng_maps[jj]->second)).second );

         }
  
         pair<int,int> combined_factor = make_pair(Re_factor, Im_factor);
         bool merged_unique = ( merged_oranges == *merged_uranges ) ?  false : true;

         combined_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, factors));
         all_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, combined_factor ));
         split_ranges_map->emplace(merged_oranges, tie(isunique, unique_ranges, factors));
  
       }
     } while(ham);
  
   } else { 
     cout << "ham " << endl;
      
     for ( auto map_it = orig_tensors_[0]->all_ranges()->begin() ; map_it != orig_tensors_[0]->all_ranges()->end(); map_it++ ){
  
       const vector<string> orig_ranges = { map_it->first } ;
       shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(vector<bool>{ get<0>(map_it->second) });
       shared_ptr<vector<shared_ptr<const vector<string>>>> unique_ranges =
       make_shared<vector<shared_ptr<const vector<string>>>>( vector<shared_ptr<const vector<string>>> { get<1>(map_it->second) });
  
       shared_ptr<pint_vec> factors = make_shared<pint_vec>( pint_vec { get<3>(map_it->second) } );
  
       split_ranges_map->emplace(orig_ranges, tie( isunique, unique_ranges, factors ) );
  
       combined_ranges->emplace( map_it->first, tie( get<0>(map_it->second), get<1>(map_it->second), get<2>(map_it->second), factors ) );
     }
   }
   
  cout << "leaving MultiTensOp::generate_ranges()" << endl;
  return tie(combined_ranges, all_ranges, split_ranges_map) ;
  
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
      if (rng_it->first[xx] != "act") {
        check=false;
        break;
      }
    }
    if(!check){
      continue;
    } else {
      enter_into_CMTP_map(*noctrs, get<3>(rng_it->second), rng_it->first );
    }
  }

  //puts_contractions, with specified ranges, into the map
  for ( int nctrs = 1 ; nctrs != (Op_dense_->num_idxs()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> ctr_lists = get_unique_pairs( Op_dense_->plus_ops(), Op_dense_->kill_ops(), nctrs );
//    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {  print_pair_vector( *ctr_vec , "ctr_vec orig" ); cout << endl; }

    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = Op_dense_->all_ranges()->begin(); rng_it != Op_dense_->all_ranges()->end(); rng_it++) {
  
        bool valid =true;
        //checks ranges for contractions match
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
            valid = false;
            break;
          }
        }
  
        //checks all uncontracted indexes are active. Should call constraint functions instead
        if (!valid) {
          continue;
        } else {
          vector<bool> unc_get(rng_it->first.size(),true);
          for( pair<int,int> ctr_pos : *ctr_vec ){
            unc_get[ctr_pos.first]  = false;
            unc_get[ctr_pos.second] = false;
          }
          for (int xx = 0; xx !=rng_it->first.size() ; xx++ ) {
            if (unc_get[xx]  && (rng_it->first[xx].substr(0,3) != "act") ) {
              valid = false;
              break;
            }
          }
        } 


        if (valid)
          enter_into_CMTP_map(*ctr_vec, get<3>(rng_it->second), rng_it->first );
        
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
  // TODO tidy up this, it seems the definition of cmlsivevec has changed.
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
    
     if( sameT_ctrs_pos[ii].size() != 0 ) {
       CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, make_shared<vector<pair<int,int>>>(sameT_ctrs_pos[ii]), ReIm_factor_vec ) ; 
     } else { 
       CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, no_ctrs, ReIm_factor_vec ); 
     }
  
     CTP_map->emplace(CTP_vec->at(ii)->name, CTP_vec->at(ii)); 
  }
  
  shared_ptr<CtrMultiTensorPart<DataType>> CMTP = make_shared<CtrMultiTensorPart<DataType> >(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos)); 
  
  CMTP_map->emplace(CMTP->myname(), CMTP); 
  CTP_map->emplace(CMTP->myname(), CMTP); 
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp::TensOp<double>;
template class MultiTensOp::MultiTensOp<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
