#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/tensop.h>

using namespace std;
using namespace WickUtils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_General::TensOp_General( std::vector<int> plus_ops, std::vector<int> kill_ops,
                                std::vector< std::shared_ptr<const std::vector<std::string>>> unique_range_blocks,
                                std::map< const std::vector<std::string>, 
                                          std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                          std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>> all_ranges ) :
                    plus_ops_(plus_ops), kill_ops_(kill_ops), unique_range_blocks_(unique_range_blocks), all_ranges_(all_ranges){                         
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_General::TensOp_General "<< endl;

  plus_ops_ptr_ = std::make_shared<const std::vector<int>>(plus_ops);
  kill_ops_ptr_ = std::make_shared<const std::vector<int>>(kill_ops);

  unique_range_blocks_ptr_ = std::make_shared<const std::vector< std::shared_ptr< const std::vector<std::string>>>>(unique_range_blocks);

  all_ranges_ptr_ = std::make_shared<const std::map< const std::vector<std::string>,
                                                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                     std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>(all_ranges);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                                 vector<bool>& aops, DataType orig_factor,
                                                 vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>), int, int > >& symmfuncs, 
                                                 vector<bool(*)(shared_ptr<vector<string>>) >& constraints,
                                                 string& Tsymm ) :
                                                 name_(name), idxs_(idxs), idx_ranges_(idx_ranges), aops_(aops), 
                                                 orig_factor_(orig_factor), symmfuncs_(symmfuncs), constraints_(constraints),
                                                 Tsymm_(Tsymm), num_idxs_(idxs_.size()) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::TensOp" <<   endl;
           
  idxs_const_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_const_ptr_ =  make_shared<const vector<bool>>(aops_);

  generate_ranges();

  cout << "TensOp::TensOp" << endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_ranges(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::generate_ranges" <<   endl;
  //set up loop utils
  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( num_idxs_, 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( num_idxs_ );
  int  num_cycles = 1;
  bool tbool = true;
  const pair<int,int> fac(1,1);

  for( int ii = 0; ii != idx_ranges_.size(); ii++ ){
    maxs->at(ii) = ( idx_ranges_[ii].size()-1 );
    num_cycles *= idx_ranges_[ii].size();
  }

  //generate all possible ranges
  vector<vector<string>> possible_ranges(0);
  for (int ii = 0 ; ii != num_cycles ; ii++){
    vector<string> new_range(idx_ranges_.size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range[jj] = idx_ranges_[jj][fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges.push_back(new_range);
    fvec_cycle( fvec, maxs );
  } 
 
  std::vector< std::shared_ptr< const std::vector<std::string>> > unique_range_blocks;

  std::map< const std::vector<std::string>, 
            std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > > all_ranges;

  // Initialize range maps and add in first range to start things off.
  const vector<string> init_range = possible_ranges[0];
  shared_ptr<const vector<string>> tmp_range = make_shared<const vector<string>>(possible_ranges[0]);
  all_ranges.emplace( init_range, tie(tbool, tmp_range , idxs_const_ptr_, fac ) );
 
  unique_range_blocks = vector< shared_ptr< const vector<string>>>(1, tmp_range );
  //Apply symmetry operations to remove unnecessary ranges   
  for (int ii = 1 ; ii!=possible_ranges.size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
      if(apply_symmetry(possible_ranges[kk], possible_ranges[ii],  all_ranges)){ //note that this sets the elements of all_ranges_
        break;
      }
      if(kk == ii-1){
        unique_range_blocks.push_back(make_shared<const vector<string>>(possible_ranges[ii]));
        break;
      }
    } 
  }

  std::vector<int> plus_ops;
  std::vector<int> kill_ops;
 
  for ( int ii =0 ; ii != aops_.size() ; ii++ ) {
      if(aops_[ii]) {
        plus_ops.push_back(ii);
      } else {
        kill_ops.push_back(ii);
      }
  }

  tensop_dense_ = make_shared<const TensOp_General>( plus_ops, kill_ops, unique_range_blocks, all_ranges );

  cout << "TensOp::generate_ranges()" <<   endl;
  return;
}
//TODO  make it so this actually applies the symmetry; does nothing for now whilst testing
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
bool TensOp::TensOp<DataType>::apply_symmetry( const vector<string>& ranges_1,  const vector<string>& ranges_2,
std::map< const std::vector<std::string>, 
            std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > >& all_ranges){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::apply_symmetry" <<   endl;

  const bool a_unique_range_block = true; 
  if (a_unique_range_block) {
    const pair<int,int> fac(1,1);
    shared_ptr<const vector<string>> ranges_2_ptr = make_shared<const vector<string>>(ranges_2);
    all_ranges.emplace(ranges_2, tie( a_unique_range_block,  ranges_2_ptr, idxs_const_ptr_, fac ) );
  }
  return false;
}; 

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
   shared_ptr<const vector<string>>   full_ranges  = make_shared<const vector<string>>(rng_it->first);
   //shared_ptr<CtrTensorPart<DataType>>  CTP       = make_shared< CtrTensorPart<DataType> >(idxs_, full_ranges, noctrs, ReIm_factors ); 
   cout << "CTP->myname() = "; // << CTP->myname() << " is going into CTP_map" <<  endl;
//   CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs.
 }
 
    //puts_contracted ranges into map
  for ( int nctrs = 1 ; nctrs != (idxs_.size()/2)+1 ; nctrs++ ){
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
          shared_ptr<const vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
        //  shared_ptr<CtrTensorPart<DataType>> CTP = make_shared<CtrTensorPart<DataType>>(idxs_, full_ranges, ctr_vec, ReIm_factors ); 
          cout << "CTP->myname() = ";// << CTP->myname() << "is going into CTP_map " << endl;
        //  CTP_map->emplace(CTP->myname(), CTP); //maybe should be added in with ctr_idxs.
        }
      }
    }
  }
  cout << "TensOp_Prep::get_ctrs_tens_ranges out" << endl;
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
MultiTensOp_Prep::MultiTensOp_Prep<DataType>::MultiTensOp_Prep( std::string name , bool spinfree,
                                                                std::vector<std::shared_ptr<TensOp_Interface::TensOp_Interface<DataType>>>& orig_tensors ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  num_idxs_ = 0;
  orig_tensors_ = orig_tensors;
 
  for ( shared_ptr<TensOp_Interface::TensOp_Interface<DataType>> tens : orig_tensors ) 
    num_idxs_ +=  tens->num_idxs();
  
  aops_       = vector<bool>(num_idxs_);
  idxs_       = vector<string>(num_idxs_);
  idx_ranges_ = vector<vector<string>>(num_idxs_);
  
  plus_ops_ = vector<int>(0);
  kill_ops_ = vector<int>(0);
 
  int xx = 0;
  for ( shared_ptr<TensOp_Interface::TensOp_Interface<DataType>> orig_tens : orig_tensors) {
 
    copy( orig_tens->idxs()->begin(), orig_tens->idxs()->end(), idxs_.begin()+xx );
    copy( orig_tens->aops()->begin(), orig_tens->aops()->end(), aops_.begin()+xx );
    copy( orig_tens->idxs()->begin(), orig_tens->idxs()->end(), idxs_.begin()+xx );
 
    // pushing back here as, in principle, we don't know how long these are (particle number may not be conserved).
    for ( int elem : *(orig_tens->plus_ops()) )
      plus_ops_.push_back(elem+xx);
 
    for ( int elem : *(orig_tens->kill_ops()) )
      kill_ops_.push_back(elem+xx);
 
    xx+=orig_tens->num_idxs();
  } 
 
  num_tensors_ = orig_tensors.size();
  
  cmlsizevec_ = vector<int>(num_tensors_);
  names_  = vector<string>(num_tensors_);
  Tsizes_ = vector<int>(num_tensors_); 
 
  for ( int ii = 0;  ii != orig_tensors.size() ; ii++ ) {
    Tsizes_[ii]  = orig_tensors[ii]->num_idxs();
    cmlsizevec_[ii]  = num_idxs_;
    num_idxs_ += Tsizes_[ii];
    names_[ii] = orig_tensors[ii]->name();
  }
  
  generate_ranges();

 // CTP_map = make_shared< map< string, shared_ptr<CtrTensorPart<DataType>> >>();
  cout << "TensOp_Prep::TensOp_Prep<DataType>::TensOp_Prep" << endl;
   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp_Prep::MultiTensOp_Prep<DataType>::generate_ranges(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp_Prep::generate_ranges()" << endl;

  vector< map< const vector<string>, tuple<bool, shared_ptr<const vector<string>>, shared_ptr<const vector<string>>, pair<int,int>> >::iterator> rng_maps;
 
  vector<int> posvec( num_tensors_, 0 ); 
 
  combined_ranges = make_shared< map< const vector<string>,
                                      tuple<bool, shared_ptr<const vector<string>>,  shared_ptr<const vector<string>>, shared_ptr<vector<pair<int,int>>> > >>();
 
  split_ranges_map = make_shared< map< const vector<string>,
                                       tuple< shared_ptr<const vector<bool>>, shared_ptr<vector<shared_ptr<const vector<string>>>>, shared_ptr<vector<pair<int,int>>> > >>();
 
  if ( num_tensors_ > 1 ) { 
 
    for (shared_ptr<TensOp_Interface::TensOp_Interface<DataType>> Ten : orig_tensors_ )
      rng_maps.push_back(Ten->all_ranges()->begin());
    
    bool ham =true; //awkward loop to generate all possible combinations of different ranges on different tensors in multitensop
    do{
      for (auto ii = rng_maps.size()-1; ii != 0; ii--){
        
        int num_unique_ranges = 0; 
        int num_unique_idxs = 0; 
        int num_orig_ranges = 0; 
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
  
         vector<vector<string>> orig_ranges(num_tensors_);
         shared_ptr<const vector<shared_ptr<const vector<string>>>> unique_ranges = make_shared<const vector<shared_ptr<const vector<string>>>>(num_tensors_);
         shared_ptr<const vector<shared_ptr<const vector<string>>>> unique_idxs   = make_shared<const vector<shared_ptr<const vector<string>>>>(num_tensors_);
         shared_ptr<pint_vec>     factors  = make_shared<pint_vec>(num_tensors_);
         shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(num_tensors_);
    
 
         for (int jj = 0 ; jj != num_tensors_ ; jj++){
           orig_ranges[jj] = vector<string>(rng_maps[jj]->first);
           num_orig_ranges += orig_ranges.size(); 
  
           isunique->at(jj) = get<0>(rng_maps[jj]->second);

           unique_ranges->at(jj) = get<1>(rng_maps[jj]->second);  
           num_unique_ranges += unique_ranges->at(jj)->size(); 

           unique_idxs->at(jj) = get<2>(rng_maps[jj]->second);  
           num_unique_idxs += unique_idxs->at(jj)->size(); 

           factors->at(jj) = get<3>(rng_maps[jj]->second) ;
         }
         
         vector<string>             merged_oranges(num_orig_ranges);
         shared_ptr<vector<string>> merged_uranges = make_shared<vector<string>>(num_unique_ranges);
         shared_ptr<vector<string>> merged_uqidxs = make_shared<vector<string>>(num_unique_idxs);
  
         vector<string>::iterator merged_oranges_it = merged_oranges.begin();
         vector<string>::iterator merged_uqidxs_it  = merged_uqidxs->begin();
         vector<string>::iterator merged_uranges_it = merged_uranges->begin();
  
         for( int qq = 0 ; qq != num_tensors_; qq++ ){ 
  
           copy( orig_ranges[qq].begin(), orig_ranges[qq].end(), merged_oranges_it );             
           copy( unique_ranges->at(qq)->begin(), unique_ranges->at(qq)->end(), merged_uranges_it ); 
           copy( unique_idxs->at(qq)->begin(), unique_idxs->at(qq)->end(), merged_uqidxs_it );     
  
           merged_oranges_it += orig_ranges[qq].size();
           merged_uranges_it += unique_ranges->at(qq)->size();
           merged_uqidxs_it += unique_idxs->at(qq)->size();
         }
        
         bool merged_unique = ( merged_oranges == *merged_uranges ) ?  false : true;
         combined_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, factors));
         split_ranges_map->emplace(merged_oranges, tie(isunique, unique_ranges, factors));
  
       }
     } while(ham);
  
   } else { 
     cout << "ham " << endl;
      
     for ( auto map_it = orig_tensors_[0]->all_ranges()->begin() ; map_it != orig_tensors_[0]->all_ranges()->end(); map_it++ ){
  
       vector<vector<string>> orig_ranges = { map_it->first } ;
       shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(vector<bool>{ get<0>(map_it->second) });
       shared_ptr<vector<shared_ptr<vector<string>>>> unique_ranges =
       make_shared<vector<shared_ptr<vector<string>>>>( vector<shared_ptr<vector<string>>> { get<1>(map_it->second) });
  
       shared_ptr<pint_vec> factors = make_shared<pint_vec>( pint_vec { get<3>(map_it->second) } );
  
       split_ranges_map->emplace(orig_ranges, tie( isunique, unique_ranges, factors ) );
  
       combined_ranges->emplace( map_it->first, tie( get<0>(map_it->second), get<1>(map_it->second), get<2>(map_it->second), factors ) );
     }
   }    
   cout << "leaving MultiTensOp_prep::generate_ranges()" << endl;
 
   return;
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::generate_ranges(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::generate_ranges()" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;

 //vector< map< const vector<string>, tuple<bool, shared_ptr<const vector<string>>, shared_ptr<const vector<string>>, pair<int,int>> >::const_iterator> rng_maps;
 //
 // vector<int> posvec( num_tensors_, 0 ); 
 //
 // combined_ranges = make_shared< map< vector<string>,
 //                                     tuple<bool, shared_ptr<vector<string>>,  shared_ptr<vector<string>>, shared_ptr<vector<pair<int,int>>> > >>();
 //
 // split_ranges = make_shared< map< vector<vector<string>>,
 //                                  tuple< shared_ptr<vector<bool>>, shared_ptr<vector<shared_ptr<vector<string>>>>, shared_ptr<vector<pair<int,int>>> > >>();
 //
 // if ( num_tensors_ > 1 ) { 
 //
 //   for (shared_ptr<TensOp::TensOp<DataType>> Ten : orig_tensors_ )
 //     rng_maps.push_back(Ten->all_ranges()->begin());
 //   
 //   bool ham =true; //awkward loop to generate all possible combinations of different ranges on different tensors in multitensop
 //   do{
 //     for (auto ii = rng_maps.size()-1; ii != 0; ii--){
 //       
 //       int num_unique_ranges = 0; 
 //        int num_unique_idxs = 0; 
  //       int num_orig_ranges = 0; 
  //       if ( (ii == posvec.size()-1) && (posvec[ii] == orig_tensors_[ii]->all_ranges()->size()-1 && (posvec[ii-1] == orig_tensors_[ii-1]->all_ranges()->size()-1  )  ) ) {
  //         ham =false;  
  //         break;
  //       } else if ( posvec[ii-1] != orig_tensors_[ii-1]->all_ranges()->size()-1 ) {
  //
  //         rng_maps[ii-1]++;
  //         posvec[ii-1]++;
  //
  //       } else if (posvec[ii-1] == orig_tensors_[ii-1]->all_ranges()->size()-1){
  //
  //         rng_maps[ii-1] = orig_tensors_[ii-1]->all_ranges()->begin();
  //         posvec[ii-1] = 0;
  //
  //         rng_maps[ii]++;
  //         posvec[ii]++;
  //
  //         continue;
  //       }     
  //
 //        vector<vector<string>> orig_ranges(num_tensors_);
  //       shared_ptr<vector<shared_ptr<vector<string>>>> unique_ranges = make_shared<vector<shared_ptr<vector<string>>>>(num_tensors_);
  //       shared_ptr<vector<shared_ptr<vector<string>>>> unique_idxs   = make_shared<vector<shared_ptr<vector<string>>>>(num_tensors_);
  //       shared_ptr<pint_vec>     factors  = make_shared<pint_vec>(num_tensors_);
  //       shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(num_tensors_);
  //   
  //       for (int jj = 0 ; jj != num_tensors_ ; jj++){
  //         orig_ranges[jj] = rng_maps[jj]->first;
  //         num_orig_ranges += orig_ranges[jj].size(); 
  //
  //         isunique->at(jj) = get<0>(rng_maps[jj]->second);
  //
  //         unique_ranges->at(jj) = get<1>(rng_maps[jj]->second);  
  //         num_unique_ranges += unique_ranges->at(jj)->size(); 
  //
  //         unique_idxs->at(jj) = get<2>(rng_maps[jj]->second);  
  //         num_unique_idxs += unique_idxs->at(jj)->size(); 
  //
  //         factors->at(jj) = get<3>(rng_maps[jj]->second) ;
  //       }
  //   
  //        split_ranges->emplace(orig_ranges, tie(isunique, unique_ranges, factors));
  //       
  //       vector<string>             merged_oranges(num_orig_ranges);
  //       shared_ptr<vector<string>> merged_uranges = make_shared<vector<string>>(num_unique_ranges);
  //       shared_ptr<vector<string>> merged_uqidxs = make_shared<vector<string>>(num_unique_idxs);
  //
  //       vector<string>::iterator merged_oranges_it = merged_oranges.begin();
  //       vector<string>::iterator merged_uqidxs_it  = merged_uqidxs->begin();
  //       vector<string>::iterator merged_uranges_it = merged_uranges->begin();
  //
  //       for( int qq = 0 ; qq != num_tensors_; qq++ ){ 
  //
  //         copy( orig_ranges[qq].begin(), orig_ranges[qq].end(), merged_oranges_it );             
  //         copy( unique_ranges->at(qq)->begin(), unique_ranges->at(qq)->end(), merged_uranges_it ); 
  //         copy( unique_idxs->at(qq)->begin(), unique_idxs->at(qq)->end(), merged_uqidxs_it );     
  //
  //         merged_oranges_it += orig_ranges[qq].size();
  //         merged_uranges_it += unique_ranges->at(qq)->size();
  //         merged_uqidxs_it += unique_idxs->at(qq)->size();
  //
  //       }
 //        
  //       bool merged_unique = ( merged_oranges == *merged_uranges ) ?  false : true;
  //       combined_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, factors));
  //
  //     }
  //   } while(ham);
  //
  // } else { 
  //    
  //   for ( auto map_it = orig_tensors_[0]->all_ranges()->begin() ; map_it != orig_tensors_[0]->all_ranges()->end(); map_it++ ){
  //
  //     vector<vector<string>> orig_ranges = { map_it->first } ;
  //     shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(vector<bool>{ get<0>(map_it->second) });
  //     shared_ptr<vector<shared_ptr<vector<string>>>> unique_ranges =
  //     make_shared<vector<shared_ptr<vector<string>>>>( vector<shared_ptr<vector<string>>> { get<1>(map_it->second) });
  //
  //     shared_ptr<pint_vec> factors = make_shared<pint_vec>( pint_vec { get<3>(map_it->second) } );
  //
  //     split_ranges->emplace(orig_ranges, tie( isunique, unique_ranges, factors ) );
  //
  //     combined_ranges->emplace( map_it->first, tie( get<0>(map_it->second), get<1>(map_it->second), get<2>(map_it->second), factors ) );
  //   }
  // }    
  // cout << "leaving MultiTensOp::generate_ranges()" << endl;
  //
  // return;
  //
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
 //   //puts uncontracted ranges into map 
  //  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);
  //  
  //  //silly, should just test {act,act,....} against contraints, and check if act in each ranges... 
  //  for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {
  //    bool check = true;
  //    for (int xx = 0; xx !=rng_it->first.size() ; xx++ ) {
  //      if (rng_it->first[xx] != "act") {
  //        check=false;
  //        break;
  //      }
  //    }
  //    if(!check){
  //      continue;
  //    } else {
  //      enter_into_CMTP_map(*noctrs, get<3>(rng_it->second), make_shared<vector<string>>(rng_it->first) );
  //    }
  //  }
  //  //puts_contractions, with specified ranges, into the map
  //  for ( int nctrs = 1 ; nctrs != (num_idxs_/2)+1 ; nctrs++ ){
  //    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> ctr_lists = get_unique_pairs(plus_ops_, kill_ops_, nctrs);
  // 
   //   for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
   //     for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {
   //
   //       bool valid =true;
   //       //checks ranges for contractions match
   //       for (int ii = 0 ; ii != ctr_vec->size(); ii++){
   //         if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
   //           valid = false;
   //           break;
   //         }
   //       }
   //
   //       if (!valid) 
   //         continue;
   //
   //        //checks all uncontracted indexes are active. Should call constraint functions instead
   //       if (valid){	 
   //         vector<bool> unc_get(rng_it->first.size(),true);
   //         for (pair<int,int> ctr_pos : *ctr_vec){
   //           unc_get[ctr_pos.first]  = false;
  //            unc_get[ctr_pos.second] = false;
   //         }
   //         for (int xx = 0; xx !=rng_it->first.size() ; xx++ ) {
   //           if (unc_get[xx]  && (rng_it->first[xx].substr(0,3) != "act") ) {
   //             valid = false;  
   //             break;
   //           }
   //         }
   //       } 
   //
   //       if (valid) { 
   //         enter_into_CMTP_map(*ctr_vec, get<3>(rng_it->second), make_shared<vector<string>>(rng_it->first) );
   //       }
   //     }
   //   }
   // }
   // return;
}  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::enter_into_CMTP_map(pint_vec ctr_pos_list, shared_ptr<vector<pair<int,int>>> ReIm_factors, shared_ptr<vector<string>> id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::enter_into_CMTP_map" << endl;
#endif 
////////////////////////////////////////////////////////////////////////////////////////
   // shared_ptr<pstr_vec> ctr_idxs_in   = make_shared<pstr_vec>(0);
   // shared_ptr<pstr_vec> ctr_ranges_in = make_shared<pstr_vec>(0);
   // shared_ptr<pstr_vec> ctr_spins_in  = make_shared<pstr_vec>(0);
   //
   // shared_ptr<vector<shared_ptr<vector<string>>>> unc_idxs_in   = make_shared<vector<shared_ptr<vector<string>>>>(0);
   // shared_ptr<vector<shared_ptr<vector<string>>>> unc_ranges_in = make_shared<vector<shared_ptr<vector<string>>>>(0);
   // shared_ptr<vector<shared_ptr<vector<bool>>>>   unc_aops_in   = make_shared<vector<shared_ptr<vector<bool>>>>(0);
   //
   // shared_ptr<vector<shared_ptr<CtrTensorPart<DataType>>>> CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart<DataType>> >> (num_tensors_); 
   // vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
   // vector<vector<pair<int,int>>> sameT_ctrs_pos( num_tensors_,  pint_vec(0));
   //
   // //seperate contractions into those on the same tensor, and those between different tensors 
   // for ( pair<int,int> ctr_pos : ctr_pos_list ) {
   //
   //   pair<int,int> ctr1;
  //    pair<int,int> ctr2;
   //
   //   for ( int ii =cmlsizevec_->size()-1 ; ii !=-1 ; ii-- ) { 
   //     if(ctr_pos.first >= cmlsizevec_->at(ii) ){
   //       ctr1 = make_pair(ii,  ctr_pos.first- cmlsizevec_->at(ii));
   //       break;
   //     }
   //   }
   //
   //   for ( int ii =cmlsizevec_->size()-1 ; ii !=-1 ; ii-- ) { 
   //     if(ctr_pos.second >= cmlsizevec_->at(ii) ){
   //       ctr2 = make_pair(ii,  ctr_pos.second- cmlsizevec_->at(ii));
   //       break;
   //     }
   //   }
   //   if (ctr1.first == ctr2.first) {
   //     sameT_ctrs_pos.at(ctr1.first).push_back(make_pair(ctr1.second, ctr2.second));  
   //   } else {
   //     diffT_ctrs_pos.push_back(make_pair(ctr1,ctr2));
   //   }
   //
   // }  
   //
   // //get_ranges for individual tensors
   // for (int ii = 0 ; ii !=num_tensors_; ii++ ){ 
   //    shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(0);
   //    for( int jj = cmlsizevec_->at(ii) ; jj != cmlsizevec_->at(ii)+orig_tensors_[ii]->idxs()->size(); jj++ )
   //     TS_id_ranges->push_back(id_ranges->at(jj));
   //   
   // 
   //   CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( orig_tensors_[ii]->idxs(), TS_id_ranges,
   //                                                          make_shared<vector<pair<int,int>>>(sameT_ctrs_pos.at(ii)),
   //                                                          make_shared<vector<pair<int,int>>>(1, ReIm_factors->at(ii)) ) ; 
   //
   //   CTP_map->emplace(CTP_vec->at(ii)->name, CTP_vec->at(ii)); 
  // }
  // 
  // shared_ptr<CtrMultiTensorPart<DataType>> CMTP = make_shared<CtrMultiTensorPart<DataType> >(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos)); 
  //
  // CMTP_map->emplace(CMTP->myname(), CMTP); 
  // CTP_map->emplace(CMTP->myname(), CMTP); 
  //
  // return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp::TensOp<double>;
template class MultiTensOp::MultiTensOp<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
