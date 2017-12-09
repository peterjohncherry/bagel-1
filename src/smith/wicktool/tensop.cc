#include <bagel_config.h>
#ifdef COMPILE_SMITH
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/tensop.h>

 //#include "wickutils.h"
 //#include "tensop.h"

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::initialize( vector<string> orig_idxs, vector<vector<string>> orig_idx_ranges,
                                vector<bool> orig_aops, int orig_factor, string orig_Tsymm, string orig_psymm  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "TensOp::initialize" <<   endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
           
   Tsymm_     = orig_Tsymm;
   psymm_     = orig_psymm;

   idxs_       = make_shared<vector<string>>(orig_idxs);
   idx_ranges_ = make_shared<vector<vector<string>>>(orig_idx_ranges);
   aops_       = make_shared<vector<bool>>(orig_aops);
   contracted  = false;

   generate_ranges( idx_ranges_, psymm_ );
   
   plus_ops_ = make_shared<vector<int>>();
   kill_ops_ = make_shared<vector<int>>();

   for ( int ii =0 ; ii != aops_->size() ; ii++ ) {
     if(aops_->at(ii)) {
       plus_ops_->push_back(ii);
     } else {
       kill_ops_->push_back(ii);
     }
   }

  plus_combs_ = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(plus_ops_->size()+1);
  kill_combs_ = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(kill_ops_->size()+1);
 
  for (int  ii = 0 ; ii != (plus_ops_->size() +1) ; ii++) {
    plus_combs_->at(ii) =  get_N_in_M_combsX( plus_ops_, ii );
    kill_combs_->at(ii) =  get_N_in_M_combsX( kill_ops_, ii );
  }

  CTP_map = make_shared< map< string, shared_ptr<CtrTensorPart<DType>> >>();

  cout << "TensOp::initialize" << endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::generate_ranges( shared_ptr<vector<vector<string>>> idx_ranges, string symm_type){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "TensOp::generate_ranges" <<   endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::generate_ranges" <<   endl;

  //set up loop utils
  shared_ptr<vector<int>> fvec = make_shared<vector<int>>(idx_ranges->size(), 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>();
  int  num_cycles = 1;
  bool tbool = true;
  pair<int,int> fac(1,1);
  all_ranges_ = make_shared<std::map< std::vector<std::string>, 
                                      std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >>();

  for( vector<string> range : *idx_ranges){
    maxs->push_back(range.size()-1);
    num_cycles *= range.size();
  }

  //generate all possible ranges
  vector<vector<string>> possible_ranges(0);
  for (int ii = 0 ; ii != num_cycles ; ii++){
    vector<string> new_range(idx_ranges->size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range[jj] = idx_ranges->at(jj)[fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges.push_back(new_range);
    fvec_cycle(fvec, maxs);
  } 
  
  //Apply symmetry operations to remove unnecessary ranges   
  orb_ranges = make_shared< map< vector<string>, pair<int,int> > >();
  orb_ranges->emplace(possible_ranges[0], make_pair(1,1) );
  
  shared_ptr<vector<string>> tmp_range =  make_shared<vector<string>>(possible_ranges[0]);
  all_ranges_->emplace(possible_ranges[0], tie(tbool, tmp_range , idxs_, fac ) );

  for (int ii = 1 ; ii!=possible_ranges.size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
       if(apply_symmetry(possible_ranges[kk], possible_ranges[ii])){
         break;
       }
       if(kk == ii-1){
       orb_ranges->emplace(possible_ranges[ii], make_pair(1,1) );
         break;
       }
     } 
  }
  cout << "TensOp::generate_ranges" <<   endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Does nothing for now, abandon symmetry until you have fixed the mapping problem
//I think it is better (safer) if it takes the arguments by value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool TensOp<DType>::apply_symmetry(vector<string> ranges_1, vector<string> ranges_2  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "TensOp::apply_symmetry" <<   endl;
   

  bool tbool = true; 
  pair<int,int> fac(1,1);
  shared_ptr<vector<string>> ranges_2_ptr = make_shared<vector<string>>(ranges_2);
  all_ranges_->emplace(ranges_2, tie(tbool, ranges_2_ptr, idxs_, fac ) );
  return false;
}; 


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool TensOp<DType>::apply_symmetry(shared_ptr<vector<string>> ranges_1, shared_ptr<vector<string>> ranges_2  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "TensOp::apply_symmetry" <<   endl;

//  for (auto symmop : symmfuncs_) {
//    if( *ranges_1 == *(get<0>(symmop)(ranges_2)) ){
//     bool fbool = false; 
//      pair<int,int> sfac(get<1>(symmop), get<2>(symmop));
//      shared_ptr<vector<string>> transformed_idxs = make_shared<vector<string>>(*(get<0>(symmop)(idxs)));
//      all_ranges_->emplace(ranges_2, tie( fbool, ranges_1, transformed_idxs, sfac));
//      get<2>(all_ranges_->at(ranges_2)) = transformed_idxs; //WHY DO I NEED TO DO THIS TO MAKE IT WORK? JUST EMPLACING IS NOT ENOUGH!!! VERY SUSPICIOUS! Well, I know why now...
//      return true;
//    }
//  }
//
// bool tbool = true; 
//  pair<int,int> fac(1,1);
//  all_ranges_->emplace(ranges_2, tie(tbool, ranges_2, idxs, fac ) );
  return false;
}; 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and ranges
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
  cout << "TensOp get_ctrs_tens_ranges" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp get_ctrs_tens_ranges" << endl;

  //puts uncontracted ranges into map 
  shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector<pair<int,int>>>(0);
  for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
    shared_ptr<vector<pair<int,int>>> ReIm_factors = make_shared<vector<pair<int,int>>>(1, get<3>(rng_it->second)); 
    shared_ptr<vector<string>>        full_ranges  = make_shared<vector<string>>(rng_it->first);
    shared_ptr<CtrTensorPart<DType>>  CTP          = make_shared< CtrTensorPart<DType> >(idxs_, full_ranges, noctrs, ReIm_factors ); 
    cout << "CTP->myname() = " << CTP->myname() << " is going into CTP_map" <<  endl;
    CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs.
  }

  //puts_contracted ranges into map
  for ( int nctrs = 1 ; nctrs != (idxs_->size()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>  ctr_lists = get_unique_pairs(plus_ops_, kill_ops_, nctrs);
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
          shared_ptr<CtrTensorPart<DType>> CTP = make_shared<CtrTensorPart<DType>>(idxs_, full_ranges, ctr_vec, ReIm_factors ); 
          cout << "CTP->myname() = " << CTP->myname() << "is going into CTP_map " << endl;
          CTP_map->emplace(CTP->myname(), CTP); //maybe should be added in with ctr_idxs.
        }
      }
    }
  }
  cout << "TensOp get_ctrs_tens_ranges out" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool TensOp<DType>::satisfies_constraints( shared_ptr<vector<string>> ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for (auto OK : constraints_) 
       if( !OK(ranges))
         return false;
    return true;
}; 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool TensOp<DType>::satisfies_constraints( vector<string> ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     shared_ptr<vector<string>> ranges_ptr = make_shared<vector<string>>(ranges);
     for (auto OK : constraints_) 
       if( !OK(ranges_ptr))
         return false;
    return true;
}; 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Performs hermitian conjugate on tensor indexes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::hconj() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  reverse(idxs_->begin(), idxs_->end());

  for( vector<string> range :  *idx_ranges_)
    reverse(range.begin(), range.end());

  for (bool op : *aops_)
    op =!op;      

//  for (auto range_it = orb_ranges->begin(); range_it!=orb_ranges->end(); range_it++ ){
//    reverse(range_it->first.begin(), range_it->first.end());
//    range_it->second.second *= -1;
//  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void MultiTensOp<DType>::initialize(std::vector<std::shared_ptr<TensOp<DType>>> orig_tensors ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::initialize" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp::initialize" << endl;

//Loosely speaking, Multitensop.has two kinds of members:
//1) Members inherited from the derived class /obtained by merging the members (e.g. the indexes and position
//of creation and annihilation operators) from the TensOps from which the MultiTensOp is constructed. 
//2) "Split" members; vectors where each element corresponds to a constituent TensOp of the MultiTensor.
//   The split members are needed to keep  track of signs resulting from symmetry operations.

  orig_tensors_ = orig_tensors;
  num_tensors = orig_tensors.size();

  int total_size = 0;
  for ( shared_ptr<TensOp<DType>> tens : orig_tensors_ ) 
     total_size +=  tens->idxs()->size();

  split_idxs = make_shared<vector<shared_ptr<vector<string>>>>(0) ;
 
  plus_ops = make_shared<vector<int>>(0);
  kill_ops = make_shared<vector<int>>(0);
  idxs     = make_shared<vector<string>>(0);
  aops     = make_shared<vector<bool>>(0);


  cmlsizevec = make_shared<vector<int>>(num_tensors);

  idx_ranges = make_shared<vector<vector<string>>>(0);
  Tsizes     = make_shared<vector<int>>(num_tensors); 
  names = make_shared<vector<string>>(num_tensors);

  int xx= 0;
  for (int ii = 0 ; ii != num_tensors; ii++) {
     split_idxs->push_back(orig_tensors_[ii]->idxs());

    for ( string elem : *(orig_tensors_[ii]->idxs()))
      idxs->push_back(elem);
    for ( int elem : *(orig_tensors_[ii]->plus_ops()))
      plus_ops->push_back(elem+xx);
    for ( int elem : *(orig_tensors_[ii]->kill_ops()))
      kill_ops->push_back(elem+xx);
    for ( vector<string> elem : *(orig_tensors_[ii]->idx_ranges()))
      idx_ranges->push_back(elem);

    for ( bool elem : *(orig_tensors_[ii]->aops()))
      aops->push_back(elem);

    names->at(ii) = orig_tensors_[ii]->name();
    Tsizes->at(ii) = orig_tensors_[ii]->idxs()->size();
    cmlsizevec->at(ii) = xx;
    xx+=orig_tensors_[ii]->idxs()->size();
  } 

  generate_ranges();
 
  plus_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(plus_ops->size()+1);
  kill_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(kill_ops->size()+1);
 
  for (int  ii = 0 ; ii != (plus_ops->size() +1) ; ii++) {
    plus_combs->at(ii) =  get_N_in_M_combsX( plus_ops, ii );
    kill_combs->at(ii) =  get_N_in_M_combsX( kill_ops, ii );
  }

  CTP_map = make_shared< map< string, shared_ptr<CtrTensorPart<DType>> >>();
  CMTP_map = make_shared< map< string, shared_ptr<CtrMultiTensorPart<DType>> >>();
  cout << "Leaving MultiTensOp::initialize" << endl;
  return;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void MultiTensOp<DType>::generate_ranges(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::generate_ranges()" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;

  vector<map< vector<string>,  tuple<bool, shared_ptr<vector<string>>, shared_ptr<vector<string>>,  pair<int,int> > >::iterator> rng_maps;
  vector<int> posvec( num_tensors, 0 ); 

  combined_ranges = make_shared< map< vector<string>,
                                      tuple<bool, shared_ptr<vector<string>>,  shared_ptr<vector<string>>, shared_ptr<vector<pair<int,int>>> > >>();

  split_ranges    = make_shared< map< vector<vector<string>>,
                                      tuple< shared_ptr<vector<bool>>, shared_ptr<vector<shared_ptr<vector<string>>>>, shared_ptr<vector<pair<int,int>>> > >>();

  if ( num_tensors > 1 ) { 

    for (shared_ptr<TensOp<DType>> Ten : orig_tensors_ )
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

        vector<vector<string>> orig_ranges(num_tensors);
        shared_ptr<vector<shared_ptr<vector<string>>>> unique_ranges = make_shared<vector<shared_ptr<vector<string>>>>(num_tensors);
        shared_ptr<vector<shared_ptr<vector<string>>>> unique_idxs   = make_shared<vector<shared_ptr<vector<string>>>>(num_tensors);
        shared_ptr<pint_vec>     factors  = make_shared<pint_vec>(num_tensors);
        shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(num_tensors);
    
        for (int jj = 0 ; jj != num_tensors ; jj++){
          orig_ranges[jj] = rng_maps[jj]->first;
          num_orig_ranges += orig_ranges[jj].size(); 

          isunique->at(jj) = get<0>(rng_maps[jj]->second);

          unique_ranges->at(jj) = get<1>(rng_maps[jj]->second);  
          num_unique_ranges += unique_ranges->at(jj)->size(); 

          unique_idxs->at(jj) = get<2>(rng_maps[jj]->second);  
          num_unique_idxs += unique_idxs->at(jj)->size(); 

          factors->at(jj) = get<3>(rng_maps[jj]->second) ;
        }
    
        split_ranges->emplace(orig_ranges, tie(isunique, unique_ranges, factors));
        
        vector<string>             merged_oranges(num_orig_ranges);
        shared_ptr<vector<string>> merged_uranges = make_shared<vector<string>>(num_unique_ranges);
        shared_ptr<vector<string>> merged_uqidxs = make_shared<vector<string>>(num_unique_idxs);

        vector<string>::iterator merged_oranges_it = merged_oranges.begin();
        vector<string>::iterator merged_uqidxs_it  = merged_uqidxs->begin();
        vector<string>::iterator merged_uranges_it = merged_uranges->begin();

        for( int qq = 0 ; qq != num_tensors; qq++ ){ 

          copy( orig_ranges[qq].begin(), orig_ranges[qq].end(), merged_oranges_it );             
          copy( unique_ranges->at(qq)->begin(), unique_ranges->at(qq)->end(), merged_uranges_it ); 
          copy( unique_idxs->at(qq)->begin(), unique_idxs->at(qq)->end(), merged_uqidxs_it );     

          merged_oranges_it += orig_ranges[qq].size();
          merged_uranges_it += unique_ranges->at(qq)->size();
          merged_uqidxs_it += unique_idxs->at(qq)->size();

        }
        
        bool merged_unique = ( merged_oranges == *merged_uranges ) ?  false : true;
        combined_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, factors));

      }
    } while(ham);

  } else { 
     
    for ( auto map_it = orig_tensors_[0]->all_ranges()->begin() ; map_it != orig_tensors_[0]->all_ranges()->end(); map_it++ ){

      vector<vector<string>> orig_ranges = { map_it->first } ;
      shared_ptr<vector<bool>> isunique = make_shared<vector<bool>>(vector<bool>{ get<0>(map_it->second) });
      shared_ptr<vector<shared_ptr<vector<string>>>> unique_ranges =
      make_shared<vector<shared_ptr<vector<string>>>>( vector<shared_ptr<vector<string>>> { get<1>(map_it->second) });

      shared_ptr<pint_vec> factors = make_shared<pint_vec>( pint_vec { get<3>(map_it->second) } );

      split_ranges->emplace(orig_ranges, tie( isunique, unique_ranges, factors ) );

      combined_ranges->emplace( map_it->first, tie( get<0>(map_it->second), get<1>(map_it->second), get<2>(map_it->second), factors ) );
    }
  }    
  cout << "leaving MultiTensOp::generate_ranges()" << endl;

  return;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and range, looks totally inefficient, but generating everything and then checking in sequence
//seems to work out faster; it's either check everything lots of times, or regenerate it lots of times. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void MultiTensOp<DType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp get_ctrs_tens_ranges" << endl;
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp::get_ctrs_tens_ranges " <<  endl;
  //puts uncontracted ranges into map 
  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);
  
  //silly, should just test {act,act,....} against contraints, and check if act in each ranges... 
  for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {
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
      enter_into_CMTP_map(*noctrs, get<3>(rng_it->second), make_shared<vector<string>>(rng_it->first) );
    }
  }

  //puts_contractions, with specified ranges, into the map
  for ( int nctrs = 1 ; nctrs != (idxs->size()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> ctr_lists = get_unique_pairs(plus_ops, kill_ops, nctrs);

    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {

        bool valid =true;
        //checks ranges for contractions match
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
            valid = false;
            break;
          }
        }

        if (!valid) 
          continue;

         //checks all uncontracted indexes are active. Should call constraint functions instead
        if (valid){	 
          vector<bool> unc_get(rng_it->first.size(),true);
          for (pair<int,int> ctr_pos : *ctr_vec){
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

        if (valid) { 
          enter_into_CMTP_map(*ctr_vec, get<3>(rng_it->second), make_shared<vector<string>>(rng_it->first) );
        }
      }
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void MultiTensOp<DType>::enter_into_CMTP_map(pint_vec ctr_pos_list, shared_ptr<vector<pair<int,int>>> ReIm_factors, shared_ptr<vector<string>> id_ranges ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "MultiTensOp::enter_into_CMTP_map" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  shared_ptr<pstr_vec>ctr_idxs_in   = make_shared<pstr_vec>(0);
  shared_ptr<pstr_vec>ctr_ranges_in = make_shared<pstr_vec>(0);
  shared_ptr<pstr_vec>ctr_spins_in  = make_shared<pstr_vec>(0);

  shared_ptr<vector<shared_ptr<vector<string>>>> unc_idxs_in   = make_shared<vector<shared_ptr<vector<string>>>>(0);
  shared_ptr<vector<shared_ptr<vector<string>>>> unc_ranges_in = make_shared<vector<shared_ptr<vector<string>>>>(0);
  shared_ptr<vector<shared_ptr<vector<bool>>>>   unc_aops_in   = make_shared<vector<shared_ptr<vector<bool>>>>(0);

  shared_ptr<vector<shared_ptr<CtrTensorPart<DType>>>> CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart<DType>> >> (num_tensors); 
  vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
  vector<vector<pair<int,int>>> sameT_ctrs_pos( num_tensors,  pint_vec(0));

 
  //seperate contractions into those on the same tensor, and those between different tensors 
  for ( pair<int,int> ctr_pos : ctr_pos_list ) {

    pair<int,int> ctr1;
    pair<int,int> ctr2;

    for ( int ii =cmlsizevec->size()-1 ; ii !=-1 ; ii-- ) { 
      if(ctr_pos.first >= cmlsizevec->at(ii) ){
        ctr1 = make_pair(ii,  ctr_pos.first- cmlsizevec->at(ii));
        break;
      }
    }

    for ( int ii =cmlsizevec->size()-1 ; ii !=-1 ; ii-- ) { 
      if(ctr_pos.second >= cmlsizevec->at(ii) ){
        ctr2 = make_pair(ii,  ctr_pos.second- cmlsizevec->at(ii));
        break;
      }
    }
    if (ctr1.first == ctr2.first) {
      sameT_ctrs_pos.at(ctr1.first).push_back(make_pair(ctr1.second, ctr2.second));  
    } else {
      diffT_ctrs_pos.push_back(make_pair(ctr1,ctr2));
    }

  }  

  //get_ranges for individual tensors
  for (int ii = 0 ; ii !=num_tensors; ii++ ){ 
     shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(0);
     for(int jj = cmlsizevec->at(ii) ; jj != cmlsizevec->at(ii)+orig_tensors_[ii]->idxs()->size(); jj++)
       TS_id_ranges->push_back(id_ranges->at(jj));
     
   
     CTP_vec->at(ii) = make_shared< CtrTensorPart<DType> >( orig_tensors_[ii]->idxs(), TS_id_ranges,
                                                            make_shared<vector<pair<int,int>>>(sameT_ctrs_pos.at(ii)),
                                                            make_shared<vector<pair<int,int>>>(1, ReIm_factors->at(ii)) ) ; 

     CTP_map->emplace(CTP_vec->at(ii)->name, CTP_vec->at(ii)); 
  }
  
  shared_ptr<CtrMultiTensorPart<DType>> CMTP = make_shared<CtrMultiTensorPart<DType> >(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos)); 

  CMTP_map->emplace(CMTP->myname(), CMTP); 
  CTP_map->emplace(CMTP->myname(), CMTP); 

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp<double>;
template class MultiTensOp<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
