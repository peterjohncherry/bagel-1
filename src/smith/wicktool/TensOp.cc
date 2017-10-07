#include <bagel_config.h>
#ifdef COMPILE_SMITH
 #include <src/smith/wicktool/WickUtils.h>
 #include <src/smith/wicktool/TensOp.h>

 //#include "WickUtils.h"
 //#include "TensOp.h"

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::initialize(vector<string> orig_idxs, vector<vector<string>> orig_idx_ranges,
                vector<bool> orig_aops, int orig_factor, string orig_Tsymm, string orig_psymm  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
cout << "TensOp::initialize" <<   endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
           
   Tsymm_     = orig_Tsymm;
   psymm_     = orig_psymm;

   idxs       = make_shared<vector<string>>(orig_idxs);
   idx_ranges = make_shared<vector<vector<string>>>(orig_idx_ranges);
   aops       = make_shared<vector<bool>>(orig_aops);
   factor     = orig_factor;
   contracted = false;

   cout << "Generating all unique combinations of index ranges" << name_ << endl;
   generate_ranges( idx_ranges, psymm_);
   
   plus_ops = make_shared<vector<int>>();
   kill_ops = make_shared<vector<int>>();

  for (auto ii =0 ; ii != aops->size() ; ii++) {
    if(aops->at(ii)) {
      plus_ops->push_back(ii);
    } else {
      kill_ops->push_back(ii);
    }
  }

  plus_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(plus_ops->size()+1);
  kill_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(kill_ops->size()+1);
 
  for (int  ii = 0 ; ii != (plus_ops->size() +1) ; ii++) {
    plus_combs->at(ii) =  get_N_in_M_combsX( plus_ops, ii );
    kill_combs->at(ii) =  get_N_in_M_combsX( kill_ops, ii );
  }

  CTP_map = make_shared< map< string, shared_ptr<CtrTensorPart<DType>> >>();

  cout << "done " << endl;
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
  auto prvec = [](shared_ptr<vector<string>> invec){ cout << "[ " ; for (auto elem : *invec) { cout << elem << " " ;} cout << "]" ;};
  auto fvec = make_shared<vector<int>>(idx_ranges->size(), 0); 
  auto maxs = make_shared<vector<int>>();
  int  num_cycles = 1;
  bool tbool = true;
  pair<int,int> fac(1,1);
  all_ranges_ = make_shared<map<shared_ptr<vector<string>>, tuple<bool, shared_ptr<vector<string>>, shared_ptr<vector<string>>, pair<int,int> > >>();

  for(auto range : *idx_ranges){
    maxs->push_back(range.size()-1);
    num_cycles *= range.size();
  }

  //generate all possible ranges
  auto possible_ranges = make_shared<vector<shared_ptr<vector<string>>>>();
  for (int ii = 0 ; ii != num_cycles ; ii++){
    auto new_range = make_shared<vector<string>>(idx_ranges->size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range->at(jj) = idx_ranges->at(jj)[fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges->push_back(new_range);
    fvec_cycle(fvec, maxs);
  } 
  
  //Apply symmetry operations to remove unnecessary ranges   
  orb_ranges = make_shared< map< shared_ptr<vector<string>>, pair<int,int> > >();
  orb_ranges->emplace(possible_ranges->at(0), make_pair(1,1) );
  all_ranges_->emplace(possible_ranges->at(0), tie(tbool, possible_ranges->at(0), idxs, fac ) );

  for (int ii = 1 ; ii!=possible_ranges->size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
       if(apply_symmetry(possible_ranges->at(kk), possible_ranges->at(ii))){
         break;
       }
       if(kk == ii-1){
         orb_ranges->emplace(possible_ranges->at(ii), make_pair(1,1) );
         break;
       }
    } 
  }
cout << "TensOp::generate_ranges" <<   endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool TensOp<DType>::apply_symmetry(shared_ptr<vector<string>> ranges_1, shared_ptr<vector<string>> ranges_2  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "TensOp::apply_symmetry" <<   endl;

  for (auto symmop : symmfuncs_) {
    if( *ranges_1 == *(get<0>(symmop)(ranges_2)) ){
      bool fbool = false; 
      pair<int,int> sfac(get<1>(symmop), get<2>(symmop));
      shared_ptr<vector<string>> transformed_idxs = make_shared<vector<string>>(*(get<0>(symmop)(idxs)));
      all_ranges_->emplace(ranges_2, tie( fbool, ranges_1, transformed_idxs, sfac));
      get<2>(all_ranges_->at(ranges_2)) = transformed_idxs; //WHY DO I NEED TO DO THIS TO MAKE IT WORK? JUST EMPLACING IS NOT ENOUGH!!! VERY SUSPICIOUS!
      return true;
    }
  }

  bool tbool = true; 
  pair<int,int> fac(1,1);
  all_ranges_->emplace(ranges_2, tie(tbool, ranges_2, idxs, fac ) );
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
  auto noctrs = make_shared<vector<pair<int,int>>>(0);
  for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
    auto ReIm_factors = make_shared<vector<pair<int,int>>>(1, get<3>(rng_it->second)); 
    auto CTP = make_shared< CtrTensorPart<DType> >(idxs, rng_it->first, noctrs, ReIm_factors ); 
    cout << "CTP->myname() = " << CTP->myname() << " is going into CTP_map" <<  endl;
    CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs.
  }

  //puts_contracted ranges into map
  for ( int nctrs = 1 ; nctrs != (idxs->size()/2)+1 ; nctrs++ ){
    auto  ctr_lists = get_unique_pairs(plus_ops, kill_ops, nctrs);
    for (auto ctr_vec : *ctr_lists) {
      for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
        bool valid =true;
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first->at(ctr_vec->at(ii).first) != rng_it->first->at(ctr_vec->at(ii).second)){
            valid = false;
            break;
          }
        }

        if (valid) {
          auto ReIm_factors = make_shared<vector<pair<int,int>>>(1, get<3>(rng_it->second)); 
          auto CTP = make_shared< CtrTensorPart<DType> >(idxs, rng_it->first, ctr_vec, ReIm_factors ); 
          cout << "CTP->myname() = " << CTP->myname() << "is going into CTP_map " << endl;
          CTP_map->emplace(CTP->myname(), CTP); //maybe should be addded in with ctr_idxs.
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
//Performs hermitian conjugate on tensor indexes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void TensOp<DType>::hconj() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  reverse(idxs->begin(), idxs->end());

  for(auto range :  *idx_ranges)
    reverse(range.begin(), range.end());

  for (bool op : *aops)
    op =!op;      

  for (auto range_it = orb_ranges->begin(); range_it!=orb_ranges->end(); range_it++ ){
    reverse(range_it->first->begin(), range_it->first->end());
    range_it->second.second *= -1;
  }
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

//Loosely speaking, MultiTensOp has two kinds of members:
//1) Members inherited from the derived class /obtained by merging the members (e.g. the indexes and position
//of creation and annihilation operators) from the TensOps from which the MultiTensOp is constructed. 
//2) "Split" members; vectors where each element corresponds to a constituent TensOp of the MultiTensor.
//   The split members are needed to keep  track of signs resulting from symmetry operations.

  orig_tensors_=orig_tensors;

  split_idxs     = make_shared<vector<shared_ptr<vector<string>>>>(0) ;
  split_plus_ops = make_shared<vector<shared_ptr<vector<int>>>>(0) ;
  split_kill_ops = make_shared<vector<shared_ptr<vector<int>>>>(0) ;
 
  plus_ops = make_shared<vector<int>>(0);
  kill_ops = make_shared<vector<int>>(0);
  idxs     = make_shared<vector<string>>(0);
  aops     = make_shared<vector<bool>>(0);

  Tsizes     = make_shared<vector<int>>(0);
  cmlsizevec = make_shared<vector<int>>(0);
  idx_ranges = make_shared<vector<vector<string>>>(0);
 
  names = make_shared<vector<string>>(0);


  int  xx= 0;
  for (auto ii = 0 ; ii != orig_tensors_.size(); ii++) {
     split_idxs->push_back(orig_tensors_[ii]->idxs);
     split_plus_ops->push_back(orig_tensors_[ii]->plus_ops);
     split_kill_ops->push_back(orig_tensors_[ii]->kill_ops);

    for ( auto elem : *(orig_tensors_[ii]->idxs)){
      idxs->push_back(elem);
    }

    for ( auto elem : *(orig_tensors_[ii]->plus_ops))
      plus_ops->push_back(elem+xx);
    for ( auto elem : *(orig_tensors_[ii]->kill_ops))
      kill_ops->push_back(elem+xx);
    for ( auto elem : *(orig_tensors_[ii]->idx_ranges))
      idx_ranges->push_back(elem);
    for ( auto elem : *(orig_tensors_[ii]->aops))
      aops->push_back(elem);

    names->push_back(orig_tensors_[ii]->name());
    Tsizes->push_back(orig_tensors_[ii]->idxs->size());
    cmlsizevec->push_back(xx);
    xx+=orig_tensors_[ii]->idxs->size();
  } 
  cout << "names = " ; for (auto elem : *names) { cout << elem << " " ; } cout << endl;

  generate_ranges();
 
  plus_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(plus_ops->size()+1);
  kill_combs = make_shared< vector< shared_ptr< vector< shared_ptr<vector<int> > > > > >(kill_ops->size()+1);
 
  for (int  ii = 0 ; ii != (plus_ops->size() +1) ; ii++) {
    plus_combs->at(ii) =  get_N_in_M_combsX( plus_ops, ii );
    kill_combs->at(ii) =  get_N_in_M_combsX( kill_ops, ii );
  }

  CTP_map = make_shared< map< string, shared_ptr<CtrTensorPart<DType>> >>();
  CMTP_map = make_shared< map< string, shared_ptr<CtrMultiTensorPart<DType>> >>();
  cout << "MultiTensOp::initialize endl;" << endl;
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

  auto prvec = [](shared_ptr<vector<string>> invec){ cout << "[ " ; for (auto elem : *invec) { cout << elem << " " ;} cout << "]" ;};

  vector<map< shared_ptr<vector<string>>,  tuple<bool, shared_ptr<vector<string>>, shared_ptr<vector<string>>,  pair<int,int> > >::iterator> rng_maps;
  vector<int> posvec; 
  vector<size_t> sizevec; 

  combined_ranges = make_shared<map< shared_ptr<vector<string>>,  tuple<bool, shared_ptr<vector<string>>,  shared_ptr<vector<string>>, shared_ptr<vector<pair<int,int>>> > >>();
  split_ranges = make_shared<map< shared_ptr<vector<shared_ptr<vector<string>>>>,
                           tuple< shared_ptr<vector<bool>>, shared_ptr<vector<shared_ptr<vector<string>>>>, shared_ptr<vector<pair<int,int>>> > >>();

  for (auto Ten : orig_tensors_) {
    rng_maps.push_back(Ten->all_ranges()->begin());
    sizevec.push_back(Ten->all_ranges()->size());
    posvec.push_back(0);
  } 

  bool ham =true;
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
        rng_maps[ii]++;
        posvec[ii-1] = 0;
        posvec[ii]++;
        continue;
      }     

      auto orig_ranges   = make_shared<vector<shared_ptr<vector<string>>>>(0);
      auto isunique      = make_shared<vector<bool>>(0);
      auto unique_ranges = make_shared<vector<shared_ptr<vector<string>>>>(0);
      auto unique_idxs   = make_shared<vector<shared_ptr<vector<string>>>>(0);
      auto factors       = make_shared<pint_vec>(0);

      for (int ii =0 ; ii != orig_tensors_.size() ; ii++){
         orig_ranges->push_back(rng_maps[ii]->first);
         isunique->push_back(get<0>(rng_maps[ii]->second))  ;
         unique_ranges->push_back(get<1>(rng_maps[ii]->second));  
         unique_idxs->push_back(get<2>(rng_maps[ii]->second));  
         factors->push_back(get<3>(rng_maps[ii]->second)) ;
      }

      split_ranges->emplace(orig_ranges, tie(isunique, unique_ranges, factors));

     auto merged_oranges = make_shared<vector<string>>(0);
      for(auto orange : *orig_ranges ) { 
        merged_oranges->insert(merged_oranges->end(), orange->begin(), orange->end());
      } 

      auto merged_uranges = make_shared<vector<string>>(0);
      for(auto urange : *unique_ranges ){
        merged_uranges->insert(merged_uranges->end(), urange->begin(), urange->end());
      }

      auto merged_uqidxs = make_shared<vector<string>>(0);
      for(auto uqidxs : *unique_idxs ){
        merged_uqidxs->insert(merged_uqidxs->end(), uqidxs->begin(), uqidxs->end());
      } 
  
      bool merged_unique = true;
      if ( *merged_oranges == *merged_uranges )
        merged_unique = false;
      
      combined_ranges->emplace(merged_oranges, tie(merged_unique, merged_uranges, merged_uqidxs, factors));
    }
  } while(ham);
 
  cout <<  "generated_ranges " <<endl;

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
  auto noctrs = make_shared<vector<pair<int,int>>>(0);
  
  //silly, should just test {act,act,....} against contraints, and check if act in each ranges... 
  for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {
    bool check = true;
    for (int xx = 0; xx !=rng_it->first->size() ; xx++ ) {
      if (rng_it->first->at(xx) != "act") {
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
  for ( int nctrs = 1 ; nctrs != (idxs->size()/2)+1 ; nctrs++ ){
    auto  ctr_lists = get_unique_pairs(plus_ops, kill_ops, nctrs);

    for (auto ctr_vec : *ctr_lists) {
      for (auto rng_it = combined_ranges->begin(); rng_it != combined_ranges->end(); rng_it++) {

        bool valid =true;
        //checks ranges for contractions match
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first->at(ctr_vec->at(ii).first) != rng_it->first->at(ctr_vec->at(ii).second)){
            valid = false;
            break;
          }
        }

        if (!valid) 
          continue;

         //checks all uncontracted indexes are active. Should call constraint functions instead
        if (valid){	 
          vector<bool> unc_get(rng_it->first->size(),true);
          for (auto ctr_pos : *ctr_vec){
            unc_get[ctr_pos.first]  = false;
            unc_get[ctr_pos.second] = false;
          }
          for (int xx = 0; xx !=rng_it->first->size() ; xx++ ) {
            if (unc_get[xx]  && (rng_it->first->at(xx).substr(0,3) != "act") ) {
              valid = false;  
              break;
            }
          }
        } 

        if (valid) { 
          enter_into_CMTP_map(*ctr_vec, get<3>(rng_it->second), rng_it->first );
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
  auto ctr_idxs_in   = make_shared<pstr_vec>(0);
  auto ctr_ranges_in = make_shared<pstr_vec>(0);
  auto ctr_spins_in = make_shared<pstr_vec>(0);

  auto unc_idxs_in   = make_shared<vector<shared_ptr<vector<string>>>>(0);
  auto unc_ranges_in = make_shared<vector<shared_ptr<vector<string>>>>(0);
  auto unc_aops_in = make_shared<vector<shared_ptr<vector<bool>>>>(0);

  auto CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart<DType>> >> (orig_tensors_.size()); 
  vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
  vector<vector<pair<int,int>>> sameT_ctrs_pos(orig_tensors_.size(),  pint_vec(0));

 
  //seperate contractions into those on the same tensor, and those between different tensors 
  for ( auto ctr_pos : ctr_pos_list ) {

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
   
//  cout << "orig_tensors names = [ " ; cout.flush(); 
//  for (int ii = 0 ; ii !=orig_tensors_.size() ; ii++ ){ 
//   cout << orig_tensors_[ii]->name() << " " ; cout.flush(); 
//  }
//  cout<< endl;

  //get_ranges for individual tensors
  for (int ii = 0 ; ii !=orig_tensors_.size() ; ii++ ){ 
     auto TS_id_ranges = make_shared<vector<string>>(0);
     for(int jj = cmlsizevec->at(ii) ; jj != cmlsizevec->at(ii)+orig_tensors_[ii]->idxs->size(); jj++){
       TS_id_ranges->push_back(id_ranges->at(jj));
     }
   
     CTP_vec->at(ii) = make_shared< CtrTensorPart<DType> >( orig_tensors_[ii]->idxs, TS_id_ranges,
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
