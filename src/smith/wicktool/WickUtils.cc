#include <bagel_config.h>
#ifdef COMPILE_SMITH
 #include <src/smith/wicktool/WickUtils.h>
 // #include "WickUtils.h"
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

namespace WickUtils {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pint_vec>>  
get_cross_pairs( shared_ptr<vector<int>> vec1 , shared_ptr<vector<int>> vec2, shared_ptr<vector<string>> id_names ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
 if  ( vec1->empty() || vec2->empty()){
   auto all_pairs =  make_shared<vector<vector<pair<int,int>>>>(0);
   return all_pairs;
 }

 auto fvec = make_shared<vector<int>>(vec2->size(),0);
 auto maxs = make_shared<vector<int>>();
 auto mins = make_shared<vector<int>>(vec2->size(),0);

 auto prvec = [](shared_ptr<vector<int>> invec){ cout << "[ " ; for (auto elem : *invec) { cout << elem << " " ;} cout << "]" ;};

 vector<int> reset_pos;
 for (int ii = 0; ii !=vec2->size();  ii++){
   maxs->push_back(vec2->size()-1-ii);
   reset_pos.push_back(ii);
 }
 
 auto all_pairs =  make_shared<vector<vector<pair<int,int>>>>(0);
 do {   
   vector<int> perm;
   vector<int> pos = reset_pos;;

   for(int ii= 0; ii !=fvec->size() ; ii++){
      perm.push_back(vec2->at(pos[fvec->at(ii)]));       
      pos.erase(pos.begin()+fvec->at(ii));
   }

   vector<pair<int,int>> pair_vec(0);
   bool check = true;
   for (auto ii =0 ;  ii !=perm.size(); ii++) {
     if ( id_names->at(vec1->at(ii))[0] == id_names->at(perm[ii])[0] ){ 
       check =false;
       break;
     }
     pair_vec.push_back(make_pair(vec1->at(ii),perm[ii])) ;
   }

   if (check) all_pairs->push_back(pair_vec);

 } while(fvec_cycle(fvec, maxs, mins));

 return all_pairs;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if (forvec->at(ii) == max->at(ii)) {
      forvec->at(ii) = 0;
      continue;
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max , shared_ptr<vector<int>> min) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (forvec->size() == 0 )
    return false;

  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if (forvec->at(ii) == max->at(ii)) {
      if (ii == 0) 
        return false;    
      forvec->at(ii) = min->at(ii);
      continue;
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Cycles vector which mimics nested for-loop, constrained so the iterators in the inner loop are always larger
//than those in the outer group. Be careful; the max vector should reflect these constraints (e.g. max[ii] < max[jj] )
// for (ii= 0 ; ii!=max->at(ii).. ii ++)
//   for (jj =ii+1 ; ...; jj++)
//     for (kk = jj+1 ; ..; kk++)
//      .....
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool constrained_fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int ii = forvec->size()-1;
  if (forvec->at(ii)+1 <= max->at(ii)){
   forvec->at(ii)++; 
   return true;
  } else {

    int kk =0;
    do {
      kk++;
      if (forvec->at(ii-kk)+1 <= max->at(ii-kk)){
        break;
      }
      if (kk == ii ){ 
        return false;}
    } while(true);

   forvec->at(ii-kk)++;
   kk--;
   if (ii-kk == 0 ) {
     if((forvec->at(0)+1 > max->at(0)))
     return false;
   } 
   
   int jj =1;
   while(kk>-1) {
   //   cout <<"forvec->at("<<ii-kk<<") = forvec->at("<<ii-kk-1<<")+1" << endl;; 
   //   cout << "in [ "; for( auto  elem  : *forvec ) cout << elem << " " ; cout << "]"<< endl;
      forvec->at(ii-kk) = forvec->at(ii-kk-1)+1; 
      kk--;;
    }

    return true;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class It_Type>
bool fvec_cycle(shared_ptr<vector<It_Type>> forvec, shared_ptr<vector<It_Type>> max , shared_ptr<vector<It_Type>> min) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (forvec->size() == 0 )
    return false;

  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if (forvec->at(ii) == max->at(ii)) {
      if (ii == 0) 
        return false;    
      forvec->at(ii) = min->at(ii);
      continue;
    } else {
      forvec->at(ii)++;
      break;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generates all possible combinations of the elements in a given vector
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class T1>
shared_ptr<vector<shared_ptr<vector<T1>>>> combgen( shared_ptr<vector<T1>> invec){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "combgen" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////

  auto all_combs = make_shared<vector<shared_ptr<vector<T1>>>>(0);
  if( invec->size() ==0 )
    return all_combs;
 
  auto fvec = make_shared<vector<int>>(invec->size(),0);
  auto maxs = make_shared<vector<int>>(invec->size());
  auto mins = make_shared<vector<int>>(invec->size(),0);
  auto origpos = make_shared<vector<int>>(invec->size());
 
  int kk=0;
  for (int ii = invec->size()-1;  ii !=-1; ii--){
    maxs->at(kk) = ii ;
    origpos->at(kk)= kk;
    kk++;
  }
 
  do {
    auto comb = make_shared<vector<T1>>();
    auto pos = make_shared<vector<int>>(*origpos);
    for (auto jj = 0 ; jj !=fvec->size() ; jj++) {
      comb->push_back(invec->at(pos->at(fvec->at(jj))));
      pos->erase(pos->begin()+fvec->at(jj));
    }
    all_combs->push_back(comb);
  } while (fvec_cycle(fvec, maxs, mins));
 
   return all_combs;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gets all possible combinations of NN elements from vec1 , all combinations will be unique iff all elements of vec1 are unique
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<vector<int>>>> get_N_in_M_combsX( shared_ptr<vector<int>> vec1, int NN ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "get_N_in_M_combsX" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  if (NN == 0 )
    return make_shared<vector<shared_ptr<vector<int>>>>(0);

//  cout <<  "getting " << NN << " elements from " ;
//  cout << "[ ";  for (auto  elem : *vec1) { cout << elem << " " ;} cout << "]"<< endl;

  auto fvec = make_shared<vector<int>>(NN,0);
  auto maxs = make_shared<vector<int>>();
  auto mins = make_shared<vector<int>>(NN,0);
 
  for (int ii = 0; ii !=fvec->size();  ii++)
    maxs->push_back(vec1->size()-1-ii);
  
  vector<int> reset_pos;
  for (int ii =0 ; ii != vec1->size(); ii++)
    reset_pos.push_back(ii);
  
  auto N_in_M_combs =  make_shared<vector<shared_ptr<vector<int>>>>();
  do {   
    auto perm = make_shared<vector<int>>() ;
    vector<int> pos = reset_pos;
    for(int ii= 0; ii !=fvec->size() ; ii++){
       perm->push_back(vec1->at(pos[fvec->at(ii)]));       
       pos.erase(pos.begin()+fvec->at(ii));
    }
    N_in_M_combs->push_back(perm);
  } while (fvec_cycle(fvec, maxs, mins));

  return N_in_M_combs;
} 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gets all possible combinations of NN elements from vec1 , all combinations will be unique iff all elements of vec1 are unique
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DT1>
shared_ptr<vector<shared_ptr<vector<DT1>>>> get_N_in_M_combsX( shared_ptr<vector<DT1>> vec1, int NN ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "get_N_in_M_combsX" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  if (NN == 0 )
    return make_shared<vector<shared_ptr<vector<DT1>>>>(0);


  auto fvec = make_shared<vector<int>>(NN,0);
  auto maxs = make_shared<vector<int>>();
  auto mins = make_shared<vector<int>>(NN,0);
 
  for (int ii = 0; ii !=fvec->size();  ii++)
    maxs->push_back(vec1->size()-1-ii);
  
  vector<int> reset_pos;
  for (int ii =0 ; ii != vec1->size(); ii++)
    reset_pos.push_back(ii);
  
  auto N_in_M_combs =  make_shared<vector<shared_ptr<vector<DT1>>>>();
  cout << " NN = " << NN << endl; 
  do {   
    auto perm = make_shared<vector<DT1>>() ;
    vector<int> pos = reset_pos;
    for(int ii= 0; ii !=fvec->size() ; ii++){
       perm->push_back(vec1->at(pos[fvec->at(ii)]));       
       pos.erase(pos.begin()+fvec->at(ii));
    }
    cout << " perm = " ; for (auto elem : *perm) { cout << elem <<  " " ; } cout << endl;
    N_in_M_combs->push_back(perm);
  } while (fvec_cycle(fvec, maxs, mins));
 
  return N_in_M_combs;
} 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out a vector of pairs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_pvec (pint_vec pvec) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "print_pvec" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  cout << " [ ";
  for (auto elem : pvec)
    cout << "[" << to_string(elem.first)  << "," << to_string(elem.second) << "] ";
  cout << "]"<< endl;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
template<class DT>
void print_vec(vector<DT> invec , string vecname){
///////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "print_vec" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  cout << vecname << " : [ " ;
  for ( auto elem : invec ) cout << elem << " " ;
  cout << "] " << endl;

  return;
} 

////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> get_unique_pairs(shared_ptr<vector<int>> ids1 , shared_ptr<vector<int>> ids2 , int num_pairs){
//////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "get_unique_pairs" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////

  auto pairs_vec = make_shared<vector<shared_ptr<vector<pair<int,int>>>>>(0);

  //stupid fudge
  if(num_pairs == 0){
     return pairs_vec;

  }  else if(num_pairs == 1){

    for (int ii = 0 ; ii !=ids1->size(); ii++){
      for (int jj = 0 ; jj !=ids2->size(); jj++){
        auto pair_comb = make_shared<vector<pair<int,int>>>(0);
        pair_comb->push_back(make_pair(ids1->at(ii), ids2->at(jj)));
        pairs_vec->push_back(pair_comb);
      }
    }
    return pairs_vec;

  } else {

    int half_sz= ids1->size()/2;;
    auto maxs = make_shared<vector<int>>(num_pairs, ids1->size());
    auto fvec = make_shared<vector<int>>(0);
    
    for (int ii = 0; ii != maxs->size();  ii++){
      maxs->at(ii) = ids1->size()-maxs->size() + ii;
      fvec->push_back(ii);
    } 
    
    //get every rearrangement here
    auto id2_combs = get_N_in_M_combsX( ids2, num_pairs );
   
    do {   
      for (auto comb : *id2_combs){
        auto pair_comb = make_shared<vector<pair<int,int>>>(0);
        for (int ii = 0; ii!=fvec->size(); ii++){
          pair_comb->push_back(make_pair(ids1->at(fvec->at(ii)), comb->at(ii)));  
        }
        pairs_vec->push_back(pair_comb); 
      }
    } while(constrained_fvec_cycle(fvec, maxs));
    
  
    return pairs_vec;
  }
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string get_name(shared_ptr<vector<string>> full_idxs, shared_ptr<vector<string>> full_id_ranges,  shared_ptr<pint_vec> all_ctrs_pos) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = "";
  for(auto idxs : *full_idxs)
    name += idxs;// + '-'; 

  name+='_'; 
  for(string id : *full_id_ranges)
    name += id; //+ '-';

  if (all_ctrs_pos->size() !=0){
    name+='_'; 
    for(pair<int,int> ctr : *all_ctrs_pos)
      name += to_string(ctr.first)+"x"+to_string(ctr.second)+ '-';
  }
  name.pop_back(); 
  return name ;
};

/////////////////////////////////////////////////////////////////////////////
string get_name_rdm(shared_ptr<vector<string>> full_idxs, shared_ptr<vector<string>> full_idx_ranges, 
                    shared_ptr<vector<pair<int,int>>> all_deltas_pos ){
/////////////////////////////////////////////////////////////////////////////
  string  name = "";
  for(string idx : *full_idxs)
    name += idx;
  name+="_"; 

  for(string idx_range : *full_idx_ranges)
    name += idx_range[0];

  if (all_deltas_pos->size() !=0){
    name+="_"; 
    for(pair<int,int> delta : *all_deltas_pos)
      name += to_string(delta.first)+to_string(delta.second);
  }
  return name;
};

/////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> get_unc_ids_from_deltas_ids(shared_ptr<vector<int>> ids , shared_ptr<vector<pair<int,int>>> deltas ){
/////////////////////////////////////////////////////////////////////////////
 
   vector<bool> contracted(ids->size(),false);
   for ( pair<int,int> del : *deltas){
     contracted[del.first] = true ;
     contracted[del.second] = true ;
   }

  vector<int> unc_ids(ids->size() - deltas->size() *2);
  
  int ii =0; int jj =0;
   while ( ii != unc_ids.size()){
      if (!contracted[jj]) {
        unc_ids[ii] = ids->at(jj);
        ii++;
      }
      jj++;
   } 

   return make_shared<vector<int>>(unc_ids);
}
/////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> get_unc_ids_from_deltas_ids_comparison(shared_ptr<vector<int>> ids , shared_ptr<vector<pair<int,int>>> deltas ){
/////////////////////////////////////////////////////////////////////////////
 
   vector<int> unc_ids(0);
   for ( int jj =0 ; jj != ids->size(); jj++){
     for ( int ii =0 ; ii != deltas->size(); ii++){
       if ( ids->at(jj) == deltas->at(ii).first || ids->at(jj) == deltas->at(ii).second ) {
         break;
       } else if ( ii == deltas->size()-1 ){
         unc_ids.push_back(ids->at(jj));
       }
     }
   }

//   cout << "unc_ids = [ " ;   for ( int pos : unc_ids) {cout << pos << " " ; } cout << "] " << endl;
   return make_shared<vector<int>>(unc_ids);
}
}//end of namespace:
#endif
