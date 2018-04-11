#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/prop/proptool/proputils.h>
 // #include "wickutils.h"
using namespace std;

  using delta_ints = std::vector<std::vector<std::pair<int,int>>>;
  using delta_strs = std::vector<std::vector<std::pair<std::string,std::string>>>;
  using delta_bools = std::vector<std::vector<std::pair<bool,bool>>>;

  using vv_ints  = std::vector< std::vector<int> >;
  using vv_strs  = std::vector< std::vector<std::string> >;
  using vv_bools = std::vector< std::vector<bool> >;

  using pint_vec = std::vector<std::pair<int,int>>;
  using pstr_vec = std::vector<std::pair<std::string,std::string>>;
  using pbool_vec = std::vector<std::pair<bool,bool>>;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pint_vec>>  
WickUtils::get_cross_pairs( shared_ptr<vector<int>> vec1 , shared_ptr<vector<int>> vec2, shared_ptr<vector<string>> id_names ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "WickUtils::get_cross_pairs" << endl;
  
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
void WickUtils::fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max ) {
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
bool WickUtils::fvec_cycle_test(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if (forvec->at(ii) == max->at(ii)) {
      if ( ii == 0 )
        return false;
      forvec->at(ii) = 0;
      continue;
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max , shared_ptr<vector<int>> min) {
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
bool WickUtils::constrained_fvec_cycle(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max) {
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
bool WickUtils::fvec_cycle(shared_ptr<vector<It_Type>> forvec, shared_ptr<vector<It_Type>> max , shared_ptr<vector<It_Type>> min) {
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
shared_ptr<vector<shared_ptr<vector<T1>>>> WickUtils::combgen( shared_ptr<vector<T1>> invec){
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
shared_ptr<vector<shared_ptr<vector<int>>>> WickUtils::get_N_in_M_combsX( shared_ptr<vector<int>> vec1, int NN ){
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
shared_ptr<vector<shared_ptr<vector<int>>>> WickUtils::get_N_in_M_combsX( shared_ptr< const vector<int>> vec1, int NN ){
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
template<class DataType>
shared_ptr<vector<shared_ptr<vector<DataType>>>> WickUtils::get_N_in_M_combsX( shared_ptr< const vector<DataType>> vec1, int NN ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "get_N_in_M_combsX" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  if (NN == 0 )
    return make_shared<vector<shared_ptr<vector<DataType>>>>(0);

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>(NN,0);
  shared_ptr<vector<int>> mins = make_shared<vector<int>>(NN,0);
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>(NN,0);
 
  for (int ii = 0; ii !=fvec->size();  ii++)
    maxs->at(ii) = vec1->size()-1-ii;
  
  vector<int> reset_pos; 
  for (int ii =0 ; ii != vec1->size(); ii++)
    reset_pos.push_back(ii);
  
  auto N_in_M_combs =  make_shared<vector<shared_ptr<vector<DataType>>>>();
  cout << " NN = " << NN << endl; 
  do {   
    auto perm = make_shared<vector<DataType>>() ;
    vector<int> pos = reset_pos;
    for(int ii= 0; ii !=fvec->size() ; ii++){
       perm->push_back(vec1->at(pos[fvec->at(ii)]));       
       pos.erase(pos.begin()+fvec->at(ii));
    }
//    cout << " perm = " ; for (auto elem : *perm) { cout << elem <<  " " ; } cout << endl;
    N_in_M_combs->push_back(perm);
  } while (fvec_cycle(fvec, maxs, mins));
 
  return N_in_M_combs;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gets all possible combinations of NN elements from vec1 , all combinations will be unique iff all elements of vec1 are unique
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<vector<shared_ptr<vector<DataType>>>> WickUtils::get_N_in_M_combsX( shared_ptr<vector<DataType>> vec1, int NN ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined DBG_WickUtils || defined DBG_all 
cout << "get_N_in_M_combsX" << endl;
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
  if (NN == 0 )
    return make_shared<vector<shared_ptr<vector<DataType>>>>(0);


  auto fvec = make_shared<vector<int>>(NN,0);
  auto maxs = make_shared<vector<int>>();
  auto mins = make_shared<vector<int>>(NN,0);
 
  for (int ii = 0; ii !=fvec->size();  ii++)
    maxs->push_back(vec1->size()-1-ii);
  
  vector<int> reset_pos;
  for (int ii =0 ; ii != vec1->size(); ii++)
    reset_pos.push_back(ii);
  
  auto N_in_M_combs =  make_shared<vector<shared_ptr<vector<DataType>>>>();
  cout << " NN = " << NN << endl; 
  do {   
    auto perm = make_shared<vector<DataType>>() ;
    vector<int> pos = reset_pos;
    for(int ii= 0; ii !=fvec->size() ; ii++){
       perm->push_back(vec1->at(pos[fvec->at(ii)]));       
       pos.erase(pos.begin()+fvec->at(ii));
    }
//    cout << " perm = " ; for (auto elem : *perm) { cout << elem <<  " " ; } cout << endl;
    N_in_M_combs->push_back(perm);
  } while (fvec_cycle(fvec, maxs, mins));
 
  return N_in_M_combs;
} 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Prints out a vector of pairs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WickUtils::print_pvec (pint_vec pvec) {
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

////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>
WickUtils::get_unique_pairs(shared_ptr<const vector<int>> ids1 , shared_ptr<const vector<int>> ids2 , int num_pairs){
//////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> pairs_vec = make_shared<vector<shared_ptr<vector<pair<int,int>>>>>(0);

  //stupid fudge
  if(num_pairs == 0){
     return pairs_vec;

  }  else if(num_pairs == 1){

    for (int ii = 0 ; ii !=ids1->size(); ii++){
      for (int jj = 0 ; jj !=ids2->size(); jj++){
        shared_ptr<vector<pair<int,int>>> pair_comb = make_shared<vector<pair<int,int>>>(0);
        pair_comb->push_back(make_pair(ids1->at(ii), ids2->at(jj)));
        pairs_vec->push_back(pair_comb);
      }
    }
    return pairs_vec;

  } else {

    int half_sz= ids1->size()/2;;
    shared_ptr<vector<int>> maxs = make_shared<vector<int>>(num_pairs, ids1->size());
    shared_ptr<vector<int>> fvec = make_shared<vector<int>>(0);
    
    for (int ii = 0; ii != maxs->size();  ii++){
      maxs->at(ii) = ids1->size()-maxs->size() + ii;
      fvec->push_back(ii);
    } 
    
    //get every rearrangement here
    auto id2_combs = get_N_in_M_combsX( ids2, num_pairs );
   
    do {   
      for (auto comb : *id2_combs){
        shared_ptr<vector<pair<int,int>>> pair_comb = make_shared<vector<pair<int,int>>>(0);
        for (int ii = 0; ii!=fvec->size(); ii++){
          pair_comb->push_back(make_pair(ids1->at(fvec->at(ii)), comb->at(ii)));  
        }
        pairs_vec->push_back(pair_comb); 
      }
    } while(constrained_fvec_cycle(fvec, maxs));
    
  
    return pairs_vec;
  }
} 
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<int,int>>>>>
WickUtils::get_unique_pairs(std::shared_ptr<std::vector<int>> ids1 , std::shared_ptr<std::vector<int>> ids2 ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////
    return get_unique_pairs( ids1, ids2 , (ids1->size() < ids2->size() ? ids1->size() : ids2->size() ) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>
WickUtils::get_unique_pairs(shared_ptr<vector<int>> ids1 , shared_ptr<vector<int>> ids2 , int num_pairs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "WickUtils::get_unique_pairs" << endl;

//  assert( !(ids2->size() < num_pairs) ); 
//  assert( !(ids1->size() < num_pairs) ); 

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

/////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> WickUtils::get_unc_ids_from_deltas_ids(shared_ptr<vector<int>> ids , shared_ptr<vector<pair<int,int>>> deltas ){
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
shared_ptr<vector<int>> WickUtils::get_unc_ids_from_deltas_ids_comparison(shared_ptr<vector<int>> ids , shared_ptr<vector<pair<int,int>>> deltas ){
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//same as  fvec_cycle, but allows skipping. Should be included everywhere to guard against max==min problem};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max , shared_ptr<vector<int>> min ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if ( max->at(ii) == min->at(ii) ) {
      if ( ii == 0 )
        return false;
    } else if (forvec->at(ii) == max->at(ii)) {
      if (ii == 0) 
        return false;    
      forvec->at(ii) = min->at(ii);
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//same as normal fvec_cycle, but will not iterate the index at fixed_index_position
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max, int fixed_index_position ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if ( ii == fixed_index_position ) {
      if ( ii == 0 )
        return false;
    } else if (forvec->at(ii) == max->at(ii)) {
      if (ii == 0) 
        return false;    
      forvec->at(ii) = 0;
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper_f2b(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max , shared_ptr<vector<int>> min ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = 0; ii!=forvec->size(); ii++) {
    if ( max->at(ii) == min->at(ii) ) {
      if ( ii == forvec->size()-1 )
        return false;
    } else if (forvec->at(ii) == max->at(ii)) {
      if (ii == forvec->size()-1) 
        return false;    
      forvec->at(ii) = min->at(ii);
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> WickUtils::reorder_vector(vector<int>& neworder , const vector<int>& origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> newvec(origvec.size());
  vector<int>::iterator newvec_it = newvec.begin();

  for( int pos : neworder )
     *newvec_it++ = origvec[pos];

  return make_shared<vector<int>>(newvec);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_civec_name(const int state_num, const int norb, const int nalpha, const int nbeta)  { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = to_string(state_num) + "_["+ to_string(norb)+"o{" + to_string(nalpha) + "a," + to_string(nbeta) + "b}]" ;
  return name ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_gamma_name( shared_ptr<const vector<string>> full_idx_ranges,  shared_ptr<const vector<bool>> aops_vec,
                                  shared_ptr<vector<int>> idxs_pos, string Bra_name, string Ket_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "WickUtils::get_gamma_name" << endl; 
  string  name;
 
  if (idxs_pos->size() == 0 ) {
     name = "ID" ;
  } else {  
    name = "<" + Bra_name + "|_(";
    for (int pos : *idxs_pos ) 
      name+=full_idx_ranges->at(pos)[0];
    
    name+='_';
    for (int pos : *idxs_pos ) {
      if(aops_vec->at(pos)){ 
        name += '1';
      } else {
        name += '0';
      }
    } 
    name += ")_|" + Ket_name + ">";
  }
  
  return name;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_gamma_name( const vector<string>& full_idx_ranges,  const vector<bool>& aops_vec,
                                  const vector<int>& idxs_pos, string Bra_name, string Ket_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "WickUtils::get_gamma_name" << endl; 
  string  name;
 
  if (idxs_pos.size() == 0 ) {
     name = "ID" ;
  } else {  
    name = "<" + Bra_name + "|_(";
    for (int pos : idxs_pos ) 
      name+=full_idx_ranges[pos][0];
    
    name+='_';
    for (int pos : idxs_pos ) {
      if(aops_vec[pos]){ 
        name += '1';
      } else {
        name += '0';
      }
    } 
    name += ")_|" + Ket_name + ">";
  }
  
  return name;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_gamma_name( shared_ptr<vector<string>> full_idx_ranges,  shared_ptr<vector<bool>> aops_vec,
                                  shared_ptr<vector<int>> idxs_pos, string Bra_name, string Ket_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_WickUtils 
cout << "WickUtils::get_gamma_name" << endl; 
#endif 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string  name;
 
  if (idxs_pos->size() == 0 ) {
     name = "ID" ;
  } else {  
    name = "<" + Bra_name + "|_(";
    for (int pos : *idxs_pos ) 
      name+=full_idx_ranges->at(pos)[0];
    
    name+='_';
    for (int pos : *idxs_pos ) {
      if(aops_vec->at(pos)){ 
        name += '1';
      } else {
        name += '0';
      }
    } 
    name += ")_|" + Ket_name + ">";
  }
  
  return name;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Check if the range block will survive without contraction with another
// range block (i.e. see if as many particles are being created in a given
// range as are destroyed) 
////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::RangeCheck(const vector<string>& id_ranges, const vector<bool>& aops ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::RangeCheck" << endl;

  vector<string> diff_rngs(1, id_ranges[0] );
  vector<int> updown(1, ( aops[0] ? 1 : -1)  );

  for ( int jj = 1;  jj !=id_ranges.size() ; jj++){
    int ii = 0;

    string rng = id_ranges[jj];

    do {
      if(rng == diff_rngs[ii]){
        if ( aops[jj]){
          updown[ii]+=1;
        } else {         
          updown[ii]-=1;
        }
        break;
      }
      if ( ii == diff_rngs.size()-1){
        diff_rngs.push_back(rng);
        if ( aops[jj]){
          updown.push_back(1);  
        } else {
          updown.push_back(-1); 

        }
        break;
      }

      ii++;
    } while (true);
  } 

  for (int ac : updown ) 
    if (ac != 0 )
      return false;
 
  return true;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get range for ctrtensorpart  ( uncontracted version ) 
////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_ctp_name( const vector<string>& idxs, const vector<string>& id_ranges ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::get_ctp_name" << endl; 

  string ctp_name = ""; 
  for ( string id : idxs ) 
    ctp_name += id;
 
  ctp_name += "_";

  for (string id : id_ranges ) 
    ctp_name += id[0];

  return ctp_name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get range for ctrtensorpart  
////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_ctp_name( const vector<string>& idxs, const vector<string>& id_ranges, 
                                const vector<pair<int,int>>& ctrs_pos ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::get_ctp_name" << endl; 

  string ctp_name = ""; 
  for ( string id : idxs ) 
    ctp_name += id;
 
  ctp_name += "_";

  for (string id : id_ranges ) 
    ctp_name += id[0];

  if (ctrs_pos.size() != 0 ){ 
    ctp_name += "_";
    vector<pair<int,int>> ctrs_buff =  ctrs_pos;
    shared_ptr<vector<pair<int,int>>> ctrs_standard = standardize_delta_ordering_generic(ctrs_buff, idxs); 
    for ( pair<int,int>& ctr : *ctrs_standard )  
      ctp_name += to_string(ctr.first)+to_string(ctr.second);
  }
 
  return ctp_name;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<int,int>>>
WickUtils::standardize_delta_ordering_generic( const vector<pair<int,int>>& deltas_pos, const vector<string>& idxs  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::Standardize_delta_ordering_generic" << endl;
//TODO must order by indexes, not just by initial position
//     If one of the indexes is X, cannot "just" not contract
//     must also account for reordering ; T_{ijkl} = ... + <I| ijklmnop | J> A_{mnop} delta_{lm}
//     T_{ijkl} = .... + < I | ijklmnop | J> A_{lmop} 
//     so must flip order of A based on T (X) position

 shared_ptr<vector<pair<int,int>>> new_deltas_pos;

 if (deltas_pos.size() > 1 ) {
   new_deltas_pos =  make_shared <vector<pair<int,int>>>( deltas_pos.size());
   vector<int> posvec(deltas_pos.size(),0);
   for (int ii = 0 ; ii != deltas_pos.size() ; ii++)
     for (int jj = 0 ; jj != deltas_pos.size() ; jj++)
       if (idxs[deltas_pos[ii].first] > idxs[deltas_pos[jj].first] )
         posvec[ii]++;

   for (int ii = 0 ; ii != deltas_pos.size() ; ii++)
     new_deltas_pos->at(posvec[ii]) = deltas_pos[ii];

 } else {
   vector<pair<int,int>> new_deltas_pos_raw(deltas_pos);
   new_deltas_pos = make_shared<vector<pair<int,int>>>( new_deltas_pos_raw ) ;

 }
  return new_deltas_pos;
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::range_to_prime(char range ) {
//////////////////////////////////////////////////////////////////////////////////
  switch (range){
    case 'v' : return 2;
    case 'V' : return 3;
    case 'c' : return 5;
    case 'C' : return 7;
    case 'a' : return 11;
    case 'A' : return 13;
    default : 
      throw std::logic_error( " unknown range " ); return 9999999;
  }
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::range_to_prime_spinfree(char range ) {
//////////////////////////////////////////////////////////////////////////////////
  switch (range){
    case 'v' : return 2;
    case 'c' : return 3;
    case 'a' : return 5;
    default : 
      throw std::logic_error( " unknown range " ); return 9999999;
  }
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::get_block_hash( const std::vector<std::string>& block  ) {
//////////////////////////////////////////////////////////////////////////////////
 
  int max_range_num = 13; // should set this in main proptool
  int range_id = 0;
  int ii = 0;
  for ( vector<string>::const_reverse_iterator b_it = block.crbegin() ; b_it != block.crend() ;  b_it++, ii++ ) 
     range_id += range_to_prime( (*b_it)[0] ) * pow( max_range_num, ii ); 
 
  assert( range_id >  -1 );

  return range_id; 
}
//////////////////////////////////////////////////////////////////////////////////////
vector<char> WickUtils::strvec_to_chrvec( const vector<string>& strvec ) {
//////////////////////////////////////////////////////////////////////////////////////
 
   vector<char> chrvec(strvec.size());
   vector<char>::iterator cv_it = chrvec.begin();
   for( vector<string>::const_iterator sv_it = strvec.begin(); sv_it != strvec.end(); sv_it++, cv_it++ )
     *cv_it = (*sv_it)[0];
  
   return chrvec  ; 
} 
//////////////////////////////////////////////////////////////////////////////////////
vector<string> WickUtils::chrvec_to_strvec( const vector<char>& chrvec ) {
//////////////////////////////////////////////////////////////////////////////////////
 
  vector<string> strvec(chrvec.size());
  vector<char>::const_iterator cv_it = chrvec.begin();
  for ( vector<string>::iterator sv_it = strvec.begin(); sv_it != strvec.end(); sv_it++, cv_it++ ) {
    string s = "";
    s += *cv_it;
    *sv_it = s;
  }

  return strvec;
} 
#endif
