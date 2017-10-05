#include <bagel_config.h>
#ifdef COMPILE_SMITH

 #include <src/smith/wicktool/CtrTensOp.h>
 #include <src/smith/wicktool/WickUtils.h>
 #include <src/smith/wicktool/gamma_generator.h>

 //#include "WickUtils.h"
 //#include "CtrTensOp.h"
 //#include "gamma_generator.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////
template<typename DType>
void TensorPart<DType>::get_name(){
////////////////////////////////////////////////////////////////////////////
 name = "";

 for(std::string id : *idxs)
   name += id;
 name+="_"; 

 for(std::string id : *id_ranges)
   name += id[0];

 return;
};

/////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrTensorPart<DType>::get_name(){
/////////////////////////////////////////////////////////////////////////////
  name = "";
  for(string id : *full_idxs)
    name += id;
  name+="_"; 

  for(string id : *full_id_ranges)
    name += id[0];

  auto ctrs_buff = make_shared<vector<pair<int,int>>>(*ctrs_pos);
  auto ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

  if (ctrs_buff_standard->size() !=0){
    name+="_"; 
    for(pair<int,int> ctr : *ctrs_buff_standard)
      name += to_string(ctr.first)+to_string(ctr.second);
  }
  return;
};

/////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrMultiTensorPart<DType>::get_name(){
/////////////////////////////////////////////////////////////////////////////
  name = "";
  for(string id : *full_idxs)
    name += id;
  name+="_"; 

  for(string id : *full_id_ranges)
    name += id[0];

  auto ctrs_buff = make_shared<vector<pair<int,int>>>(*all_ctrs_pos);
  auto ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

  if (ctrs_buff_standard->size() !=0){
    name+="_"; 
    for(pair<int,int> ctr : *ctrs_buff_standard)
      name += to_string(ctr.first)+to_string(ctr.second);
  }
  return;
};

/////////////////////////////////////////////////////////////////////////////
template<class DType>
string CtrTensorPart<DType>::get_next_name(shared_ptr<vector<pair<int,int>>> new_ctrs_pos){
/////////////////////////////////////////////////////////////////////////////
  string new_name = "";
  for(string id : *full_idxs)
    new_name += id;
  new_name+="_"; 

  for(string id : *full_id_ranges)
    new_name += id[0];

  if (new_ctrs_pos->size() >=1){
    new_name+="_"; 
    for(int ii=0;  ii!= new_ctrs_pos->size(); ii++){
      new_name += to_string(new_ctrs_pos->at(ii).first)+to_string(new_ctrs_pos->at(ii).second);
    }
  }

  return new_name;
}
/////////////////////////////////////////////////////////////////////////////
template<class DType>
string CtrMultiTensorPart<DType>::get_next_name(shared_ptr<vector<pair<int,int>>> new_ctrs_pos){
/////////////////////////////////////////////////////////////////////////////
  string new_name = "";
  for(string id : *full_idxs)
    new_name += id;
  new_name+="_"; 

  for(string id : *full_id_ranges)
    new_name += id[0];

  if (new_ctrs_pos->size() >=1){
    new_name+="_"; 
    for(int ii=0;  ii!= new_ctrs_pos->size(); ii++){
      new_name += to_string(new_ctrs_pos->at(ii).first)+to_string(new_ctrs_pos->at(ii).second);
    }
  }

  return new_name;
};
//////////////////////////////////////////////////////////////////////////////
template<typename DType>
void CtrTensorPart<DType>::get_ctp_idxs_ranges(){
//////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::get_ctp_idxs_ranges" << endl; 
#endif 
//////////////////////////////////////////////////////////////////////////////

  vector<bool> get_unc(full_idxs->size(), true);
  for (int ii =0; ii<ctrs_pos->size() ; ii++){
    if ( ctrs_pos->at(ii).first == ctrs_pos->at(ii).second)
      break;
    get_unc[ctrs_pos->at(ii).first] = false;
    get_unc[ctrs_pos->at(ii).second] = false;
   }

  bool survive_indep = true;
  unc_pos = make_shared<vector<int>>(0);
  id_ranges = make_shared<vector<string>>(0);
  idxs = make_shared<vector<string>>(0);
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      id_ranges->push_back(full_id_ranges->at(ii));
      idxs->push_back(full_idxs->at(ii));
      unc_pos->push_back(ii);
    }
  } 
  
  unc_rel_pos = make_shared<map<int,int>>();
  for( int ii =0 ; ii != unc_pos->size(); ii++) 
    unc_rel_pos->emplace(unc_pos->at(ii), ii);

//  vprint(*unc_pos, "unc_pos"); 
//  vprint(*idxs, "idxs"); 
//  vprint(*id_ranges, "id_ranges"); 
 
  return; 
}
//////////////////////////////////////////////////////////////////////////////
template<typename DType>
void CtrMultiTensorPart<DType>::get_cmtp_idxs_ranges(){
//////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::get_ctp_idxs_ranges" << endl; 
#endif 
//////////////////////////////////////////////////////////////////////////////

  vector<bool> get_unc(full_idxs->size(), true);
  for (int ii =0; ii<all_ctrs_pos->size() ; ii++){
    if ( all_ctrs_pos->at(ii).first == all_ctrs_pos->at(ii).second)
      break;
    get_unc[all_ctrs_pos->at(ii).first] = false;
    get_unc[all_ctrs_pos->at(ii).second] = false;
   }

  bool survive_indep = true;
  all_unc_pos = make_shared<vector<int>>(0);
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      all_unc_pos->push_back(ii);
    }
  } 
  
  all_unc_rel_pos = make_shared<map<int,int>>();
  for( int ii =0 ; ii != all_unc_pos->size(); ii++) 
    all_unc_rel_pos->emplace(all_unc_pos->at(ii), ii);

//  vprint(*unc_pos, "unc_pos"); 
//  vprint(*idxs, "idxs"); 
//  vprint(*id_ranges, "id_ranges"); 
 
  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrTensorPart<DType>::FullContract(shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                        shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list ){
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::FullContract" << endl; 
#endif 
/////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ctp_loc = Tmap->find(name);
  
 // search for tensor in map if not there, search for ones with fewer contraction (rubbish at present; involves repetitive searching)
  if ( (ctp_loc != Tmap->end())  ){ 
    contracted = true;
    return;    
  }

  if(ctp_loc == Tmap->end()){
    if(ctrs_pos->size() == 0) { 
      auto new_ctrs_pos = make_shared<vector<pair<int,int>>>(0);
      auto new_ReIm_factors = make_shared<vector<pair<int,int>>>(0); 
      auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, new_ctrs_pos, new_ReIm_factors ); 
      string new_name = get_next_name(new_ctrs_pos);
      Tmap->emplace(new_name, new_CTP);
    }
   

    if(ctrs_pos->size() > 1) { 
      auto new_ctrs_pos = make_shared<vector<pair<int,int>>>(ctrs_pos->begin(), ctrs_pos->end()-1);
      string new_name = get_next_name(new_ctrs_pos);
      auto new_ctp_loc = Tmap->find(new_name);
      
      // if find one less, contract it.
      if ( new_ctp_loc != Tmap->end() ){ 
        pair<int,int> ctrs_rel_pos = make_pair(new_ctp_loc->second->unc_rel_pos->at(ctrs_pos->back().first), new_ctp_loc->second->unc_rel_pos->at(ctrs_pos->back().second));
        ACompute_list->push_back(make_shared<CtrOp_same_T> (new_name, name, ctrs_pos->back(), ctrs_rel_pos, "same_T D" ));
        Tmap->emplace(name , make_shared<CtrTensorPart>(*this));

      } else { 
        auto new_ReIm_factors = make_shared<vector<pair<int,int>>>(ReIm_factors->begin(), ReIm_factors->end()-1); 
        auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, new_ctrs_pos, new_ReIm_factors ); 
        new_CTP->FullContract(Tmap,  ACompute_list);
        Tmap->emplace(name , make_shared<CtrTensorPart>(*this));
        pair<int,int> ctrs_rel_pos = make_pair(new_CTP->unc_rel_pos->at(ctrs_pos->back().first), new_CTP->unc_rel_pos->at(ctrs_pos->back().second));
        ACompute_list->push_back(make_shared<CtrOp_same_T> (new_name, name, ctrs_pos->back(), ctrs_rel_pos, "same_T C" ));
      }

    } else if(ctrs_pos->size() == 1 ){ 
      auto unc_ctrs_pos = make_shared<vector<pair<int,int>>>(1, make_pair(0,0) );
      auto unc_ReIm_factors = make_shared<vector<pair<int,int>>>(0); 
      auto unc_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, unc_ctrs_pos, unc_ReIm_factors ); 

      // Need to create a list of uncontracted tensors
      unc_CTP->CTdata = make_shared<DType>(); 
      Tmap->emplace(unc_CTP->name, unc_CTP);
      required_Tblocks->push_back(unc_CTP->name);
      pair<int,int> ctrs_rel_pos = make_pair(unc_CTP->unc_rel_pos->at(ctrs_pos->back().first), unc_CTP->unc_rel_pos->at(ctrs_pos->back().second));
      ACompute_list->push_back(make_shared<CtrOp_same_T> (unc_CTP->name, name, ctrs_pos->back(), ctrs_rel_pos, "same_T B" ));
      Tmap->emplace(name, make_shared<CtrTensorPart>(*this));
    }
  }
  contracted = true;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Very similar to single tensor version; only difference in case where nothing relevant tensors are found.
// This only works for Two tensors, must test for irreducibility when >2 tensors....
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrMultiTensorPart<DType>::FullContract(shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                             shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::FullContract" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "CtrMultiTensorPart<DType>::FullContract" << endl; 
  auto ctp_loc = Tmap->find(name);

  if (ctp_loc != Tmap->end() ) {
    contracted = true;
    return;    
  }

  if((ctp_loc == Tmap->end())){ ;
    if(cross_ctrs_pos->size() > 1) {  

      auto new_ReIm_factors = make_shared<vector<pair<int,int>>>(ReIm_factors->begin(), ReIm_factors->end()-1);
      auto new_all_ctrs_pos = make_shared<vector<pair<int,int>>>(all_ctrs_pos->begin(), all_ctrs_pos->end()-1);
      auto new_cross_ctrs_pos = make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(cross_ctrs_pos->begin(), cross_ctrs_pos->end()-1);

      string new_name = get_next_name(new_all_ctrs_pos);
      auto new_ctp_loc = Tmap->find(new_name);
  
      // if find one less, contract it.
      if( new_ctp_loc != Tmap->end()  ){  
        auto new_CTP = make_shared<CtrTensorPart<DType>>(full_idxs, full_id_ranges, all_ctrs_pos, ReIm_factors ); 
        pair<int,int> ctr_rel_pos = make_pair( all_ctrs_pos->back().first, all_ctrs_pos->back().second );
        ACompute_list->push_back(make_shared<CtrOp_same_T> (new_name, new_CTP->name, all_ctrs_pos->back(), ctr_rel_pos, "same_T A" ));
      
      } else {  
        auto new_CMTP = make_shared<CtrMultiTensorPart<DType>>(CTP_vec, new_cross_ctrs_pos); 
        new_CMTP->FullContract(Tmap, ACompute_list);
//        auto new_CTP = make_shared<CtrTensorPart<DType>>(new_CMTP->full_idxs, new_CMTP->full_id_ranges, all_ctrs_pos, ReIm_factors ); 
//        new_CTP->FullContract(Tmap,  ACompute_list);
      }
     
    } else if(cross_ctrs_pos->size() == 1 ){
        int T1loc, T2loc;
        T1loc = cross_ctrs_pos->back().first.first;
        T2loc = cross_ctrs_pos->back().second.first;
 
      if ( CTP_vec->size() == 2) {

        for (auto CTP : *CTP_vec) 
          CTP->FullContract(Tmap, ACompute_list);
          
        //auto new_CTP = Binary_Contract_diff_tensors(CTP_vec->at(T1loc)->myname(), CTP_vec->at(T2loc)->myname(), all_ctrs_pos->back(), Tmap,  ACompute_list);
        auto new_CTP = Binary_Contract_diff_tensors(cross_ctrs_pos->back(), all_ctrs_pos->back(), Tmap,  ACompute_list);
        Tmap->emplace(new_CTP->name, new_CTP);

      } else {
        cout << "USING MT BINARY CONTRACT DIFF TENSORS" << endl;
        auto new_CMTP = Binary_Contract_diff_tensors_MT(CTP_vec->at(cross_ctrs_pos->back().first.first)->myname(),CTP_vec->at(cross_ctrs_pos->back().second.first)->myname(),
                                                        all_ctrs_pos->back(), Tmap, ACompute_list);
        new_CMTP->FullContract(Tmap, ACompute_list);
      }
    }
  }
  contracted = true;
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<CtrMultiTensorPart<DType>>
CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT(string T1name, string T2name, pair<int,int> ctr_todo,
                                                           shared_ptr< map<string, shared_ptr<CtrTensorPart<DType>>> > Tmap,
                                                           shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT" << endl; 
  

   auto T1T2_ctrd =  Binary_Contract_diff_tensors(cross_ctrs_pos->back(), ctr_todo, Tmap, ACompute_list );
   Tmap->emplace(T1T2_ctrd->name, T1T2_ctrd);

   auto new_CTP_vec = make_shared<vector<shared_ptr<CtrTensorPart<DType>>>>(0);
   for ( auto CTP : *CTP_vec ) {  
     if (CTP->name == T1name || CTP->name == T2name)
       continue;
     new_CTP_vec->push_back(CTP);
   }
   
   CTP_vec->push_back(T1T2_ctrd);
   auto new_cross_ctrs_pos = make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(cross_ctrs_pos->begin(), cross_ctrs_pos->end()-1);

   auto new_CMTP = make_shared< CtrMultiTensorPart<DType> >(new_CTP_vec, new_cross_ctrs_pos); 

   return new_CMTP;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<CtrTensorPart<DType>>
 CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors( pair<pair<int,int>, pair<int,int>> cross_ctr,// These two arguments should be 
                                                          pair<int,int> ctr_todo,                      // equivalent ! 
                                                          shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                                          shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors" << endl; 
   cout << "name = " << name << endl; 

   shared_ptr<CtrTensorPart<DType>> T1;
   shared_ptr<CtrTensorPart<DType>> T2;

   int T1pos = cross_ctr.first.first;
   int T2pos = cross_ctr.second.first;
   int T1ctr = cross_ctr.first.second;
   int T2ctr = cross_ctr.second.second;

   //Swapping tensors round to maintain consistent ordering
   if (cross_ctr.first.first < cross_ctr.second.first) {
     T1 = CTP_vec->at(cross_ctr.first.first);
     T2 = CTP_vec->at(cross_ctr.second.first);
   } else { 
     T2 = CTP_vec->at(cross_ctr.first.first);
     T1 = CTP_vec->at(cross_ctr.second.first);
     int buff = T1pos;
     T1pos = T2pos;
     T2pos = buff;
     buff = T1ctr ;
     T1ctr = T2ctr;
     T2ctr = buff;
   }
   string T1name = T1->name; 
   string T2name = T2->name; 

   cout << "T1name = "<< T1name << "    "  << "T1pos = "<< T1pos << "    " << "T1ctr = "<< T1ctr << endl;
   cout << "T2name = "<< T2name << "    "  << "T2pos = "<< T2pos << "    " << "T2ctr = "<< T2ctr << endl;

   auto full_id_ranges = make_shared<vector<string>>(T1->full_id_ranges->begin(), T1->full_id_ranges->end()) ;
   full_id_ranges->insert(full_id_ranges->end(), T2->full_id_ranges->begin(), T2->full_id_ranges->end()); 

   auto full_idxs = make_shared<vector<string>>(T1->full_idxs->begin(), T1->full_idxs->end()) ;
   full_idxs->insert(full_idxs->end(), T2->full_idxs->begin(), T2->full_idxs->end()) ;

   auto full_ctrs = make_shared<vector<pair<int,int>>>(0);
   cout << " getting cml pos " << endl;
   cout << " T1shift =  Tsizes_cml->at("<<T1pos<<")"<< endl;
   int T1shift =  Tsizes_cml->at(T1pos);
   cout << " T2shift =  Tsizes_cml->at("<<T2pos<<")"<< endl;
   int T2shift =  Tsizes_cml->at(T2pos);
   cout << " got cmlpos"<< endl;

   for (auto ctr : *T1->ctrs_pos)
      full_ctrs->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));

   for (auto ctr : *T2->ctrs_pos)
      full_ctrs->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));
      
   cout << "("<< T1ctr+T1shift <<","<< T2ctr+T2shift << ")  ?= (" << ctr_todo.first <<"," << ctr_todo.second << ")" << endl; 

   full_ctrs->push_back( make_pair(T1ctr+T1shift, T2ctr+T2shift) );

   cout << " getting relative positions " ;
   cout << "T1->unc_rel_pos->at("<< T1ctr<<") " << endl;
   int T1_ctr_rel_pos = T1->unc_rel_pos->at(T1ctr);
   cout << "T2->unc_rel_pos->at("<< T2ctr<<") " << endl;
   int T2_ctr_rel_pos = T2->unc_rel_pos->at(T2ctr);

   auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, full_ctrs, make_shared<vector<pair<int,int>>>(1, ctr_todo )); 
   ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, new_CTP->name,  ctr_todo.first, ctr_todo.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
   auto new_CTData = make_shared<DType>();
   new_CTP->contracted = true;

   cout << "BCDT contracting " << T1name << " and " << T2name << " over (" << ctr_todo.first << "," << ctr_todo.second << ") to get " << new_CTP->name << endl;
   return new_CTP;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensorPart<double>;
template class CtrTensorPart<double>;
template class CtrMultiTensorPart<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
    
#endif
