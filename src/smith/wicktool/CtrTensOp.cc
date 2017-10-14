#include <bagel_config.h>
#ifdef COMPILE_SMITH

 #include <src/smith/wicktool/CtrTensOp.h>
 #include <src/smith/wicktool/WickUtils.h>
 #include <src/smith/wicktool/gamma_generator.h>

 //#include "WickUtils.h"
 //#include "CtrTensOp.h"
 //#include "gamma_generator.h"

using namespace std;

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
string CtrTensorPart<DType>::get_next_name(shared_ptr<vector<pair<int,int>>> new_ctrs_pos){
/////////////////////////////////////////////////////////////////////////////
  string new_name = "";
  for(string id : *full_idxs)
    new_name += id;
  new_name+="_"; 

  for(string id : *full_id_ranges)
    new_name += id[0];

  auto ctrs_buff = make_shared<vector<pair<int,int>>>(*new_ctrs_pos);
  auto ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

  if (ctrs_buff_standard->size() !=0){
    new_name+="_"; 
    for(pair<int,int> ctr : *ctrs_buff_standard)
      new_name += to_string(ctr.first)+to_string(ctr.second);
  }

  return new_name;
}
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
  unc_id_ranges = make_shared<vector<string>>(0);
  unc_idxs = make_shared<vector<string>>(0);
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      unc_id_ranges->push_back(full_id_ranges->at(ii));
      unc_idxs->push_back(full_idxs->at(ii));
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
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrTensorPart<DType>::FullContract( shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                         shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                         shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::FullContract new" << endl; 
#endif 
/////////////////////////////////////////////////////////////////////////////////////////////////////
cout << endl <<  "CtrTensorPart<DType>::FullContract NEWVER : CTP name =  " << name << endl;

  while ( ctrs_todo->size() != 0 ){ 

    string CTP_in_name = get_next_name(ctrs_done);
    ctrs_done->push_back(ctrs_todo->back());
    string CTP_out_name = get_next_name(ctrs_done);
    ctrs_todo->pop_back();

    shared_ptr<CtrTensorPart<DType>> CTP_in;
    if ( Tmap->find(CTP_in_name) == Tmap->end()) {
      shared_ptr<vector<pair<int,int>>> ctrs_pos_in = make_shared<vector<pair<int,int>>>(*ctrs_todo);
      shared_ptr<vector<pair<int,int>>> new_ReIm_factors = make_shared<vector<pair<int,int>>>(1, make_pair(1,1));
      CTP_in = make_shared< CtrTensorPart<DType> >( full_idxs, full_id_ranges, ctrs_pos_in, new_ReIm_factors );
      Tmap->emplace(CTP_in->name,  CTP_in);
    } else {
      CTP_in = Tmap->at(CTP_in_name);
    }
    CTP_in->dependents.emplace(name);
  
    pair<int,int> ctrs_rel_pos_in = make_pair(CTP_in->unc_rel_pos->at(ctrs_done->back().first), CTP_in->unc_rel_pos->at(ctrs_done->back().second));   

    if ( ACompute_map->find(CTP_in_name) == ACompute_map->end()) {                                                       
      shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list_new = make_shared<vector<shared_ptr<CtrOp_base> >>(0);   
      CTP_in->FullContract(Tmap, ACompute_list_new, ACompute_map);                                                       
    }

    shared_ptr<CtrTensorPart<DType>> CTP_out;
    if ( Tmap->find(CTP_out_name) == Tmap->end()) {
      shared_ptr<vector<pair<int,int>>> ctrs_pos_in = make_shared<vector<pair<int,int>>>(*ctrs_done);
      shared_ptr<vector<pair<int,int>>> new_ReIm_factors = make_shared<vector<pair<int,int>>>(1, make_pair(1,1));
      CTP_out = make_shared< CtrTensorPart<DType> >( full_idxs, full_id_ranges, ctrs_done, new_ReIm_factors );
      Tmap->emplace(CTP_out_name,  CTP_out); 

    } else {
      CTP_out = Tmap->at(CTP_out_name);
    }
    CTP_out->dependencies.emplace(name);

    cout << "CTP Contract " << CTP_in_name << " over  (" << ctrs_done->back().first << ","<< ctrs_done->back().second << ") to get " << CTP_out_name ; cout.flush();
    cout << " added to " << name << "'s Acompute_list"<<  endl;
    ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done->back(), ctrs_rel_pos_in, "same_T new" ));
    shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list_out =  make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list);
    ACompute_map->emplace(CTP_out_name, ACompute_list_out);
    cout << "Acompute_list->size() =" << ACompute_list->size() << endl;
    dependencies.emplace(CTP_in_name);
  } 
  ACompute_map->emplace(myname(), ACompute_list);
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrMultiTensorPart<DType>::FullContract(shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                             shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list ,
                                             shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::FullContract" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << endl << "CtrMultiTensorPart<DType>::FullContract NEWVER :   CMTP name = " << myname() << endl;

   if (ctrs_pos->size() > 0 ) {
    
      cout << "CTP_vec->size() = " << CTP_vec->size() << endl;
      cout << "cross_ctrs_pos->size() = " <<  cross_ctrs_pos->size() << endl;
     if ( (CTP_vec->size() == 2) && ( cross_ctrs_pos->size() > 0 ) ) {
     
       shared_ptr<CtrTensorPart<DType>> new_CTP = Binary_Contract_diff_tensors(cross_ctrs_pos->back(), ctrs_pos->back(), Tmap,  ACompute_list, ACompute_map);
    
       if (Tmap->find(new_CTP->name) == Tmap->end())
         Tmap->emplace(new_CTP->name, new_CTP);
         
       if ( cross_ctrs_pos->size() >1 )   
         new_CTP->FullContract(Tmap, ACompute_list, ACompute_map);
       
     } else if ( cross_ctrs_pos->size() == 0 ) {
      
       for ( shared_ptr<CtrTensorPart<DType>> inner_CTP : *CTP_vec ) 
         inner_CTP->FullContract(Tmap, ACompute_list, ACompute_map);

     } else {
       cout << "USING MT BINARY CONTRACT DIFF TENSORS" << endl;
       auto new_CMTP = Binary_Contract_diff_tensors_MT(CTP_vec->at(cross_ctrs_pos->back().first.first)->myname(),CTP_vec->at(cross_ctrs_pos->back().second.first)->myname(),
                                                       ctrs_pos->back(), Tmap, ACompute_list, ACompute_map);
       new_CMTP->FullContract(Tmap, ACompute_list, ACompute_map);
     }
  }
  ACompute_map->emplace(name, ACompute_list);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<CtrTensorPart<DType>>
 CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors( pair<pair<int,int>, pair<int,int>> cross_ctr,// These two arguments should be 
                                                          pair<int,int> abs_ctr,                      // equivalent ! 
                                                          shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                                          shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                          shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors  NEWVER : " << name << endl; 

   shared_ptr<CtrTensorPart<DType>> T1;
   shared_ptr<CtrTensorPart<DType>> T2;

   int T1pos = cross_ctr.first.first;
   int T2pos = cross_ctr.second.first;
   int T1ctr = cross_ctr.first.second;
   int T2ctr = cross_ctr.second.second;

   CTP_vec->at(T1pos)->FullContract(Tmap, ACompute_list, ACompute_map);
   CTP_vec->at(T2pos)->FullContract(Tmap, ACompute_list, ACompute_map);
   int num_internal_ctrs = CTP_vec->at(T1pos)->ctrs_pos->size() + CTP_vec->at(T2pos)->ctrs_pos->size();

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

   auto full_id_ranges = make_shared<vector<string>>(T1->full_id_ranges->begin(), T1->full_id_ranges->end()) ;
   full_id_ranges->insert(full_id_ranges->end(), T2->full_id_ranges->begin(), T2->full_id_ranges->end()); 

   auto full_idxs = make_shared<vector<string>>(T1->full_idxs->begin(), T1->full_idxs->end()) ;
   full_idxs->insert(full_idxs->end(), T2->full_idxs->begin(), T2->full_idxs->end()) ;

   shared_ptr<vector<pair<int,int>>> ctrs_done = make_shared<vector<pair<int,int>>>(0);
   
   int T1shift =  Tsizes_cml->at(T1pos);
   int T2shift =  Tsizes_cml->at(T2pos);
   for (auto ctr : *T1->ctrs_pos)
     ctrs_done->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));
   for (auto ctr : *T2->ctrs_pos)
     ctrs_done->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));

   shared_ptr<vector<pair<int,int>>> ctrs_todo = make_shared<vector<pair<int,int>>>(0);
   for (int ii = 0 ;  ii != cross_ctrs_pos->size() ; ii++) 
     if ( (cross_ctrs_pos->at(ii).first.first + cross_ctrs_pos->at(ii).second.first) == (T1pos + T2pos) )
       ctrs_todo->push_back(make_pair(Tsizes_cml->at(cross_ctrs_pos->at(ii).first.first)  +  cross_ctrs_pos->at(ii).first.second,
                                     Tsizes_cml->at(cross_ctrs_pos->at(ii).second.first) +  cross_ctrs_pos->at(ii).second.second));

    
   auto full_ctrs = make_shared<vector<pair<int,int>>>(0);
   full_ctrs->insert(full_ctrs->end(), ctrs_done->begin(), ctrs_done->end());
   full_ctrs->insert(full_ctrs->end(), ctrs_todo->begin(), ctrs_todo->end());
 
   ctrs_done->push_back(ctrs_todo->back());
   ctrs_todo->pop_back();
  
   int T1_ctr_rel_pos = T1->unc_rel_pos->at(T1ctr);
   int T2_ctr_rel_pos = T2->unc_rel_pos->at(T2ctr);

   auto CTP_new = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, full_ctrs, make_shared<vector<pair<int,int>>>(1, abs_ctr )); 
   auto CTP_intermediate = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, ctrs_done, make_shared<vector<pair<int,int>>>( 1, abs_ctr )); 

   CTP_new->ctrs_todo = ctrs_todo;
   CTP_new->ctrs_done = ctrs_done;
   
   ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done),  abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
   auto new_CTData = make_shared<DType>();
    
   if (Tmap->find(CTP_new->name) == Tmap->end()){
     Tmap->emplace(CTP_new->name, CTP_new);
   }
   
   if (Tmap->find(CTP_intermediate->name) == Tmap->end()) {
     Tmap->emplace(CTP_intermediate->name, CTP_intermediate); 
   }
    
   cout << "BCDT contracting " << T1name << " and " << T2name << " over (" << abs_ctr.first << "," << abs_ctr.second << ") to get " << get_next_name(CTP_new->ctrs_done) << endl;
   return CTP_new;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<CtrMultiTensorPart<DType>>
CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT(string T1name, string T2name, pair<int,int> ctr,
                                                           shared_ptr< map<string, shared_ptr<CtrTensorPart<DType>>> > Tmap,
                                                           shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                           shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT" << endl; 

   auto T1T2_ctrd =  Binary_Contract_diff_tensors(cross_ctrs_pos->back(), ctr, Tmap, ACompute_list, ACompute_map );
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class CtrTensorPart<double>;
template class CtrMultiTensorPart<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
    
#endif
