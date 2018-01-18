#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/ctrtensop.h>
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/gamma_generator.h>

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

  shared_ptr<vector<pair<int,int>>> ctrs_buff = make_shared<vector<pair<int,int>>>(*ctrs_pos);
  shared_ptr<vector<pair<int,int>>> ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

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

  shared_ptr<vector<pair<int,int>>> ctrs_buff = make_shared<vector<pair<int,int>>>(*new_ctrs_pos);
  shared_ptr<vector<pair<int,int>>> ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

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
  for (int ii =0; ii != ctrs_pos->size() ; ii++){
    get_unc[ctrs_pos->at(ii).first] = false;
    get_unc[ctrs_pos->at(ii).second] = false;
  }

  bool survive_indep = true;
  int num_unc_ids =  get_unc.size() - ctrs_pos->size()*2;

  unc_id_ranges = make_shared<vector<string>>( num_unc_ids );
  unc_idxs = make_shared<vector<string>>( num_unc_ids );
  unc_pos = make_shared<vector<int>>( num_unc_ids );

  int jj = 0;
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      unc_id_ranges->at(jj) = full_id_ranges->at(ii);
      unc_idxs->at(jj)      = full_idxs->at(ii);
      unc_pos->at(jj)       = ii;
      jj++;
    }
  } 

  unc_rel_pos = make_shared<map<int,int>>();
  for( int ii =0 ; ii != unc_pos->size(); ii++) 
    unc_rel_pos->emplace(unc_pos->at(ii), ii);
 
  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a hack and should be done in a better way, maybe incorporate into constructor
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
pair<int,int> CtrTensorPart<DType>::get_pre_contract_ctr_rel_pos(pair<int,int>& ctr_pos ) { 
//////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "CtrTensorPart::get_pre_contract_ctr_rel_pos" << endl;

  vector<bool> get_unc(full_idxs->size(), true);
  for (int ii = 0 ; ii != ctrs_done->size()-1 ; ii++ )  {
    get_unc[ctrs_done->at(ii).first] = false;
    get_unc[ctrs_done->at(ii).second] = false;
  }
 
  int jj = 0;
  map<int,int> unc_rel_pos;
  for ( int ii = 0 ; ii != get_unc.size(); ii++ ) 
    if (get_unc[ii]){
      unc_rel_pos[ii] = jj;
      jj++;
    }
 
  return make_pair( unc_rel_pos[ctr_pos.first], unc_rel_pos[ctr_pos.second] ); 

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrTensorPart<DType>::FullContract( shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                         shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                         shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::FullContract" << endl; 
#endif 
/////////////////////////////////////////////////////////////////////////////////////////////////////
cout << endl <<  "CtrTensorPart<DType>::FullContract : CTP name =  " << name << endl;

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
 
    pair<int,int> ctrs_rel_pos_in = get_pre_contract_ctr_rel_pos( ctrs_done->back() ) ;

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

    if ( full_idxs->at(ctrs_done->back().first)[0] == 'X' || full_idxs->at(ctrs_done->back().second)[0] == 'X' ) { 
      if ( full_idxs->at(ctrs_done->back().first)[0] != 'X' || full_idxs->at(ctrs_done->back().second)[0] != 'X' ) { // TODO replace with CtrOp_single_id
        ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name << "'s Acompute_list"<<  endl;
      } else { //TODO replace with CtrOp_exc_ids
        ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name << "'s Acompute_list"<<  endl;
      }
    } else  {  
      cout << "CTP Contract " << CTP_in_name << " over  (" << ctrs_done->back().first << ","<< ctrs_done->back().second << ") to get " << CTP_out_name ; cout.flush();
      ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name << "'s Acompute_list"<<  endl;
    }

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
    
     cout << "CTP_vec->size() = " << CTP_vec->size() <<  "     cross_ctrs_pos->size() = " <<  cross_ctrs_pos->size() << endl;
     if ( (CTP_vec->size() == 2) && ( cross_ctrs_pos->size() > 0 ) ) {
     
       shared_ptr<CtrTensorPart<DType>> new_CTP = Binary_Contract_diff_tensors(cross_ctrs_pos->back(), ctrs_pos->back(), Tmap,  ACompute_list, ACompute_map);
    
       if (Tmap->find(new_CTP->name) == Tmap->end())
         Tmap->emplace(new_CTP->name, new_CTP);
         
       if ( cross_ctrs_pos->size() > 1 )   
         new_CTP->FullContract(Tmap, ACompute_list, ACompute_map);
       
     } else if ( cross_ctrs_pos->size() == 0 ) {
      
       for ( shared_ptr<CtrTensorPart<DType>> inner_CTP : *CTP_vec ) 
         inner_CTP->FullContract(Tmap, ACompute_list, ACompute_map);

     } else {
       cout << "USING MT BINARY CONTRACT DIFF TENSORS" << endl;
       shared_ptr<CtrMultiTensorPart> new_CMTP = Binary_Contract_diff_tensors_MT( CTP_vec->at(cross_ctrs_pos->back().first.first)->myname(),
                                                                                  CTP_vec->at(cross_ctrs_pos->back().second.first)->myname(),
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
   cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors : " << name << endl; 

   int T1pos, T2pos, T1ctr, T2ctr;

   //Swapping tensors round to maintain consistent ordering
   if (cross_ctr.first.first < cross_ctr.second.first) {
     T1pos = cross_ctr.first.first;
     T2pos = cross_ctr.second.first;
     T1ctr = cross_ctr.first.second;
     T2ctr = cross_ctr.second.second;
   } else { 
     T1pos = cross_ctr.second.first;
     T2pos = cross_ctr.first.first;
     T1ctr = cross_ctr.second.second;
     T2ctr = cross_ctr.first.second;
   }
   shared_ptr<CtrTensorPart<DType>> T1 = CTP_vec->at(T1pos);
   shared_ptr<CtrTensorPart<DType>> T2 = CTP_vec->at(T2pos);

   string T1name = T1->name; 
   string T2name = T2->name; 

   T1->FullContract(Tmap, ACompute_list, ACompute_map);
   T2->FullContract(Tmap, ACompute_list, ACompute_map);

   shared_ptr<vector<string>> full_id_ranges = make_shared<vector<string>>(T1->full_id_ranges->begin(), T1->full_id_ranges->end()) ;
   full_id_ranges->insert(full_id_ranges->end(), T2->full_id_ranges->begin(), T2->full_id_ranges->end()); 

   shared_ptr<vector<string>> full_idxs = make_shared<vector<string>>(T1->full_idxs->begin(), T1->full_idxs->end()) ;
   full_idxs->insert(full_idxs->end(), T2->full_idxs->begin(), T2->full_idxs->end()) ;

   shared_ptr<vector<pair<int,int>>> ctrs_done = make_shared<vector<pair<int,int>>>(0);
   
   int T1shift =  Tsizes_cml->at(T1pos);
   int T2shift =  Tsizes_cml->at(T2pos);
   for ( pair<int,int> ctr : *T1->ctrs_pos)
     ctrs_done->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));
   for ( pair<int,int> ctr : *T2->ctrs_pos)
     ctrs_done->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));

   shared_ptr<vector<pair<int,int>>> ctrs_todo = make_shared<vector<pair<int,int>>>(0);
   for (int ii = 0 ;  ii != cross_ctrs_pos->size() ; ii++) 
     if ( (cross_ctrs_pos->at(ii).first.first + cross_ctrs_pos->at(ii).second.first) == (T1pos + T2pos) )
       ctrs_todo->push_back(make_pair(Tsizes_cml->at(cross_ctrs_pos->at(ii).first.first)  +  cross_ctrs_pos->at(ii).first.second,
                                      Tsizes_cml->at(cross_ctrs_pos->at(ii).second.first) +  cross_ctrs_pos->at(ii).second.second));

    
   shared_ptr<vector<pair<int,int>>> full_ctrs = make_shared<vector<pair<int,int>>>( ctrs_done->begin(), ctrs_done->end() );
   full_ctrs->insert(full_ctrs->end(), ctrs_todo->begin(), ctrs_todo->end());
 
   ctrs_done->push_back(ctrs_todo->back());
   ctrs_todo->pop_back();
  
   int T1_ctr_rel_pos = T1->unc_rel_pos->at(T1ctr);
   int T2_ctr_rel_pos = T2->unc_rel_pos->at(T2ctr);

   shared_ptr<CtrTensorPart<DType>> CTP_new = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, full_ctrs, make_shared<vector<pair<int,int>>>(1, abs_ctr ) ); 
   CTP_new->ctrs_todo = ctrs_todo;
   CTP_new->ctrs_done = ctrs_done;
   Tmap->emplace(CTP_new->name, CTP_new);
    

    if ( full_idxs->at(abs_ctr.first)[0] == 'X' || full_idxs->at(abs_ctr.second)[0] == 'X' ) {
      if ( full_idxs->at(abs_ctr.first)[0] != 'X' || full_idxs->at(abs_ctr.second)[0] != 'X' ) { // TODO replace with CtrOp_single_id
        ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done),
                                                            abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
      } else {  // TODO replace with CtrOp_exc_ids 
        ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done),
                                                            abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
      }
    } else  {  
       ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done),
                                                           abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
    }

   ACompute_map->emplace( CTP_new->get_next_name(CTP_new->ctrs_done), make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list));
   
   //Add intermediate compute list into map seperately
   shared_ptr<CtrTensorPart<DType>> CTP_intermediate = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, ctrs_done, make_shared<vector<pair<int,int>>>( 1, abs_ctr )); 
   Tmap->emplace(CTP_intermediate->name, CTP_intermediate); 
    
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

   shared_ptr<CtrTensorPart<DType>> T1T2_ctrd =  Binary_Contract_diff_tensors(cross_ctrs_pos->back(), ctr, Tmap, ACompute_list, ACompute_map );
   Tmap->emplace(T1T2_ctrd->name, T1T2_ctrd);

   shared_ptr<vector<shared_ptr<CtrTensorPart<DType>>>> new_CTP_vec = make_shared<vector<shared_ptr<CtrTensorPart<DType>>>>(0);
   for ( shared_ptr<CtrTensorPart<DType>> CTP : *CTP_vec ) {  
     if (CTP->name == T1name || CTP->name == T2name)
       continue;
     new_CTP_vec->push_back(CTP);
   }
   
   CTP_vec->push_back(T1T2_ctrd);
   shared_ptr<vector<pair<pair<int,int>, pair<int,int>>>> new_cross_ctrs_pos = make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(cross_ctrs_pos->begin(), cross_ctrs_pos->end()-1);
   shared_ptr<CtrMultiTensorPart<DType>> new_CMTP = make_shared< CtrMultiTensorPart<DType> >(new_CTP_vec, new_cross_ctrs_pos ); 

   return new_CMTP;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class CtrTensorPart<double>;
template class CtrMultiTensorPart<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
    
#endif
