#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator.h>
#include <cassert>
using namespace std;
  
///////////////////////////////////////////////////////////////////////////
void CtrTensorPart_Base::get_name(){
///////////////////////////////////////////////////////////////////////////
  name_ = "";
  for(string id : *full_idxs_)
    name_ += id;
  name_ +="_"; 

  for(string id : *full_id_ranges_)
    name_ += id[0];

  shared_ptr<vector<pair<int,int>>> ctrs_buff = make_shared<vector<pair<int,int>>>(*ctrs_pos_);
  shared_ptr<vector<pair<int,int>>> ctrs_buff_standard = GammaGenerator::Standardize_delta_ordering_generic(ctrs_buff ) ;

  if (ctrs_buff_standard->size() !=0){
    name_ +="_"; 
    for(pair<int,int> ctr : *ctrs_buff_standard)
      name_ += to_string(ctr.first)+to_string(ctr.second);
  }
  return;
};

/////////////////////////////////////////////////////////////////////////////
string CtrTensorPart_Base::get_next_name(shared_ptr<vector<pair<int,int>>> new_ctrs_pos ){
/////////////////////////////////////////////////////////////////////////////
  string new_name = "";
  for(string id : *full_idxs_)
    new_name += id;
  new_name+="_"; 

  for(string id : *full_id_ranges_)
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
void CtrTensorPart_Base::get_ctp_idxs_ranges(){
//////////////////////////////////////////////////////////////////////////////
//cout << "CtrTensorPart_Base::get_ctp_idxs_ranges() " << name_  << endl;

  int counter = 0;
  vector<bool> get_unc(full_idxs_->size(), true);

  for (int ii =0; ii != ctrs_pos_->size() ; ii++){
    get_unc[ctrs_pos_->at(ii).first] = false;
    get_unc[ctrs_pos_->at(ii).second] = false;
  }

  int num_unc_ids =  get_unc.size() - ctrs_pos_->size()*2;

  unc_id_ranges_ = make_shared<vector<string>>( num_unc_ids );
  unc_idxs_ = make_shared<vector<string>>( num_unc_ids );
  unc_pos_ = make_shared<vector<int>>( num_unc_ids );

  int jj = 0;
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      unc_id_ranges_->at(jj) = full_id_ranges_->at(ii); 
      unc_idxs_->at(jj)      = full_idxs_->at(ii);      
      unc_pos_->at(jj)       = ii;                      
      jj++;
    }
  } 

  unc_rel_pos_ = make_shared<map<int,int>>();
  for( int ii = 0 ; ii != unc_pos_->size(); ii++)
    unc_rel_pos_->emplace(unc_pos_->at(ii), ii);

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a hack and should be done in a better way, maybe incorporate into constructor
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
pair<int,int> CtrTensorPart<DataType>::get_pre_contract_ctr_rel_pos(pair<int,int>& ctr_pos ) { 
//////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "CtrTensorPart::get_pre_contract_ctr_rel_pos" << endl;

  vector<bool> get_unc(full_idxs_->size(), true);
  for ( int ii = 0 ; ii != ctrs_done_->size()-1 ; ii++ )  {
    get_unc[ctrs_done_->at(ii).first] = false;
    get_unc[ctrs_done_->at(ii).second] = false;
  }
 
  int jj = 0;
  map<int,int> unc_rel_pos_;
  for ( int ii = 0 ; ii != get_unc.size(); ii++ ) 
    if (get_unc[ii]){
      unc_rel_pos_[ii] = jj;
      jj++;
    }
 
  return make_pair( unc_rel_pos_[ctr_pos.first], unc_rel_pos_[ctr_pos.second] ); 

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void CtrTensorPart<DataType>::FullContract( shared_ptr<map<string,shared_ptr<CtrTensorPart_Base> >> Tmap,
                                            shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                            shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DataType>::FullContract" << endl; 
#endif 
/////////////////////////////////////////////////////////////////////////////////////////////////////
cout << endl <<  "CtrTensorPart<DataType>::FullContract : CTP name =  " << name_ << endl;
int counter = 0 ;
  while ( ctrs_todo_->size() != 0 ){ 

    string CTP_in_name = get_next_name(ctrs_done_);
    ctrs_done_->push_back(ctrs_todo_->back());
    string CTP_out_name = get_next_name(ctrs_done_);
    ctrs_todo_->pop_back();

    shared_ptr<CtrTensorPart_Base> CTP_in;
    if ( Tmap->find(CTP_in_name) == Tmap->end()) {
      shared_ptr<vector<pair<int,int>>> new_ctrs_pos = make_shared<vector<pair<int,int>>>(*ctrs_todo_);
      shared_ptr<vector<pair<int,int>>> new_ReIm_factors = make_shared<vector<pair<int,int>>>(1, make_pair(1,1));

      CTP_in = make_shared< CtrTensorPart<DataType> >( full_idxs_, full_id_ranges_, new_ctrs_pos, new_ReIm_factors );
      Tmap->emplace(CTP_in->name(),  CTP_in);
    } else {
      CTP_in = Tmap->at(CTP_in_name);
      assert(CTP_in);
    }
    CTP_in->dependents_.emplace(CTP_out_name);
    CTP_in->dependents_.emplace(name_);
 
    pair<int,int> ctrs_rel_pos_in = get_pre_contract_ctr_rel_pos( ctrs_done_->back() ) ;

    if ( ACompute_map->find(CTP_in_name) == ACompute_map->end()) {                                                       
      shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list_new = make_shared<vector<shared_ptr<CtrOp_base> >>(0);   
      CTP_in->FullContract(Tmap, ACompute_list_new, ACompute_map);                                                       
    }

    shared_ptr<CtrTensorPart_Base> CTP_out;
    if ( Tmap->find(CTP_out_name) == Tmap->end()) {

      shared_ptr<vector<pair<int,int>>> new_ctrs_pos = make_shared<vector<pair<int,int>>>(*ctrs_done_);
      shared_ptr<vector<pair<int,int>>> new_ReIm_factors = make_shared<vector<pair<int,int>>>(1, make_pair(1,1));
      CTP_out = make_shared< CtrTensorPart<DataType> >( full_idxs_, full_id_ranges_, new_ctrs_pos, new_ReIm_factors );
      Tmap->emplace(CTP_out_name,  CTP_out); 

    } else {
      CTP_out = Tmap->at(CTP_out_name);
    }
    CTP_out->dependencies_.emplace(name_);

    // TODO This does nothing at the moment, but you want these ifs here when you do the proper treatment of excitation ops, so keep it!
    if ( full_idxs_->at(ctrs_done_->back().first)[0] == 'X' || full_idxs_->at(ctrs_done_->back().second)[0] == 'X' ) { 
      if ( full_idxs_->at(ctrs_done_->back().first)[0] != 'X' || full_idxs_->at(ctrs_done_->back().second)[0] != 'X' ) { // TODO replace with CtrOp_single_id
        ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done_->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name_ << "'s Acompute_list"<<  endl;
      } else { //TODO replace with CtrOp_exc_ids
        ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done_->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name_ << "'s Acompute_list"<<  endl;
      }
    } else  {  
      ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done_->back(), ctrs_rel_pos_in, "same_T new" )); cout << " added to " << name_ << "'s Acompute_list"<<  endl;
    }
    cout << "CTP Contract " << CTP_in_name << " over  (" << ctrs_done_->back().first << ","<< ctrs_done_->back().second << ") to get " << CTP_out_name ; cout.flush();

    shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list_out =  make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list);
    ACompute_map->emplace(CTP_out_name, ACompute_list_out);

    cout << "Acompute_list->size() =" << ACompute_list->size() << endl;
    dependencies_.emplace(CTP_in_name);
  } 
  ACompute_map->emplace(name_, ACompute_list);
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void CtrMultiTensorPart<DataType>::FullContract( shared_ptr<map<string,shared_ptr<CtrTensorPart_Base> >> Tmap,
                                                 shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list ,
                                                 shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DataType>::FullContract" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << endl << "CtrMultiTensorPart<DataType>::FullContract NEWVER :   CMTP name = " << name_ << endl;

  if ( get_compute_list_from_reordered_tens_ ) {  
              
    Tmap->at(reordered_tens_name_)->FullContract( Tmap, ACompute_list, ACompute_map) ; 
    ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >> (*(ACompute_map->at( reordered_tens_name_ ))); 

    ACompute_list->push_back( make_shared<CtrOp_reorder> ( reordered_tens_name_, name_, reordering_, "reordering" ));

  } else if (ctrs_pos_->size() > 0 ) {
    
    cout << "CTP_vec->size() = " << CTP_vec_->size() <<  "     cross_ctrs_pos_->size() = " <<  cross_ctrs_pos_->size() << endl;
    if ( (CTP_vec_->size() == 2) && ( cross_ctrs_pos_->size() > 0 ) ) {
    
      shared_ptr<CtrTensorPart<DataType>> new_CTP = Binary_Contract_diff_tensors(cross_ctrs_pos_->back(), ctrs_pos_->back(), Tmap,  ACompute_list, ACompute_map);
    
      if (Tmap->find(new_CTP->name_) == Tmap->end())
        Tmap->emplace(new_CTP->name_, new_CTP);
        
      if ( cross_ctrs_pos_->size() > 1 )   
        new_CTP->FullContract(Tmap, ACompute_list, ACompute_map);
      
    } else if ( cross_ctrs_pos_->size() == 0 ) {
     
      for ( auto& inner_CTP : *CTP_vec_ ) 
        inner_CTP->FullContract(Tmap, ACompute_list, ACompute_map);

    } else {
      throw logic_error("USING MT BINARY CONTRACT DIFF TENSORS MT!! This shouldn't happen anymore ... Aborting !! ");
      shared_ptr<CtrMultiTensorPart> new_CMTP = Binary_Contract_diff_tensors_MT( CTP_vec_->at(cross_ctrs_pos_->back().first.first)->name_,
                                                                                 CTP_vec_->at(cross_ctrs_pos_->back().second.first)->name_,
                                                                                 ctrs_pos_->back(), Tmap, ACompute_list, ACompute_map);
      new_CMTP->FullContract(Tmap, ACompute_list, ACompute_map);
    }
  }
  ACompute_map->emplace(name_, ACompute_list);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<CtrTensorPart<DataType>>
 CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors( pair<pair<int,int>, pair<int,int>> cross_ctr,// These two arguments should be 
                                                          pair<int,int> abs_ctr,                      // equivalent ! 
                                                          shared_ptr<map<string,shared_ptr<CtrTensorPart_Base > >> Tmap,
                                                          shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                          shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors : " << name_ << endl; 

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
   shared_ptr<CtrTensorPart_Base> T1 = CTP_vec_->at(T1pos);
   shared_ptr<CtrTensorPart_Base> T2 = CTP_vec_->at(T2pos);

   string T1name = T1->name_; 
   string T2name = T2->name_; 

   T1->FullContract(Tmap, ACompute_list, ACompute_map);
   T2->FullContract(Tmap, ACompute_list, ACompute_map);

   shared_ptr<vector<string>> full_id_ranges_ = make_shared<vector<string>>(T1->full_id_ranges_->begin(), T1->full_id_ranges_->end()) ;
   full_id_ranges_->insert(full_id_ranges_->end(), T2->full_id_ranges_->begin(), T2->full_id_ranges_->end()); 

   shared_ptr<vector<string>> full_idxs_ = make_shared<vector<string>>(T1->full_idxs_->begin(), T1->full_idxs_->end()) ;
   full_idxs_->insert(full_idxs_->end(), T2->full_idxs_->begin(), T2->full_idxs_->end()) ;

   shared_ptr<vector<pair<int,int>>> ctrs_done_ = make_shared<vector<pair<int,int>>>(0);
   
   int T1shift =  Tsizes_cml->at(T1pos);
   int T2shift =  Tsizes_cml->at(T2pos);
   for ( pair<int,int> ctr : *T1->ctrs_pos_)
     ctrs_done_->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));
   for ( pair<int,int> ctr : *T2->ctrs_pos_)
     ctrs_done_->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));

   shared_ptr<vector<pair<int,int>>> ctrs_todo_ = make_shared<vector<pair<int,int>>>(0);
   for (int ii = 0 ;  ii != cross_ctrs_pos_->size() ; ii++) 
     if ( (cross_ctrs_pos_->at(ii).first.first + cross_ctrs_pos_->at(ii).second.first) == (T1pos + T2pos) )
       ctrs_todo_->push_back(make_pair( Tsizes_cml->at( cross_ctrs_pos_->at(ii).first.first )  +  cross_ctrs_pos_->at(ii).first.second,
                                        Tsizes_cml->at( cross_ctrs_pos_->at(ii).second.first ) +  cross_ctrs_pos_->at(ii).second.second ));
    
   shared_ptr<vector<pair<int,int>>> full_ctrs = make_shared<vector<pair<int,int>>>( ctrs_done_->begin(), ctrs_done_->end() );
   full_ctrs->insert(full_ctrs->end(), ctrs_todo_->begin(), ctrs_todo_->end());
 
   ctrs_done_->push_back(ctrs_todo_->back());
   ctrs_todo_->pop_back();
  
   int T1_ctr_rel_pos = T1->unc_rel_pos_->at(T1ctr);
   int T2_ctr_rel_pos = T2->unc_rel_pos_->at(T2ctr);

   shared_ptr<CtrTensorPart<DataType>> CTP_new = make_shared< CtrTensorPart<DataType> >(full_idxs_, full_id_ranges_, full_ctrs, make_shared<vector<pair<int,int>>>(1, abs_ctr ) ); 
   CTP_new->ctrs_todo_ = ctrs_todo_;
   CTP_new->ctrs_done_ = ctrs_done_;
   Tmap->emplace(CTP_new->name_, CTP_new);
    
    // TODO This does nothing at the moment, but you want these ifs here when you do the proper treatment of excitation ops, so keep it!
    if ( full_idxs_->at(abs_ctr.first)[0] == 'X' || full_idxs_->at(abs_ctr.second)[0] == 'X' ) {
      if ( full_idxs_->at(abs_ctr.first)[0] != 'X' || full_idxs_->at(abs_ctr.second)[0] != 'X' ) { // TODO replace with CtrOp_single_id
        ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done_),
                                                            abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
      } else {  // TODO replace with CtrOp_exc_ids 
        ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done_),
                                                            abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
      }
    } else  {  
       ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done_),
                                                           abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
    }

   ACompute_map->emplace( CTP_new->get_next_name(CTP_new->ctrs_done_), make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list));
   
   //Add intermediate compute list into map seperately
   shared_ptr<CtrTensorPart<DataType>> CTP_intermediate = make_shared< CtrTensorPart<DataType> >(full_idxs_, full_id_ranges_, ctrs_done_, make_shared<vector<pair<int,int>>>( 1, abs_ctr )); 
   Tmap->emplace(CTP_intermediate->name_, CTP_intermediate); 
    
   cout << "BCDT contracting " << T1name << " and " << T2name << " over (" << abs_ctr.first << "," << abs_ctr.second << ") to get " << get_next_name(CTP_new->ctrs_done_) << endl;
   return CTP_new;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<CtrMultiTensorPart<DataType>>
CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors_MT( string T1name, string T2name, pair<int,int> ctr,
                                                               shared_ptr< map<string, shared_ptr<CtrTensorPart_Base>> > Tmap,
                                                               shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                               shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors_MT" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors_MT" << endl; 
   throw logic_error( " Should not be using this rountine anymore !!! Kept only for checking; Aborting!! " ); 
   for ( pair<pair<int,int> , pair<int,int>>  cross_ctr : *cross_ctrs_pos_ ) {
     int T1_pos = cross_ctr.first.first;
     int T2_pos = cross_ctr.second.first;
     int T1_ctr_pos = cross_ctr.second.first;
     int T2_ctr_pos = cross_ctr.second.second;
     cout << " {[" << T1_pos << ":" << T1_ctr_pos << "].[ " << T2_pos << ":" << T2_ctr_pos << "]}" ; cout.flush();
   } cout << endl;

   shared_ptr<CtrTensorPart<DataType>> T1T2_ctrd =  Binary_Contract_diff_tensors(cross_ctrs_pos_->back(), ctr, Tmap, ACompute_list, ACompute_map );
   Tmap->emplace(T1T2_ctrd->name_, T1T2_ctrd);

   shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> new_CTP_vec = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>(0);
   for ( shared_ptr<CtrTensorPart_Base> CTP : *CTP_vec_ ) {  
     if (CTP->name_ == T1name || CTP->name_ == T2name)
       continue;
     new_CTP_vec->push_back(CTP);
   }
   
   CTP_vec_->push_back(T1T2_ctrd);
   shared_ptr<vector<pair<pair<int,int>, pair<int,int>>>> new_cross_ctrs_pos_ = make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(cross_ctrs_pos_->begin(), cross_ctrs_pos_->end()-1);
   shared_ptr<CtrMultiTensorPart<DataType>> new_CMTP = make_shared< CtrMultiTensorPart<DataType> >(new_CTP_vec, new_cross_ctrs_pos_ ); 

   return new_CMTP;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class CtrTensorPart<double>;
template class CtrTensorPart<complex<double>>;
template class CtrMultiTensorPart<double>;
template class CtrMultiTensorPart<complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
