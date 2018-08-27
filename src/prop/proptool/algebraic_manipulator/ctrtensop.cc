#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/proputils.h>
#include <cassert>
using namespace std;
/////////////////////////////////////////////////////////////////////////////////////////////
string CtrTensorPart_Base::get_next_name(shared_ptr<vector<pair<int,int>>> new_ctrs_pos ){
/////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CTR_TENSOR_PART
cout << "CtrTensorPart_Base::get_next_name" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////
  
  string new_name = op_info_->op_state_name_canonical();

  new_name+="_"; 

  for(string id : *full_id_ranges_)
    new_name += id[0];

  shared_ptr<vector<pair<int,int>>> ctrs_buff = make_shared<vector<pair<int,int>>>(*new_ctrs_pos);
  shared_ptr<vector<pair<int,int>>> ctrs_buff_standard = WickUtils::standardize_delta_ordering_generic(*ctrs_buff, *full_idxs_ ) ;

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
#ifdef __DEBUG_PROPTOOL_CTR_TENSOR_PART
cout << "CtrTensorPart_Base::get_ctp_idxs_ranges() " << name_  << endl;
#endif ///////////////////////////////////////////////////////////////////////

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
#ifdef __DEBUG_PROPTOOL_CTR_TENSOR_PART
cout << "CtrTensorPart::get_pre_contract_ctr_rel_pos" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void CtrTensorPart<DataType>::build_contraction_sequence( shared_ptr<map<string,shared_ptr<CtrTensorPart_Base> >> Tmap,
                                                          shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                          shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CTR_TENSOR_PART
cout << endl <<  "CtrTensorPart<DataType>::build_contraction_sequence : CTP name =  " << name_ << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( ctrs_pos_->size() == 0 ) {

    ACompute_list->push_back( make_shared<CtrOp_Get>( name_, "get" ) ); 

    return;

  } else {  
 
    while ( ctrs_todo_->size() != 0 ){ 
  
      string CTP_in_name = get_next_name(ctrs_done_);
      ctrs_done_->push_back(ctrs_todo_->back());
      string CTP_out_name = get_next_name(ctrs_done_);
      ctrs_todo_->pop_back();
    
      shared_ptr<CtrTensorPart_Base> CTP_in;
      if ( Tmap->find(CTP_in_name) == Tmap->end()) {
        shared_ptr<vector<pair<int,int>>> new_ctrs_pos = make_shared<vector<pair<int,int>>>(*ctrs_todo_);
  
        CTP_in = make_shared< CtrTensorPart<DataType> >( full_idxs_, full_id_ranges_, new_ctrs_pos, op_info_ );
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
        CTP_in->build_contraction_sequence(Tmap, ACompute_list_new, ACompute_map);                                                       
      }
  
      shared_ptr<CtrTensorPart_Base> CTP_out;
      if ( Tmap->find(CTP_out_name) == Tmap->end()) {
  
        shared_ptr<vector<pair<int,int>>> new_ctrs_pos = make_shared<vector<pair<int,int>>>(*ctrs_done_);
        CTP_out = make_shared< CtrTensorPart<DataType> >( full_idxs_, full_id_ranges_, new_ctrs_pos, op_info_ );
        Tmap->emplace(CTP_out_name,  CTP_out); 
  
      } else {
        CTP_out = Tmap->at(CTP_out_name);
      }
      CTP_out->dependencies_.emplace(name_);
  
      ACompute_list->push_back( make_shared<CtrOp_same_T> (CTP_in_name, CTP_out_name, ctrs_done_->back(), ctrs_rel_pos_in, "same_T new" ));

#ifdef __DEBUG_PROPTOOL_CTR_TENSOR_PART
      cout << " added to " << name_ << "'s Acompute_list"<<  endl;
      cout << "CTP Contract " << CTP_in_name << " over  (" << ctrs_done_->back().first << ","<< ctrs_done_->back().second << ") to get " << CTP_out_name ; cout.flush(); 
#endif
  
      shared_ptr<vector<shared_ptr<CtrOp_base>>> ACompute_list_out =  make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list);
      ACompute_map->emplace(CTP_out_name, ACompute_list_out);
  
      dependencies_.emplace(CTP_in_name);
    } 
    ACompute_map->emplace(name_, ACompute_list);
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void CtrMultiTensorPart<DataType>::build_contraction_sequence( shared_ptr<map<string,shared_ptr<CtrTensorPart_Base> >> Tmap,
                                                               shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list ,
                                                               shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CTR_MULTITENSOR_PART
cout << endl << "CtrMultiTensorPart<DataType>::build_contraction_sequence :   CMTP name = " << name_ << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( get_compute_list_from_reordered_tens_ ) {  // TODO Should not need this anymore now op_info class has been extended; all Acontrib should be in canonical order on entry.
              
    Tmap->at(reordered_tens_name_)->build_contraction_sequence( Tmap, ACompute_list, ACompute_map) ; 
    ACompute_list = make_shared<vector<shared_ptr<CtrOp_base> >> (*(ACompute_map->at( reordered_tens_name_ ))); 
    ACompute_list->push_back( make_shared<CtrOp_reorder> ( reordered_tens_name_, name_, reordering_, "reordering" ));

  } else if (ctrs_pos_->size() > 0 ) {
    
    if ( (CTP_vec_->size() == 2) && ( cross_ctrs_pos_->size() > 0 ) ) {

      shared_ptr<CtrTensorPart<DataType>> new_CTP = Binary_Contract_diff_tensors(cross_ctrs_pos_->back(), ctrs_pos_->back(), Tmap,  ACompute_list, ACompute_map);
      if (Tmap->find(new_CTP->name_) == Tmap->end())
        Tmap->emplace(new_CTP->name_, new_CTP);
      

      if ( cross_ctrs_pos_->size() > 1 )
        new_CTP->build_contraction_sequence(Tmap, ACompute_list, ACompute_map);
      
    } else if ( cross_ctrs_pos_->size() == 0 ) {
     
      shared_ptr<vector<string>> direct_product_tensor_list = make_shared<vector<string>>( CTP_vec_->size()); 
      vector<string>::iterator dptl_it = direct_product_tensor_list->begin();
      for ( vector<shared_ptr<CtrTensorPart_Base>>::iterator cv_it = CTP_vec_->begin(); cv_it != CTP_vec_->end() ;  cv_it++, dptl_it++ ){ 
        (*cv_it)->build_contraction_sequence(Tmap, ACompute_list, ACompute_map);
        *dptl_it = (*cv_it)->name();
      }
      
      ACompute_list->push_back ( make_shared<CtrOp_DirectProduct> ( direct_product_tensor_list, name_, "cartesian") );

    } else {
      throw logic_error("CtrMultiTensorPart<DataType>::build_contraction_sequence ; Should always meet one of the above conditions... Aborting !! ");
    }
  }
  ACompute_map->emplace(name_, ACompute_list);
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<CtrTensorPart<DataType>>
 CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors( pair<pair<int,int>, pair<int,int>> cross_ctr,// These two arguments should be 
                                                          pair<int,int> abs_ctr,                      // equivalent ! 
                                                          shared_ptr<map<string,shared_ptr<CtrTensorPart_Base > >> Tmap,
                                                          shared_ptr<vector<shared_ptr<CtrOp_base> >> ACompute_list,
                                                          shared_ptr<map<string, shared_ptr<vector<shared_ptr<CtrOp_base>> > >> ACompute_map ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CTR_MULTITENSOR_PART
cout << "CtrMultiTensorPart<DataType>::Binary_Contract_diff_tensors : " << name_ << endl; 
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

   T1->build_contraction_sequence(Tmap, ACompute_list, ACompute_map);
   T2->build_contraction_sequence(Tmap, ACompute_list, ACompute_map);

   shared_ptr<vector<string>> full_id_ranges_ = make_shared<vector<string>>(T1->full_id_ranges_->begin(), T1->full_id_ranges_->end()) ;
   full_id_ranges_->insert(full_id_ranges_->end(), T2->full_id_ranges_->begin(), T2->full_id_ranges_->end()); 

   shared_ptr<vector<string>> full_idxs_ = make_shared<vector<string>>(T1->full_idxs_->begin(), T1->full_idxs_->end()) ;
   full_idxs_->insert(full_idxs_->end(), T2->full_idxs_->begin(), T2->full_idxs_->end()) ;

   shared_ptr<vector<pair<int,int>>> ctrs_done_ = make_shared<vector<pair<int,int>>>(0);
   
   int T1shift =  Tsizes_cml->at(T1pos);
   int T2shift =  Tsizes_cml->at(T2pos);

   for ( pair<int,int> ctr : *T1->ctrs_pos_){ 
     ctrs_done_->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));
   }

   for ( pair<int,int> ctr : *T2->ctrs_pos_){
     ctrs_done_->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));
   }

   shared_ptr<vector<pair<int,int>>> ctrs_todo_ = make_shared<vector<pair<int,int>>>(0);
   for (int ii = 0 ;  ii != cross_ctrs_pos_->size() ; ii++) {
     if ( (cross_ctrs_pos_->at(ii).first.first + cross_ctrs_pos_->at(ii).second.first) == (T1pos + T2pos) )
       ctrs_todo_->push_back(make_pair( Tsizes_cml->at( cross_ctrs_pos_->at(ii).first.first )  +  cross_ctrs_pos_->at(ii).first.second,
                                        Tsizes_cml->at( cross_ctrs_pos_->at(ii).second.first ) +  cross_ctrs_pos_->at(ii).second.second ));
   }

   shared_ptr<vector<pair<int,int>>> full_ctrs = make_shared<vector<pair<int,int>>>( ctrs_done_->begin(), ctrs_done_->end() );
   full_ctrs->insert(full_ctrs->end(), ctrs_todo_->begin(), ctrs_todo_->end());
 
   ctrs_done_->push_back(ctrs_todo_->back());
   ctrs_todo_->pop_back();
  
   int T1_ctr_rel_pos = T1->unc_rel_pos_->at(T1ctr);
   int T2_ctr_rel_pos = T2->unc_rel_pos_->at(T2ctr);
   shared_ptr<CtrTensorPart<DataType>> CTP_new = make_shared< CtrTensorPart<DataType> >(full_idxs_, full_id_ranges_, full_ctrs, op_info_ ); 
   CTP_new->ctrs_todo_ = ctrs_todo_;
   CTP_new->ctrs_done_ = ctrs_done_;
   Tmap->emplace(CTP_new->name_, CTP_new);
   
   ACompute_list->push_back(make_shared<CtrOp_diff_T>( T1name, T2name, CTP_new->get_next_name(CTP_new->ctrs_done_),
                                                       abs_ctr.first, abs_ctr.second, T1_ctr_rel_pos, T2_ctr_rel_pos, "diff_T_prod"));
 
   ACompute_map->emplace( CTP_new->get_next_name(CTP_new->ctrs_done_), make_shared<vector<shared_ptr<CtrOp_base>>>(*ACompute_list));
   
   //Add intermediate compute list into map seperately
   shared_ptr<CtrTensorPart<DataType>> CTP_intermediate = make_shared< CtrTensorPart<DataType> >(full_idxs_, full_id_ranges_, ctrs_done_, op_info_ ); 
   Tmap->emplace(CTP_intermediate->name_, CTP_intermediate); 
    
   cout << "BCDT contracting " << T1name << " and " << T2name << " over (" << abs_ctr.first << "," << abs_ctr.second << ") to get " << get_next_name(CTP_new->ctrs_done_) << endl;
   return CTP_new;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class CtrTensorPart<double>;
template class CtrTensorPart<complex<double>>;
template class CtrMultiTensorPart<double>;
template class CtrMultiTensorPart<complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
