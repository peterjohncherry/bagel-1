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
void CtrMultiTensorPart<DType>::get_name_orig(){
/////////////////////////////////////////////////////////////////////////////
  name = "";
  for(string id : *full_idxs)
    name += id;
  name+="_"; 

  for(string id : *full_id_ranges)
    name += id[0];

  if (all_ctrs_pos->size() !=0){
    name+="_"; 
    for(pair<int,int> ctr : *all_ctrs_pos)
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
void CtrMultiTensorPart<DType>::get_name_readable(){
/////////////////////////////////////////////////////////////////////////////
  name = "";
  for(auto idxs : *full_idxs)
    name += idxs;// + '-'; 

  name+='_'; 
  for(string id : *full_id_ranges)
    name += id[0]; //+ '-';

  for (auto CTP : *CTP_vec){
    if (CTP->ctrs_pos->size() !=0){
      name+='_'; 
      for(pair<int,int> ctr : *CTP->ctrs_pos)
         name +="("+ CTP->full_idxs->at(ctr.first)+"."
                     + CTP->full_idxs->at(ctr.second)+")";
    }
  }

  if (cross_ctrs_pos->size() !=0){
    name+='_'; 
    for(pair<pair<int,int>,pair<int,int>> ctr : *cross_ctrs_pos){
      name +="["+ CTP_vec->at(ctr.first.first)->full_idxs->at(ctr.first.second)+"."
                  + CTP_vec->at(ctr.second.first)->full_idxs->at(ctr.second.second)+"]";
    }
  }
  return ;
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
//  vprint(*unc_pos, "unc_pos"); 
//  vprint(*idxs, "idxs"); 
//  vprint(*id_ranges, "id_ranges"); 
 
  return; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void CtrTensorPart<DType>::FullContract(shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                        shared_ptr<vector<tuple<string,string,pair<int,int>, string> >> ACompute_list ){
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
      
      new_CTP->CTdata = make_shared<DType>(); //CHANGE TO GET UNC TENSOR DATA
      Tmap->emplace(new_name, new_CTP);
    }
   

    if(ctrs_pos->size() > 1) { 
      auto new_ctrs_pos = make_shared<vector<pair<int,int>>>(ctrs_pos->begin(), ctrs_pos->end()-1);
      string new_name = get_next_name(new_ctrs_pos);
      auto new_ctp_loc = Tmap->find(new_name);
      
      // if find one less, contract it.
      if ( new_ctp_loc != Tmap->end() ){ 

        CTdata = new_ctp_loc->second->Binary_Contract_same_tensor(ctrs_pos->back(), ACompute_list);
        ACompute_list->push_back(tie(new_name, new_name , ctrs_pos->back(), name));
        Tmap->emplace(name , make_shared<CtrTensorPart>(*this));

      } else { 

        auto new_ReIm_factors = make_shared<vector<pair<int,int>>>(ReIm_factors->begin(), ReIm_factors->end()-1); 
        auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, new_ctrs_pos, new_ReIm_factors ); 
        new_CTP->FullContract(Tmap, ACompute_list);
        CTdata = new_CTP->Binary_Contract_same_tensor(ctrs_pos->back(), ACompute_list);
        Tmap->emplace(name , make_shared<CtrTensorPart>(*this));
        ACompute_list->push_back(tie(new_name, new_name , ctrs_pos->back(), name));

      }

    } else if(ctrs_pos->size() == 1 ){ 

      auto unc_ctrs_pos = make_shared<vector<pair<int,int>>>(1, make_pair(0,0) );
      auto unc_ReIm_factors = make_shared<vector<pair<int,int>>>(0); 
      auto unc_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, unc_ctrs_pos, unc_ReIm_factors ); 

      // Need to create a list of uncontracted tensors
      unc_CTP->CTdata = make_shared<DType>(); 
      Tmap->emplace(unc_CTP->name, unc_CTP);

      required_Tblocks->push_back(unc_CTP->name);

      ACompute_list->push_back(tie(unc_CTP->name, unc_CTP->name, ctrs_pos->back(), name));
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
                                             shared_ptr<vector<tuple<string,string,pair<int,int>,string> >> ACompute_list ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::FullContract" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
        auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, all_ctrs_pos, ReIm_factors ); 
        
        new_CTP->CTdata = new_ctp_loc->second->Binary_Contract_same_tensor(all_ctrs_pos->back(), ACompute_list);
        ACompute_list->push_back(tie(new_name, new_name, all_ctrs_pos->back(), new_CTP->name));
      
      } else {  
        auto new_CMTP = make_shared< CtrMultiTensorPart<DType> >(CTP_vec, new_cross_ctrs_pos); 
        new_CMTP->FullContract(Tmap,ACompute_list);

        auto new_CTP = make_shared< CtrTensorPart<DType> >(new_CMTP->full_idxs, new_CMTP->full_id_ranges, all_ctrs_pos, ReIm_factors ); 
        new_CTP->FullContract(Tmap, ACompute_list);
      }
     
    } else if(cross_ctrs_pos->size() == 1 ){

      if ( CTP_vec->size() == 2) {

        for (auto CTP : *CTP_vec) {
          CTP->FullContract(Tmap,ACompute_list);
        }
        
        auto new_CTP = Binary_Contract_diff_tensors(CTP_vec->at(cross_ctrs_pos->back().first.first)->myname(),
                                                    CTP_vec->at(cross_ctrs_pos->back().second.first)->myname(),
                                                    all_ctrs_pos->back(), Tmap , ACompute_list);
        Tmap->emplace(new_CTP->name, new_CTP);

      } else {
        auto new_CMTP = Binary_Contract_diff_tensors_MT(CTP_vec->at(cross_ctrs_pos->back().first.first)->myname(),
                                                        CTP_vec->at(cross_ctrs_pos->back().second.first)->myname(),
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
                                                           shared_ptr<vector<tuple<string,string,pair<int,int>,string > >> ACompute_list){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors_MT" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   auto T1T2_ctrd =  Binary_Contract_diff_tensors(T1name, T2name, ctr_todo, Tmap, ACompute_list);
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
 CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors(string T1name, string T2name, pair<int,int> ctr_todo,
                                                         shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                                         shared_ptr<vector<tuple<string,string,pair<int,int>,string> >> ACompute_list){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrMultiTensorPart
cout << "CtrMultiTensorPart<DType>::Binary_Contract_diff_tensors" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   auto T1 = Tmap->at(T1name);
   auto T2 = Tmap->at(T2name);
    
   int T1pos = 0;  
   for (   ; T1pos != CTP_vec->size() ; T1pos++) { if (T1name == CTP_vec->at(T1pos)->name) { break;}}
   
   int T2pos = 0;  
   for (   ; T2pos != CTP_vec->size() ; T2pos++) { if (T2name == CTP_vec->at(T2pos)->name) { break;}}
 
   auto full_id_ranges = make_shared<vector<string>>(0);
   full_id_ranges->insert(full_id_ranges->begin(), T1->full_id_ranges->begin(), T1->full_id_ranges->end()) ;
   full_id_ranges->insert(full_id_ranges->begin(), T2->full_id_ranges->begin(), T2->full_id_ranges->end()); 

   auto full_idxs = make_shared<vector<string>>(0);
   full_idxs->insert(full_idxs->begin(), T1->full_idxs->begin(), T1->full_idxs->end()) ;
   full_idxs->insert(full_idxs->begin(), T2->full_idxs->begin(), T2->full_idxs->end()) ;

   auto full_ctrs = make_shared<vector<pair<int,int>>>(0);
   int T1shift =  Tsizes_cml->at(T1pos);
   int T2shift =  Tsizes_cml->at(T2pos);

   for (auto ctr : *T1->ctrs_pos)
      full_ctrs->push_back( make_pair(ctr.first+T1shift, ctr.second+T1shift));

   for (auto ctr : *T2->ctrs_pos)
      full_ctrs->push_back( make_pair(ctr.first+T2shift, ctr.second+T2shift));

   full_ctrs->push_back(ctr_todo);
  
   auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, full_ctrs, make_shared<vector<pair<int,int>>>(1, make_pair(1,1) )); 

   auto new_CTData = make_shared<DType>();
   ACompute_list->push_back(tie(T1name, T2name, ctr_todo, new_CTP->name));

   new_CTP->CTdata= new_CTData;
   new_CTP->contracted = true;

   //cout << "BCDT contracting " << T1name << " and " << T2name << " over (" << ctr_todo.first << "," << ctr_todo.second << ") to get " << new_CTP->name << endl;
   return new_CTP;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<CtrTensorPart<DType>> CtrTensorPart<DType>::Binary_Contract_same_tensor( string T1name , pair<int,int> ctr_todo ,
                                                                                    shared_ptr<map<string,shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                                                                    shared_ptr<vector<tuple<string,string,pair<int,int>,string> >> ACompute_list ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::Binary_Contract_same_tensor (returns CTP)" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "BCST contracting  "; cout.flush(); cout << name << " over (" << ctr_todo.first << "," << ctr_todo.second << ") to get " ; cout.flush();
  auto T1 = Tmap->at(T1name);

  auto new_ctrs_pos = T1->ctrs_pos;
  auto new_ReIm_factors = T1->ReIm_factors;
  new_ctrs_pos->push_back(ctr_todo);
  new_ReIm_factors->push_back(make_pair(1,1));;

  auto new_CTP = make_shared< CtrTensorPart<DType> >(full_idxs, full_id_ranges, new_ctrs_pos, new_ReIm_factors); 

  CTdata = make_shared<DType>();

  return new_CTP;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<DType> CtrTensorPart<DType>::Binary_Contract_same_tensor(pair<int,int> ctr_todo, std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>,string> >> ACompute_list ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_CtrTensorPart
cout << "CtrTensorPart<DType>::Binary_Contract_same_tensor (returns data)" << endl; 
#endif 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "BCST contracting  "; cout.flush(); cout << name <<  " over (" << ctr_todo.first << "," << ctr_todo.second << ") to get " ; cout.flush();

//  string c1_idx = full_idxs->at(ctr_todo.first);  cout << "c1_idx = " <<c1_idx << endl;
//  string c2_idx = full_idxs->at(ctr_todo.second); cout << "c2_idx = " <<c2_idx << endl;
//  while (TP::idxs->at(relctr.first++)  != c1_idx) { cout << "TP::idxs->at("<<relctr.first-1 << ") = " << TP::idxs->at(relctr.first-1) << endl;   }
//  while (TP::idxs->at(relctr.second++) != c2_idx) { cout << "TP::idxs->at("<<relctr.second-1 << ") = " << TP::idxs->at(relctr.second-1) << endl;  }
//  BC_same_T(TP::id_ranges, relctr, inTdata);
  auto bob = make_shared<DType>();
  CTdata = make_shared<DType> ();
  return bob;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<vector<int>> CtrTensorPart<DType>::unc_id_ordering_with_ctr_at_front(int ctr_pos){
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "unc_id_ordering_with_ctr_at_front" << endl;
  cout << "unc_pos = [" ;for (int elem : *unc_pos ) {cout << elem << " " ; } cout << "]     ctr_pos = " << ctr_pos << endl;

  vector<int> new_order(unc_pos->size(),0);
  for (int pos : *unc_pos){
    cout << "pos = "<< pos << "   ctr_pos = " << ctr_pos <<  "     rel_pos = " << new_order.front() << endl;
    if ( pos  == ctr_pos ) {
      cout << "found  ctr_pos"<< endl;
      break;
    }
    new_order.front() = new_order.front()+1;
  }

  cout << "new_orderA = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
  int ii = 1;
  while ( ii < new_order.front()+1){
    new_order[ii] = ii-1;
    ii++;
  }

  cout << "new_orderB = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
  while ( ii < new_order.size()  ){
    new_order[ii] = ii;
    ii++;
  }
  
  cout << "new_orderC = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
  return make_shared<vector<int>>(new_order);
   
}     
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<vector<int>> CtrTensorPart<DType>::unc_id_ordering_with_ctr_at_back(int ctr_pos){
///////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "unc_id_ordering_with_ctr_at_back" << endl;
  cout << "unc_pos = [" ; for (int elem : *unc_pos ) {cout << elem << " " ; } cout << "]     ctr_pos = " << ctr_pos << endl;
  vector<int> new_order(unc_pos->size(),0);
  for (int pos : *unc_pos){
    cout << "pos = "<< pos << "   ctr_pos = " << ctr_pos <<  "     rel_pos = " << new_order.front() << endl;
    if ( pos  == ctr_pos ){
      cout << "found  ctr_pos"<< endl;
      break;
    }
    new_order.back() = new_order.back()+1;
  }
 
  cout << "new_orderA = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
  int ii = 0;
  while ( ii < new_order.back()){
    new_order[ii] = ii;
    cout <<" first loop ii =  " << ii << endl;
    ii++;
  }

  cout << "new_orderB = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
  while ( ii < new_order.size()-1 ){
    cout <<" second loop ii =  " << ii << endl;
    new_order[ii] = ii+1;
    ii++;
  }
  
  cout << "new_orderC = [ " ; cout.flush(); for ( int ii : new_order ) { cout << ii << " " ; cout.flush(); }; cout << "]" << endl;
 
  return make_shared<vector<int>>(new_order);
}     

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensorPart<double>;
template class CtrTensorPart<double>;
template class CtrMultiTensorPart<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
    
#endif
