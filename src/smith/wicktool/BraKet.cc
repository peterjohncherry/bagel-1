#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/BraKet.h>
//#include "BraKet.h"

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;


      using pint_vec = std::vector<std::pair<int,int>>;
      using pstr_vec = std::vector<std::pair<std::string,std::string>>;
      
      using CombinedGammaMap = std::map<std::vector<std::string>, /*spins of gamma_ids*/ 
                                        std::tuple< std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* contractions from reordering to gamma */  
                                                    std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* spins of contractions from reordering to gamma */  
                                                    std::shared_ptr<std::vector<int>>  /* signs from reordering to gamma */ > >;

      using RelCombinedGammaMap = std::map< std::pair<std::vector<std::string>, std::pair<int,int>>,  /*spins of gamma_ids and original spinsector*/ 
                                            std::tuple< std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* contractions from reordering to gamma */  
                                                        std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* spins of contractions from reordering to gamma */  
                                                        std::shared_ptr<std::vector<int>>  /* signs from reordering to gamma */ > >;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
BraKet<DType>::BraKet(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Sub_Ops = make_shared<vector<shared_ptr<TensOp<DType>>>>(0);
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::add_Op(string op_name,
                    shared_ptr<vector<string>> op_idxs,
                    shared_ptr<vector<bool>> op_aops, 
                    shared_ptr<vector<vector<string>>> op_idx_ranges,
                    vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                    vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                    pair<double,double> factor, string Tsymmetry, bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(hconj){
    reverse(op_idxs->begin(), op_idxs->end());  
    reverse(op_aops->begin(), op_aops->end());  
    reverse(op_idx_ranges->begin(), op_idx_ranges->end());  
    factor.second = -factor.second;// change this 
  }

  //NOTE: change to use proper factor
  int tmpfac =1;
  shared_ptr<TensOp<DType>>  New_Op = make_shared<TensOp<DType>>(op_name, Symmetry_Funcs, Constraint_Funcs);

  New_Op->initialize(*op_idxs, *op_idx_ranges, *op_aops, tmpfac, Tsymmetry);
  New_Op->get_ctrs_tens_ranges();

  Sub_Ops->push_back(New_Op);

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_TotalOp(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::Build_TotalOp" << endl;
  string MT_name = "";
  for (auto Tens : *Sub_Ops)
    MT_name += Tens->name_;
  Total_Op = make_shared<MultiTensOp<DType>>( MT_name , true);
  Total_Op->initialize( *Sub_Ops );
  Total_Op->get_ctrs_tens_ranges() ;

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::initialize(int nel, int norb, int num_idxs, bool spinfree) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  nact_orb = norb;
  nact_el = nel;
  nidxs = num_idxs;

  GammaMap = make_shared<CombinedGammaMap>();
  set_spin_info();  
  RelGammaMap = make_shared<RelCombinedGammaMap>();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::set_spin_info() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  auto prvec = [&ss](shared_ptr<vector<bool>> invec){ cout << "[ " ; for (auto elem : *invec) { cout << ss(elem) << " " ;} cout << "]" ;};

  spin_max = (nact_orb < nact_el) ? nact_orb : nact_el;
  spin_max = (nact_el - nact_orb) > 0 ?  (nact_el - nact_orb)  : 0;
  
  spin_sectors = spin_sector_gen(nact_el , nact_orb);  
  spin_paths = make_shared<map< pair<int, pair<int,int>> , shared_ptr<vector<shared_ptr<vector<bool>>>> >>();
  for (int ii = 0; ii != (nidxs/2)+1 ; ii++){
    auto all_sc = all_spin_combinations(ii);
    for (pair<int,int> spin_sec : *spin_sectors){
      bool newsector = true;
      for (auto spin_comb : *all_sc) {
        if(sector_constrained_spin_combs( spin_sec.first, spin_sec.second, spin_comb)){
          if(!newsector){
            spin_paths->at(make_pair(ii, spin_sec))->push_back(spin_comb);
          } else {
            auto sscombs_vec = make_shared<vector<shared_ptr<vector<bool>>>>(1,spin_comb); 
            spin_paths->emplace(make_pair(ii, spin_sec), sscombs_vec);
            newsector = false;
          }
        }
      }
    }     
  }
  return; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::bfunc(int xx, shared_ptr<vector<bool>> qvec, int nel, int length, shared_ptr<vector<shared_ptr<vector<bool>>>> all_sc) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  auto xvec = make_shared<vector<bool>>(*qvec);
  if (nel == 0){
    all_sc->push_back(xvec);
    return ;
  } else {
    nel -=1;
    for (int ii = xx; ii != length-nel; ii++){ 
      xvec->at(ii)=true;
      bfunc(ii+1, xvec, nel, length, all_sc);
      xvec->at(ii)=false;
    }
    nel= nel-1;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generates all true/false combinations for a vector which is _half_ the length of the total number of electrons.
//Thish is then used to build spin sector constrained transition lists. Should probably rename
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<vector<shared_ptr<vector<bool>>>>  BraKet<DType>::all_spin_combinations(int length){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };

  auto xvec = make_shared<vector<bool>>(length, false);
  auto all_sc = make_shared<vector<shared_ptr<vector<bool>>>>(0) ;
  for (int nel = 0 ; nel != length+1 ; nel++){
     bfunc( 0, xvec, nel, length, all_sc); 
  }

  auto svecs = make_shared<vector<shared_ptr<vector<bool>>>>(0) ;
  for (auto plus_vec : *all_sc){
    for (auto kill_vec : *all_sc){
      auto buffvec = make_shared<vector<bool>>();
      buffvec->insert(buffvec->end(), plus_vec->begin(), plus_vec->end() );
      buffvec->insert(buffvec->end(), kill_vec->begin(), kill_vec->end() );
      svecs->push_back(buffvec);
    }
  }
  all_sc->swap(*svecs);
  return all_sc ;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Checks that the spin transition pathway never leaves the active spin_sectors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool BraKet<DType>::sector_constrained_spin_combs(int nalpha, int nbeta , shared_ptr<vector<bool>> gamma_comb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  for (int ii = gamma_comb->size()-1; ii != -1; ii-=2 ){
    if (gamma_comb->at(ii)){
      nalpha++;
    } else {
      nbeta++;
    }
    if( (nalpha < 0) || (nbeta < 0) ){
      return false; 
    }

    if (gamma_comb->at(ii-1)){
      nalpha--;
    }else {
      nbeta--;
    }
    if ((nbeta > nact_orb) || (nalpha > nact_orb)){
      return false; 
    }
  }
  
  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Checks that the spin transition pathway never leaves the active spin_sectors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
bool BraKet<DType>::sector_constrained_spin_combs(int nalpha, int nbeta , shared_ptr<vector<string>> gamma_comb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  for (int ii = gamma_comb->size()-1; ii != -1; ii-=2 ){
    if (gamma_comb->at(ii)=="A"){
      nalpha++;
    } else {
      nbeta++;
    }
    if( (nalpha < 0) || (nbeta < 0) ){
      return false; 
    }

    if (gamma_comb->at(ii-1)=="A"){
      nalpha--;
    }else {
      nbeta--;
    }
    if ((nbeta > nact_orb) || (nalpha > nact_orb)){
      return false; 
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//generates all possible spin sectors, should not be needed in Bagel
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<vector<pair<int,int>>> BraKet<DType>::spin_sector_gen(int nact_el, int nact_orb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  int max_spin = nact_el > nact_orb ? nact_orb : nact_el;
  int min_spin = nact_el - nact_orb > 0 ? nact_el- nact_orb : 0;

  auto spin_sectors = make_shared<vector<pair<int,int>>>(0);
  for (int ii = 0 ; ii != max_spin-min_spin+1; ii++)
    spin_sectors->push_back(make_pair(min_spin+ii, max_spin-ii));

  return spin_sectors;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Temporary fix, should either decide on having bools or strings for spins
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::get_string_spin_paths(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };

  spin_paths_str =  make_shared<map< pair< int, pair<int,int> > , shared_ptr<vector<shared_ptr<vector<string>>>> >>();

  for (auto mapit= spin_paths->begin(); mapit!=spin_paths->end(); mapit++){
   auto boolspinvecs = mapit->second;
   auto strspinvecs = make_shared<vector<shared_ptr<vector<string>>>>(0);
   for (auto vec : *boolspinvecs){ 
     auto strspins = make_shared<vector<string>>(0);
     for (auto bspin : *vec) {
       strspins->push_back(ss(bspin)); 
     }
     strspinvecs->push_back(strspins);
    }
   spin_paths_str->emplace(mapit->first, strspinvecs);
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Gamma_WithSpin(shared_ptr<vector<bool>> aops, shared_ptr<vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_WithSpin" << endl;
  auto rdmd  = make_shared<RDMderiv>(); 

  get_string_spin_paths();  
  for (auto spin_sec : *spin_sectors){
    auto spin_lists = spin_paths_str->at(make_pair(nidxs/2, spin_sec));
    for (auto spins : *spin_lists) {
      rdmd->initialize(aops, idxs, spins);
      rdmd->norm_order();

      for ( int kk =0 ; kk != rdmd->allops->size() ; kk++) {
        auto gammas = make_shared<alt_RDMderiv>(); 
        gammas->initialize(rdmd->allops->at(kk), rdmd->allids->at(kk),  rdmd->alldeltas->at(kk),  rdmd->allsigns->at(kk));
        gammas->alt_order();
      
        for ( int ll=0 ; ll != gammas->allids->size(); ll++ ) {
          auto gamma_ranges = gammas->get_ranges(gammas->allids->at(ll)) ;
          auto delta_spins = gammas->strip_delta_spins(gammas->alldeltas->at(ll)) ;
       
          //Change this: spin stripping isn't working, and I also need index ranges for map
           auto deltas_spin_free      = make_shared<vector<pair<string,string>>>(0);
          for (auto  delta : *gammas->alldeltas->at(ll)) {
            delta = make_pair(delta.first.substr(0,delta.first.size()-1), delta.second.substr(0,delta.second.size()-1)) ;
            deltas_spin_free->push_back(delta);
          } 
          
          auto loc_in_map = RelGammaMap->find(make_pair(*gamma_ranges,spin_sec));  
          if ( loc_in_map == RelGammaMap->end() ) {
            auto tmpsigns = make_shared<vector<int>>( 1, gammas->allsigns->at(ll) );
            auto deltas_init       = make_shared<vector<shared_ptr<vector<pair<string,string>>>>>( 1, deltas_spin_free );
            auto deltas_spins_init = make_shared<vector<shared_ptr<vector<pair<string,string>>>>>( 1, delta_spins );
            RelGammaMap->emplace(make_pair(*gamma_ranges,spin_sec), tie( deltas_init, deltas_spins_init, tmpsigns ) );
          } else {
            ( get<0>(loc_in_map->second) )->push_back(deltas_spin_free);
            ( get<1>(loc_in_map->second) )->push_back(delta_spins);
            ( get<2>(loc_in_map->second) )->push_back(gammas->allsigns->at(ll));
          }
        }
      }
    }
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Gamma_SpinFree_New(shared_ptr<vector<bool>> aops, shared_ptr<vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_SpinFree_New" << endl;
  auto rdmd  = make_shared<RDMderiv>(); 
  cout << "aops = " ; for (bool aop : *aops) { cout << aop << " " ; } cout << endl;
  auto aops_buff  = make_shared<vector<bool>>(*aops );

  auto G_to_A_map = make_shared< unordered_map<string, shared_ptr< unordered_map<string, pair<int,int> > >>>(); 
  for (auto range_map_it = Total_Op->combined_ranges->begin() ;  range_map_it !=Total_Op->combined_ranges->end(); range_map_it++){
    auto GGen = make_shared<GammaGenerator>(aops_buff, idxs, G_to_A_map); 
    GGen->add_gamma(range_map_it->first, 1) ;
    GGen->norm_order();
    GGen->alt_order();
  }
 
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){
    cout << endl;
    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;
    for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
      cout <<  A_map_it->first  << "  (" << A_map_it->second.first  << "," <<  A_map_it->second.second << ")" << endl;
    }
  }

  return; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Gamma_SpinFree(shared_ptr<vector<bool>> aops, shared_ptr<vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_SpinFree" << endl;
  auto rdmd  = make_shared<RDMderiv>(); 
  auto spins = make_shared<vector<string>>(vector<string> {"A","A","A","A","A","A","A","A"} );
  cout << "aops = " ; for (bool aop : *aops) { cout << aop << " " ; } cout << endl;
  auto aops_buff  = make_shared<vector<bool>>(*aops );

  rdmd->initialize(aops_buff, idxs, spins);
  aops_buff  = make_shared<vector<bool>>(*aops );
  rdmd->norm_order();
 
  for ( int kk =0 ; kk != rdmd->allops->size() ; kk++) {
    auto gammas = make_shared<alt_RDMderiv>(); 
    gammas->initialize(rdmd->allops->at(kk), rdmd->allids->at(kk),  rdmd->alldeltas->at(kk),  rdmd->allsigns->at(kk));
    gammas->alt_order();

    for ( int ll=0 ; ll != gammas->allids->size(); ll++ ) {
      auto gamma_ranges = make_shared<vector<string>>(gammas->allids->at(ll)->size(), "act") ;
      auto delta_spins = gammas->strip_delta_spins(gammas->alldeltas->at(ll)) ;
   
      //Change this: spin stripping isn't working, and I also need index ranges for map
      auto deltas_spin_free      = make_shared<vector<pair<string,string>>>(0);
      for (auto  delta : *gammas->alldeltas->at(ll)) {
        delta = make_pair(delta.first.substr(0,delta.first.size()-1), delta.second.substr(0,delta.second.size()-1)) ;
        deltas_spin_free->push_back(delta);
      } 
      
      auto loc_in_map = GammaMap->find(*gamma_ranges);  
      if ( loc_in_map == GammaMap->end() ) {
        auto tmpsigns = make_shared<vector<int>>( 1, gammas->allsigns->at(ll) );
        auto deltas_init       = make_shared<vector<shared_ptr<vector<pair<string,string>>>>>( 1, deltas_spin_free );
        auto deltas_spins_init = make_shared<vector<shared_ptr<vector<pair<string,string>>>>>( 1, delta_spins );
        GammaMap->emplace(*gamma_ranges, tie( deltas_init, deltas_spins_init, tmpsigns ) );

      } else {
        ( get<0>(loc_in_map->second) )->push_back(deltas_spin_free);
        ( get<1>(loc_in_map->second) )->push_back(delta_spins);
        ( get<2>(loc_in_map->second) )->push_back(gammas->allsigns->at(ll));
        
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pairs up each gamma term (defined by a spin pathway,and spin sector)  with relevant CtrMultiTensorParts .
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Tensor_Contraction_list_CMTP(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  BK_Compute_List_CMTP = make_shared<map< vector<string>, shared_ptr<vector<pair<shared_ptr<vector<string>>,pair<int,int>>>> >>();
  
  for (auto mapit = GammaMap->begin(); mapit != GammaMap->end(); mapit++) {
    for (int kk = 0 ; kk!= (get<0>(mapit->second))->size(); kk++){
 
      auto Acontrib_loc = Total_Op->CMTP_gamma_contribs->find(tie(mapit->first, *(get<0>(mapit->second))->at(kk), *(get<1>(mapit->second))->at(kk)));  
      if ( Acontrib_loc == Total_Op->CMTP_gamma_contribs->end() ) 
        continue;
       
      auto gamma_factor = make_pair( (get<2>(mapit->second))->at(kk),(get<2>(mapit->second))->at(kk) );
      auto BraKet_Contrib = make_pair(Acontrib_loc->second, gamma_factor);
      
      auto loc_in_map = BK_Compute_List_CMTP->find(mapit->first);  
      if ( loc_in_map == BK_Compute_List_CMTP->end() ) {
        auto init = make_shared<vector<pair< shared_ptr<vector<string>>, pair<int,int> >>>(1,BraKet_Contrib);
        BK_Compute_List_CMTP->emplace(mapit->first, init );

      } else {
        loc_in_map->second->push_back(BraKet_Contrib);
      }
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<std::vector<double>>;
template class BraKet<Tensor_<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
