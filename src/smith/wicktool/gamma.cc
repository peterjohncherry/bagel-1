#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/gamma.h>

//#include "WickUtils.h"
//#include "gamma.h"

using namespace std; 
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RDMderiv_new::initialize(shared_ptr<vector<bool>> ac_init,
                              shared_ptr<vector<string>> ids_init, 
                              shared_ptr<pint_vec> deltas_pos_init,
                              int sign){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  full_ids = ids_init; 
  full_aops = ac_init;
  spinfree = false; 

  vector<bool>  unc_get(ids_init->size(), true);
  for(int ii =0 ; ii != deltas_pos_init->size() ; ii++) {
    unc_get[deltas_pos_init->at(ii).first] = false;
    unc_get[deltas_pos_init->at(ii).second] = false;
  }

  shared_ptr<vector<int>> unc_pos(0);
  shared_ptr<vector<string>> unc_ids(0);
  for (int ii = 0 ; ii != unc_get.size(); ii++ ) {
    if(unc_get[ii]){
     unc_pos->push_back(ii);
     unc_ids->push_back(full_ids->at(ii));
    }
  }

  ids_pos_all = make_shared<vector<shared_ptr<vector<int>>>>(1, unc_pos);
  deltas_pos_all = make_shared<vector<shared_ptr<pint_vec>>>(1, deltas_pos_init); 
  signs_all  = make_shared<vector<int>>(1,sign);

  auto opname = unc_ids->at(0)[0];
  op_order = make_shared<map< char, int>>();
  int ii =0;

  op_order->emplace (opname, ii);
  for ( auto elem : *unc_ids ){
    if ( opname != elem[0] ){
      opname = elem[0];
      op_order->emplace(opname, ++ii);
    }
  }

  ii=0;
  idx_order = make_shared<map< string, int>>();
  for ( auto elem : *unc_ids ){
     idx_order->emplace(elem, ii);
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RDMderiv_new::initialize(shared_ptr<vector<bool>> ac_init,
                              shared_ptr<vector<string>> ids_init, 
                              shared_ptr<vector<string>> id_ranges){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  full_ids = ids_init; 
  full_id_ranges = id_ranges;
  full_aops = ac_init;
  
  spinfree = false; 
  auto  ids_pos_init = make_shared<vector<int>>(id_ranges->size());
  for (int ii = 0 ; ii != id_ranges->size() ; ii++) ids_pos_init->at(ii) = ii;
 
  auto nodel = make_shared<pint_vec>(0);
  ids_pos_all = make_shared<vector<shared_ptr<vector<int>>>>(1, ids_pos_init);
  deltas_pos_all = make_shared<vector<shared_ptr<pint_vec>>>(1, nodel); 
  signs_all  = make_shared<vector<int>>(1,1);

  //neeeded to keep ordering of contractions consistent 
  auto opname = full_ids->at(0)[0];
  op_order = make_shared<map< char, int>>();
  int ii =0;
  op_order->emplace (opname, ii);
  for ( auto elem : *full_ids ){
    if ( opname != elem[0] ){
      opname = elem[0];
      op_order->emplace(opname, ++ii);
    }
  }

  ii=0;
  idx_order = make_shared<map< string, int>>();
  for ( auto elem : *full_ids ){
     idx_order->emplace(elem, ii);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
void RDMderiv_new::swap(shared_ptr<vector<int>> ids_pos,                                                                            
                        shared_ptr<pint_vec> deltas_pos, int ii, int jj, int kk ){                                                  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
  cout << "RDMderiv_new::swap" << endl;                                                                                             
                                                                                                                                    
  int idx_buff = ids_pos->at(ii);                                                                                                   
                                                                                                                                    
  ids_pos->at(ii) = ids_pos->at(jj);                                                                                                
  ids_pos->at(jj) = idx_buff;                                                                                                       
                                                                                                                                    
                                                                                                                                    
  if ( full_id_ranges->at(ids_pos->at(jj)) == full_id_ranges->at(ids_pos->at(ii)) &&                                                
       full_aops->at( ids_pos->at(ii) ) !=  full_aops->at( ids_pos->at(jj) ) )  {                                                   
    cout << "adding contraction" <<endl;                                                                                            
    auto new_deltas_tmp = make_shared<pint_vec>(*deltas_pos);                                                                       
    new_deltas_tmp->push_back(make_pair(ids_pos->at(jj), ids_pos->at(ii)));                                                         
    auto new_deltas = Standardize_delta_ordering( new_deltas_tmp ) ;                                                                
    deltas_pos_all->push_back(new_deltas);                                                                                          
                                                                                                                                    
    auto new_ids_pos =  make_shared<vector<int>>(*ids_pos);                                                                         
    new_ids_pos->erase(new_ids_pos->begin()+(jj-1));                                                                                
    new_ids_pos->erase(new_ids_pos->begin()+(jj-1));                                                                                
    ids_pos_all->push_back(new_ids_pos);                                                                                            
    int new_sign = signs_all->at(kk);                                                                                               
    signs_all->push_back(new_sign);                                                                                                 
  }                                                                                                                                 
  signs_all->at(kk) = -1 * signs_all->at(kk);                                                                                       
  return ;                                                                                                                          
}                                                                                                                                   
/////////////////////////////////////////////////////////////////////////                                                           
void RDMderiv_new::generic_reordering(shared_ptr<vector<int>> new_order ){                                                          
/////////////////////////////////////////////////////////////////////////                                                           
  int kk = 0;
  while ( kk != ids_pos_all->size()){
    auto ids_pos = ids_pos_all->at(kk);
    auto deltas  = deltas_pos_all->at(kk);

    for (int ii = ids_pos->size()-1 ; ii != -1; ii-- ){
      if (ids_pos->at(ii) == new_order->at(ii))
        continue;
 
      while( ids_pos->at(ii) != new_order->at(ii) ){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
           if (ids_pos->at(jj) == new_order->at(ii)){
             swap( ids_pos, deltas, jj, jj+1, kk);
             break;
           }
        }
      }
    }
    cout << endl;
    kk++;
  }
  return;
}

/////////////////////////////////////////////////////                                                              
void RDMderiv_new::norm_order(){                                                                                   
//////////////////////////////////////////////////////                                                               
  cout << "RDMderiv_New::norm_order() " << endl;                                                                   
  int kk = 0;                                                                                                      
  cout << "full_aops = "  ; for (bool aops_ : *full_aops ) { cout << aops_ << " " ; } cout << endl;                
  while ( kk != ids_pos_all->size()){                                                                                 
    auto ids_pos = ids_pos_all->at(kk);        
    auto deltas_pos  = deltas_pos_all->at(kk); 
    int  num_pops = (ids_pos->size()/2)-1;     
                                                                                                                   
    cout << "full_aops->size() = "; cout.flush(); cout  << full_aops->size() << endl;                                
                                                                                                                     
    for (int ii = ids_pos->size()-1 ; ii != -1; ii--){            
      cout << "ii = " << ii << endl;                              
      if ( ii > num_pops ) {                                      
        if (!full_aops->at(ids_pos->at(ii)))                      
          continue;                                               
                                                                  
        while(full_aops->at( ids_pos->at(ii) )){                  
          for ( int jj = (ii-1); jj != -1 ; jj--) {               
            if (!full_aops->at( ids_pos->at(jj) )){               
              swap( ids_pos, deltas_pos, jj, jj+1, kk);           
              break;                                              
            }                                                     
          }                                                       
        }                                                         
      } else if (ii<=num_pops){                                   
                                                                  
          if (full_aops->at(ids_pos->at(ii)))                     
            continue;                                             
                                                                  
        while(!full_aops->at(ids_pos->at(ii))){                   
          for ( int jj = (ii-1); jj != -1 ; jj--) {               
             if(full_aops->at(ids_pos->at(jj)) ){                 
               swap( ids_pos, deltas_pos, jj, jj+1, kk);          
               break;                                             
             }                                                    
           }                                                      
        }                                                         
      }                                                           
    }                                                             
    kk++;                                                         
  }                                                               
  return;
}
///////////////////////////////////////////////////////                           
void RDMderiv_new::alt_order(){
///////////////////////////////////////////////////////  

  auto even = []( int pos) { return ( 0 == pos % 2);}; 
  int kk = 0;
  while ( kk != ids_pos_all->size()){
    auto ids_pos = ids_pos_all->at(kk);
    auto deltas_pos  = deltas_pos_all->at(kk);

    for (int ii = ids_pos->size()-1 ; ii != -1; ii-- ){
      if (!even(ii)) {
        if (!full_aops->at(ids_pos->at(ii)))
          continue;
 
        while(full_aops->at(ids_pos->at(ii))){ 
          for ( int jj = (ii-1); jj != -1 ; jj--) {
            if (full_aops->at(ids_pos->at(ii)) ){
              swap( ids_pos, deltas_pos, jj, jj+1, kk);
              break;
            }
          }
        }
      } else if (even(ii)){
        if(full_aops->at(ids_pos->at(ii))) 
            continue;
        while(!full_aops->at(ids_pos->at(ii))){ 
          for ( int jj = (ii-1); jj != -1 ; jj--) {
            if (!full_aops->at(ids_pos->at(ii)) ) {
              swap( ids_pos, deltas_pos,  jj, jj+1, kk);
              break;
            }
          }
        }        
      }
    }
    kk++;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec> RDMderiv_new::Standardize_delta_ordering(shared_ptr<pint_vec> deltas_pos ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   bool(RDMderiv_new::*orderfunc)(pair<string,string>, pair<string,string>) = &RDMderiv_new::ordering;

   auto dtmp = make_shared<pint_vec>(*deltas_pos);

   if ( full_aops->at( dtmp->back().second )){
      dtmp->back() = make_pair( deltas_pos->back().second, deltas_pos->back().first );	
   }

   //sort(dtmp->begin(), dtmp->end() );
   // note that by construction only the last element of dtmp could be in the wrong place.
   int ii = 0;
   while (ii < dtmp->size()-1) {
     if ( op_order->at(full_ids->at(dtmp->back().first)[0])  < op_order->at(full_ids->at(dtmp->at(ii).first)[0]) ){

       auto dtmp2 = make_shared<pint_vec>();
       
       for (int qq= 0 ; qq != ii ; qq++)
          dtmp2->push_back(dtmp->at(qq));
       
       dtmp2->push_back(dtmp->back());
       
       for (int qq= ii ; qq !=dtmp->size()-1 ; qq++)
          dtmp2->push_back(dtmp->at(qq));
         
       dtmp = dtmp2;
       break;
      }
      ii++;
  }
  deltas_pos = dtmp;
  return deltas_pos;   
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RDMderiv::initialize(shared_ptr<vector<bool>> ac_init,
                          shared_ptr<vector<string>> ids_init, 
                          shared_ptr<pstr_vec> deltas_init,
                          int sign){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  orig_ids = ids_init; 
  aops = ac_init;
  spinfree = false; 
  allops = make_shared<vector<shared_ptr<vector<bool>>>>(1, ac_init);
  allids = make_shared<vector<shared_ptr<vector<string>>>>(1, ids_init);
  alldeltas = make_shared<vector<shared_ptr<pstr_vec>>>(1,deltas_init); 
  allsigns  = make_shared<vector<int>>(1,sign);

  auto opname = orig_ids->at(0)[0];
  op_order = make_shared<map< char, int>>();
  int ii =0;

  op_order->emplace (opname, ii);
  for ( auto elem : *orig_ids ){
    if ( opname != elem[0] ){
      opname = elem[0];
      op_order->emplace(opname, ++ii);
    }
  }

  ii=0;
  idx_order = make_shared<map< string, int>>();
  for ( auto elem : *orig_ids ){
     idx_order->emplace(elem, ii);
  }


  aop_map = make_shared<map< string, bool>>();
  for ( int jj = 0 ;  jj !=orig_ids->size(); jj++ )
     aop_map->emplace(orig_ids->at(jj), aops->at(jj));

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RDMderiv::initialize(shared_ptr<vector<bool>> ac_init,
                          shared_ptr<vector<string>> ids_init, 
                          shared_ptr<vector<string>> ids_spins){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  orig_ids = ids_init; 

  spinfree = false; 
 
  auto spin_ids = make_shared<vector<string>>(*ids_init);
  for(int ii = 0; ii != ids_spins->size(); ii++) 
    spin_ids->at(ii)+=ids_spins->at(ii);

  auto nodel = make_shared<pstr_vec>(0);
  allops = make_shared<vector<shared_ptr<vector<bool>>>>(1, ac_init);
  allids = make_shared<vector<shared_ptr<vector<string>>>>(1, spin_ids);
  alldeltas = make_shared<vector<shared_ptr<pstr_vec>>>(1,nodel); 
  allsigns  = make_shared<vector<int>>(1,1);
 
  auto opname = orig_ids->at(0)[0];
  op_order = make_shared<map< char, int>>();
  int ii =0;
  op_order->emplace (opname, ii);
  for ( auto elem : *orig_ids ){
    if ( opname != elem[0] ){
      opname = elem[0];
      op_order->emplace(opname, ++ii);
    }
  }

  ii=0;
  idx_order = make_shared<map< string, int>>();
  for ( auto elem : *orig_ids ){
     idx_order->emplace(elem, ii);
  }

  aop_map = make_shared<map< string, bool>>();
  for ( int jj = 0 ;  jj !=orig_ids->size(); jj++ )
    aop_map->emplace(spin_ids->at(jj), ac_init->at(jj));

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RDMderiv::swap(shared_ptr<vector<bool>> ac, shared_ptr<vector<string>> ids, shared_ptr<pstr_vec> dlist,
                  int ii, int jj, int kk ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool a_buff = ac->at(ii);
  string idx_buff = ids->at(ii);;

  ac->at(ii) = ac->at(jj);
  ids->at(ii) = ids->at(jj);
  ac->at(jj) = a_buff;
  ids->at(jj) = idx_buff;

  if (spinfree || (ids->at(jj).back() == ids->at(ii).back())){
    auto new_deltas_tmp = make_shared<pstr_vec>(*dlist); 
    new_deltas_tmp->push_back(make_pair(ids->at(jj), ids->at(ii)));
    auto new_deltas = Standardize_delta_ordering( new_deltas_tmp ) ;
    alldeltas->push_back(new_deltas);
    
    auto new_ac =  make_shared<vector<bool>>(*ac);   
    new_ac->erase(new_ac->begin()+(jj-1)); 
    new_ac->erase(new_ac->begin()+(jj-1)); 
    allops->push_back(new_ac);
    
    auto new_ids =  make_shared<vector<string>>(*ids); 
    new_ids->erase(new_ids->begin()+(jj-1)); 
    new_ids->erase(new_ids->begin()+(jj-1)); 
    allids->push_back(new_ids); 
    int new_sign = allsigns->at(kk); 
    allsigns->push_back(new_sign);
  }
  allsigns->at(kk) = -1 * allsigns->at(kk);
  return ;
}
/////////////////////////////////////////////////////////////////////////
void RDMderiv::generic_reordering(shared_ptr<vector<bool>> new_order ){
/////////////////////////////////////////////////////////////////////////  
  int kk = 0;
  while ( kk != allops->size()){
    auto ac  = allops->at(kk);
    auto ids = allids->at(kk);
    auto deltas  = alldeltas->at(kk);
    int num_pops = (ids->size()/2)-1;

    for (int ii = ids->size()-1 ; ii != -1; ii-- ){
      if (ac->at(ii) == new_order->at(ii))
        continue;
 
      while( ac->at(ii) != new_order->at(ii) ){
        for ( int jj = (ii-1); jj != -1 ; jj--) {
           if (ac->at(jj) == new_order->at(ii)){
             swap(ac, ids, deltas, jj, jj+1, kk);
             break;
           }
        }
      }
    }
    cout << endl;
    kk++;
  }
  return;
}
//////////////////////////////////////////////////////
void RDMderiv::norm_order(){
//////////////////////////////////////////////////////  

  int kk = 0;
  while ( kk != allops->size()){
    auto ac  = allops->at(kk);
    auto ids = allids->at(kk);
    auto deltas  = alldeltas->at(kk);
    int num_pops = (ids->size()/2)-1;

    for (int ii = ids->size()-1 ; ii != -1; ii-- ){
      if (ii > num_pops) {
        if (!ac->at(ii))
          continue;
 
        while(ac->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (!ac->at(jj)){
               swap(ac, ids, deltas, jj, jj+1, kk);
               break;
             }
          }
        }
      } else if (ii<=num_pops){
          if (ac->at(ii))
            continue;
        while(!ac->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (ac->at(jj)){
               swap(ac, ids, deltas, jj, jj+1, kk);
               break;
             }
           }
        }        
      }
    }
    kk++;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pstr_vec> RDMderiv::Standardize_delta_ordering(shared_ptr<pstr_vec> delta_ids ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   bool(RDMderiv::*orderfunc)(pair<string,string>, pair<string,string>) = &RDMderiv::ordering;

   auto dtmp = make_shared<pstr_vec>(*delta_ids);

   if ( aop_map->at( dtmp->back().second )){
      dtmp->back() = make_pair( delta_ids->back().second, delta_ids->back().first );	
   }

   //sort(dtmp->begin(), dtmp->end() );

   // note that by construction only the last element of dtmp could be in the wrong place.
   int ii = 0;
   while (ii < dtmp->size()-1) {
     if (op_order->at(dtmp->back().first[0]) < op_order->at(dtmp->at(ii).first[0])){

       auto dtmp2 = make_shared<pstr_vec>();
       
       for (int qq= 0 ; qq != ii ; qq++)
          dtmp2->push_back(dtmp->at(qq));
       
       dtmp2->push_back(dtmp->back());
       
       for (int qq= ii ; qq !=dtmp->size()-1 ; qq++)
          dtmp2->push_back(dtmp->at(qq));
         
       dtmp = dtmp2;
       break;
      }
      ii++;
  }
   delta_ids = dtmp;
   return delta_ids;   
}
///////////////////////////////////////////////////////
void alt_RDMderiv::alt_order(){
///////////////////////////////////////////////////////  

  auto even = []( int pos) { return ( 0 == pos % 2);}; 
  int kk = 0;
  while ( kk != allops->size()){
    auto ac  = allops->at(kk);
    auto ids = allids->at(kk);
    auto deltas  = alldeltas->at(kk);

    for (int ii = ids->size()-1 ; ii != -1; ii-- ){
      if (!even(ii)) {
        if (!ac->at(ii))
          continue;
 
        while(ac->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (!ac->at(jj)){
               swap(ac, ids, deltas, jj, jj+1, kk);
               break;
             }
          }
        }
      } else if (even(ii)){
          if (ac->at(ii))
            continue;
        while(!ac->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (ac->at(jj)){
               swap(ac, ids, deltas,  jj, jj+1, kk);
               break;
             }
           }
        }        
      }
    }
    kk++;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<string>> alt_RDMderiv::get_ranges(shared_ptr<vector<string>> spin_ids ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// must be a better way... to_string is returning numbers, and I need ranges not spins for map
 
   auto sonv = []( char inch) {
      if (inch == 'A' ) {
        return "actA";
      } else if (inch == 'B' ) {
        return "actB";
      }
   };

   auto spin_vec = make_shared<vector<string>>(0);
   for (auto id : *spin_ids)
      spin_vec->push_back(sonv(id.back()));

  return spin_vec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pstr_vec> alt_RDMderiv::strip_delta_spins(shared_ptr<pstr_vec> deltas_spin_ids ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   // must be a better way... to_string is returning numbers
   auto sonv = []( char inch) {
      if (inch == 'A' ) {
        return "A";
      } else if (inch == 'B' ) {
        return "B";
      }
   };

   auto delta_spins = make_shared<pstr_vec>(0);
   auto new_delta_ids = make_shared<pstr_vec>(0);
   for (auto delta_id : *deltas_spin_ids){
    delta_spins->push_back( make_pair( sonv(delta_id.first.back()), sonv(delta_id.second.back()) ) );
    delta_id.first.pop_back();
    delta_id.second.pop_back();
    
    new_delta_ids->push_back( make_pair( delta_id.first, delta_id.second ) );
   }

   deltas_spin_ids = new_delta_ids;
  return delta_spins;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec> alt_RDMderiv::get_spinsector_path(shared_ptr<pstr_vec> gamma_spin_ids , pair<int,int> spinsector ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note the spin sector is that of the  Bra, hence daggered ops will annihilate


  pair<int,int> sector = spinsector;
  auto ssector_path = make_shared<pint_vec>(0);
  for (auto exc : *gamma_spin_ids) {
    if (exc.first == exc.second ){
      ssector_path->push_back(sector);
    } else if( exc.first == "A"){
      ssector_path->push_back(make_pair(--sector.first, ++sector.second));
    } else {
      ssector_path->push_back(make_pair(++sector.first, --sector.second));
    }
  }

  return ssector_path;
}


#endif
