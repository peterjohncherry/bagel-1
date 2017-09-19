#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/gamma_generator.h>

using namespace std; 
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaGenerator::GammaGenerator(shared_ptr<vector<bool>> ac_init, shared_ptr<vector<string>> ids_init){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  orig_aops = ac_init;
  orig_ids = ids_init; 

  gamma_vec = make_shared<vector<shared_ptr<GammaMat>>>(0);
  final_gamma_vec = make_shared<vector<shared_ptr<GammaMat>>>(0);
  G_to_A_map = make_shared<std::map<std::string, shared_ptr<vector<std::string>>>>();

  //neeeded to keep ordering of contractions consistent 
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

}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::add_gamma(shared_ptr<vector<string>> full_id_ranges_in, int my_sign_in) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(RangeCheck(full_id_ranges_in)){
    auto ids_pos_init =  make_shared<vector<int>>(full_id_ranges_in->size());
    for ( int ii = 0; ii != full_id_ranges_in->size(); ii++)
        ids_pos_init->at(ii) = ii; 
    
    int my_sign_in = 1;
    auto deltas_pos_in = make_shared<vector<pair<int,int>>>(0);
    
    auto new_gamma = make_shared<GammaMat>(full_id_ranges_in, orig_aops, ids_pos_init, deltas_pos_in, my_sign_in); 
  }

  return; 

}                     
/////////////////////////////////////////////////////////////////////////////////////////////////////////                                                              
//void GammaGenerator::norm_order(shared_ptr<vector<shared_ptr<GammaMat>>> gamma_vec ){                                                                                   
void GammaGenerator::norm_order( ){                                                                                   
//////////////////////////////////////////////////////////////////////////////////////////////////////////                                                              
  int kk = 0;                                                                                                      

  while ( kk != gamma_vec->size()){                                               
    shared_ptr<vector<bool>> full_aops = gamma_vec->at(kk)->full_aops;        
    shared_ptr<vector<int>>  ids_pos = gamma_vec->at(kk)->ids_pos;        
    shared_ptr<vector<string>> full_id_ranges = gamma_vec->at(kk)->full_id_ranges;                  
    shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 

    int  num_pops = ( ids_pos->size()/2 )-1;     
 
    string Aname_init = get_Aname(orig_ids, full_id_ranges, deltas_pos );   
    string Gname_init = get_gamma_name( full_id_ranges, orig_aops, ids_pos );
                                            
    for (int ii = ids_pos->size()-1 ; ii != -1; ii--){            
      if ( ii > num_pops ) {                                      
        if (!full_aops->at(ids_pos->at(ii)))                      
          continue;                                               
                                                                  
        while(full_aops->at( ids_pos->at(ii) )){                  
          for ( int jj = (ii-1); jj != -1 ; jj--) {               
            if (!full_aops->at( ids_pos->at(jj) )){               
              swap( jj, jj+1, kk, gamma_vec);           
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
               swap( jj, jj+1, kk, gamma_vec);          
               break;                                             
             }                                                    
           }                                                      
        }                                                         
      }                                                           
    }
    if (!gamma_survives(ids_pos, full_id_ranges)){
       Contract_remaining_indexes(kk);
    } else {
      final_gamma_vec->push_back(gamma_vec->at(kk));
    } 
    kk++;                                                         
  }                                                               
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
void GammaGenerator::Contract_remaining_indexes(int kk){
//////////////////////////////////////////////////////////////////////////////////////  

  shared_ptr<vector<bool>> full_aops = gamma_vec->at(kk)->full_aops;        
  shared_ptr<vector<int>>  ids_pos = gamma_vec->at(kk)->ids_pos;        
  shared_ptr<vector<string>> full_id_ranges = gamma_vec->at(kk)->full_id_ranges;                  
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 


  //vector of different index ranges, and vectors containing list of index positions
  //associated with each of these ranges
  vector<string> diff_rngs(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> make_ops_pos(0);
  shared_ptr<vector<shared_ptr<vector<int>>>> kill_ops_pos(0);

  // get the positions of the creation and annihilation operators associated with
  // each different range
  for ( int jj = 0;  jj !=full_id_ranges->size() ; jj++){
    int ii = 0;
 
    string rng = full_id_ranges->at(jj);
    do {

      if( rng == diff_rngs[ii] ){
        if ( orig_aops->at(jj) ){
          make_ops_pos->at(ii)->push_back(jj);
        } else {
          kill_ops_pos->at(ii)->push_back(jj);
        }
        break;
      }

      if ( ii == diff_rngs.size()-1 ){
        diff_rngs.push_back(rng);
        if ( orig_aops->at(jj) ){
          auto new_make_vec = make_shared<vector<int>>(1,jj);
          make_ops_pos->push_back(new_make_vec);
        } else {
          auto new_kill_vec = make_shared<vector<int>>(1,jj);
          kill_ops_pos->push_back(new_kill_vec);
        }
        break;
      }
      ii++;
    } while (true);
  } 


//  vector<vector<pair<int,int>>> new_contractions(make_ops_pos->size()); 

  // first index  = range.
  // second index = pair vector defining possible way of contracting indexes with that range.
  // third index  = pair defining contraction of indexes with the same range.
  vector<shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>> new_contractions(diff_rngs.size());


  //CHECK THE CONTRACTIONS ARE STILL IN STANDARD ORDERING AND NOT GETTING DUPLICATES
  for (int ii =0 ; ii != new_contractions.size(); ii++)
    new_contractions[ii] = get_unique_pairs( make_ops_pos->at(ii),  kill_ops_pos->at(ii), make_ops_pos->size());
  
  auto forvec = make_shared<vector<int>>(0,diff_rngs.size()) ;
  auto min = make_shared<vector<int>>(0,diff_rngs.size()) ;
  auto max = make_shared<vector<int>>(0,diff_rngs.size()) ;

  for ( int ii = 0 ; ii != diff_rngs.size(); ii++ ) 
    max->at(ii) = new_contractions[ii]->size()-1;

  do {
    
    auto new_deltas_pos = make_shared<vector<pair<int,int>>>(*deltas_pos);
    for ( int qq = 0; qq != forvec->size(); qq++)
      new_deltas_pos->insert(new_deltas_pos->begin(), new_contractions[qq]->at(forvec->at(qq))->begin(), new_contractions[qq]->at(forvec->at(qq))->begin()); 

    auto new_ids_pos = make_shared<vector<int>>(*ids_pos);

    for ( int qq = 0; qq != new_deltas_pos->size(); qq++ ){ 
      for ( int xx = 0; xx != new_ids_pos->size(); xx++){
        if (new_deltas_pos->at(qq).first == new_ids_pos->at(xx)){
           new_ids_pos->erase(new_ids_pos->begin()+xx); 
        } else if (new_deltas_pos->at(qq).second == new_ids_pos->at(xx)) {
          new_ids_pos->erase(new_ids_pos->begin()+xx); 
        }
      }
    }

    auto new_gamma = make_shared<GammaMat>(full_id_ranges, orig_aops, new_ids_pos, new_deltas_pos, gamma_vec->at(kk)->my_sign); 

    final_gamma_vec->push_back(new_gamma);

  } while (fvec_cycle( forvec, max , min)) ;
 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
//void GammaGenerator::alt_order(shared_ptr<vector<shared_ptr<GammaMat>>> gamma_vec ){
void GammaGenerator::alt_order(){
//////////////////////////////////////////////////////////////////////////////////////  

  auto even = []( int pos) { return ( 0 == pos % 2);}; 
  int kk = 0;
  while ( kk != gamma_vec->size()){
    shared_ptr<vector<bool>> full_aops = gamma_vec->at(kk)->full_aops;        
    shared_ptr<vector<int>>  ids_pos = gamma_vec->at(kk)->ids_pos;        
    shared_ptr<vector<string>> full_id_ranges = gamma_vec->at(kk)->full_id_ranges;                  
    shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 
    int  num_pops = ( ids_pos->size()/2 )-1;     
 
    string Aname_init = get_Aname( orig_ids, full_id_ranges, deltas_pos );   
    string Gname_init = get_gamma_name( full_id_ranges, orig_aops, ids_pos );

    for (int ii = orig_ids->size()-1 ; ii != -1; ii-- ){
      if (!even(ii)) {
        if (!full_aops->at(ii))
          continue;
 
        while(full_aops->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (!full_aops->at(jj)){
               swap(jj, jj+1, kk, gamma_vec);
               break;
             }
          }
        }
      } else if (even(ii)){
          if (full_aops->at(ii))
            continue;
        while(!full_aops->at(ii)){
          for ( int jj = (ii-1); jj != -1 ; jj--) {
             if (full_aops->at(jj)){
               swap(jj, jj+1, kk, gamma_vec);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
//Swaps indexes round, flips sign, and if ranges are the same puts new density matrix in the list. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
void GammaGenerator::swap( int ii, int jj, int kk, shared_ptr<vector<shared_ptr<GammaMat>>> gamma_vec  ){                                                  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
  
  shared_ptr<vector<bool>> full_aops = gamma_vec->at(kk)->full_aops;        
  shared_ptr<vector<int>>  ids_pos = gamma_vec->at(kk)->ids_pos;        
  shared_ptr<vector<string>> full_id_ranges = gamma_vec->at(kk)->full_id_ranges;                  
  shared_ptr<vector<pair<int,int>>> deltas_pos = gamma_vec->at(kk)->deltas_pos; 


  int idx_buff = ids_pos->at(ii);                                                                                                   
  ids_pos->at(ii) = ids_pos->at(jj);                                                                                                
  ids_pos->at(jj) = idx_buff;                                                                                                       
                                                                                                                                    
  if ( (full_id_ranges->at(ids_pos->at(jj)) == full_id_ranges->at(ids_pos->at(ii))) &&                                                
       (orig_aops->at(ids_pos->at(ii)) !=  orig_aops->at(ids_pos->at(jj))) )  {                                                   

    auto new_deltas_tmp = make_shared<pint_vec>(*deltas_pos);                                                                       
    new_deltas_tmp->push_back(make_pair(ids_pos->at(jj), ids_pos->at(ii)));                                                         
    auto new_deltas = Standardize_delta_ordering( new_deltas_tmp ) ;                                                                

    int new_sign = gamma_vec->at(kk)->my_sign;                                                                                               

    auto new_ids_pos = make_shared<vector<int>>();
    for( int qq = 0 ; qq !=ids_pos->size() ; qq++)
      if ( (qq !=ii) && (qq!=jj))
        new_ids_pos->push_back(ids_pos->at(qq));

    auto new_gamma = make_shared<GammaMat>(full_id_ranges, orig_aops, new_ids_pos, new_deltas, new_sign); 
    gamma_vec->push_back(new_gamma);

  }                                                                                                                                 
  gamma_vec->at(kk)->my_sign *= -1;
  return ;                                                                                                                          
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string GammaGenerator::get_Aname(shared_ptr<vector<string>> full_idxs, shared_ptr<vector<string>> full_idx_ranges, 
                                 shared_ptr<vector<pair<int,int>>> all_ctrs_pos ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string  name = "";
  for(string idx : *full_idxs)
    name += idx;
  name+="_"; 

  for(string idx_range : *full_idx_ranges)
    name += idx_range[0];

  if (all_ctrs_pos->size() !=0 ){
    name+="_"; 
    for(pair<int,int> delta : *all_ctrs_pos)
      name += to_string(delta.first)+to_string(delta.second);
  }
  return name;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string GammaGenerator::get_gamma_name(shared_ptr<vector<string>> full_idx_ranges,  shared_ptr<vector<bool>> aops_vec,
                                      shared_ptr<vector<int>> idxs_pos ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string  name = "";
  
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
  
  return name;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<pint_vec> GammaGenerator::Standardize_delta_ordering(shared_ptr<pint_vec> deltas_pos ) {
//////////////////////////////////////////////////////////////////////////////////////////////////

   bool(GammaGenerator::*orderfunc)(pair<string,string>, pair<string,string>) = &GammaGenerator::ordering;

   auto dtmp = make_shared<pint_vec>(*deltas_pos);

   if ( orig_aops->at( dtmp->back().second )){
      dtmp->back() = make_pair( deltas_pos->back().second, deltas_pos->back().first );	
   }

   //sort(dtmp->begin(), dtmp->end() );
   // note that by construction only the last element of dtmp could be in the wrong place.
   int ii = 0;
   while (ii < dtmp->size()-1) {
     if ( op_order->at(orig_ids->at(dtmp->back().first)[0])  < op_order->at(orig_ids->at(dtmp->at(ii).first)[0]) ){

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

///////////////////////////////////////////////////////////////////////////////                                           
//Most basic constraint of all active
//should instead take function
///////////////////////////////////////////////////////////////////////////////                                   
bool GammaGenerator::gamma_survives(shared_ptr<vector<int>> ids_pos, shared_ptr<vector<string>> id_ranges) {
///////////////////////////////////////////////////////////////////////////////

   for (int pos : *ids_pos )
     if (id_ranges->at(pos) != "act")
        return false;
 
   return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
bool GammaGenerator::RangeCheck(shared_ptr<vector<string>> full_id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////

   vector<string> diff_rngs(0);
   vector<int> updown(0);

   for ( int jj = 0;  jj !=full_id_ranges->size() ; jj++){
     int ii = 0;
 
     string rng = full_id_ranges->at(jj);
     do {

       if(rng == diff_rngs[ii]){
         if ( orig_aops->at(jj)){
           updown[jj]+=1;
         } else {
           updown[jj]-=1;
         }
         break;
       }

       if ( ii == diff_rngs.size()-1){
         diff_rngs.push_back(rng);
         if ( orig_aops->at(jj)){
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
#endif
//      cout <<  "------------ KEEP ------------ "  << endl;
//      cout << "Start Aname  alt = " << Aname_init << endl;
//      cout << "Finish Aname alt = " << get_Aname(orig_ids, full_id_ranges, deltas_pos ) << endl;
//      cout << "Start Gname  alt = " << Gname_init << endl;
//      cout << "Finish Gname alt = " << get_gamma_name( full_id_ranges, orig_aops, ids_pos ) << endl << endl << endl;
//    }                                                       

