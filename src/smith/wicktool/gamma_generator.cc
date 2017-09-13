#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/gamma_generator.h>

using namespace std; 
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
GammaGenerator::GammaGenerator(shared_ptr<vector<bool>> ac_init,
                               shared_ptr<vector<string>> ids_init){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
  orig_aops = ac_init;
  orig_ids = ids_init; 

  Gamma_vec = make_shared<vector<shared_ptr<GammaMat>>>(0);
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
                      
/////////////////////////////////////////////////////////////////////////////////////////////////////////                                                              
void GammaGenerator::norm_order(shared_ptr<vector<shared_ptr<GammaMat>>> gamma_vec ){                                                                                   
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
//    if (gamma_survives(ids_pos, full_id_ranges)){
//      cout <<  "------------ KEEP ------------ "  << endl;
      cout << "Start Aname  alt = " << Aname_init << endl;
      cout << "Finish Aname alt = " << get_Aname(orig_ids, full_id_ranges, deltas_pos ) << endl;
      cout << "Start Gname  alt = " << Gname_init << endl;
      cout << "Finish Gname alt = " << get_gamma_name( full_id_ranges, orig_aops, ids_pos ) << endl << endl << endl;
//    }                                                       
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

/////////////////////////////////////////////////////////////////////////////
string GammaGenerator::get_Aname(shared_ptr<vector<string>> full_idxs, shared_ptr<vector<string>> full_idx_ranges, 
                                 shared_ptr<vector<pair<int,int>>> all_ctrs_pos ){
/////////////////////////////////////////////////////////////////////////////
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
