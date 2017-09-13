#include <src/smith/wicktool/WickUtils.h>
#ifndef __SRC_SMITH_Gamma_Generator_H
#define __SRC_SMITH_Gamma_Generator_H

//#include "WickUtils.h"
using namespace WickUtils;

class GammaMat {

   public :
     std::shared_ptr<std::vector<std::string>> full_id_ranges ;
     std::shared_ptr<std::vector<bool>> full_aops ;
     std::shared_ptr<std::vector<int>> ids_pos ;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos ;
     int my_sign ;

//     GammaMat(full_id_ranges_in, full_aops_in, ids_pos_in, deltas_pos_in, my_sign_in){};
     //full_id_ranges(full_id_ranges_in), full_aops(full_aops_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), my_sign(my_sign_in) {}; 
     GammaMat(){};
     ~GammaMat(){};

};

class GammaGenerator{ 
 
  public : 

    // variables
    std::shared_ptr<std::vector<bool>> orig_aops ;
    std::shared_ptr<std::vector<std::string>> orig_ids ;
    
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::string>>>> G_to_A_map ;
    std::shared_ptr<std::vector<std::shared_ptr<GammaMat>>> Gamma_vec;

    std::shared_ptr<std::map< char, int>>   op_order ;
    std::shared_ptr<std::map< std::string, int>>  idx_order ;
    
    // functions
    GammaGenerator(std::shared_ptr<std::vector<bool>> aops_init, std::shared_ptr<std::vector<std::string>> ids_init);
    ~GammaGenerator(){};
                                                                                                                                                                          
    void norm_order(std::shared_ptr<std::vector<std::shared_ptr<GammaMat>>> gamma_vec );                                                                                  
    
    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaMat>>> gamma_vec  );
    
    std::shared_ptr<pint_vec> Standardize_delta_ordering(std::shared_ptr<pint_vec> deltas_pos ) ;
 
    bool ordering(std::pair<std::string,std::string> a, std::pair<std::string,std::string> b) {
                  return( (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false ); };
   
    std::string get_Aname(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges, 
                          std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos );
   
    std::string get_gamma_name(std::shared_ptr<std::vector<std::string>> full_idx_ranges,  std::shared_ptr<std::vector<bool>> full_aops,
                               std::shared_ptr<std::vector<int>> idxs_pos );

    bool RangeCheck(std::shared_ptr<std::vector<std::string>> full_id_ranges) ;
 };
 #endif
