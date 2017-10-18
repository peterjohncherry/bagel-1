 #ifndef __SRC_SMITH_Gamma_Generator_H
 #define __SRC_SMITH_Gamma_Generator_H

 #include <src/smith/wicktool/WickUtils.h>
 #include "WickUtils.h"
 #include <unordered_map>
using namespace WickUtils;


class GammaInfo {

   public :
     std::shared_ptr<std::vector<std::string>> id_ranges ;
     std::shared_ptr<std::vector<bool>> aops ;
     std::shared_ptr<std::vector<bool>> spins ; 
     std::shared_ptr<std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>>> ssector_sequence;
     std::pair<int,int> factor ;
     std::string name;     
     std::shared_ptr<std::vector<std::string>> one_el_gammas;

     GammaInfo(std::shared_ptr<std::vector<std::string>> id_ranges_in, std::shared_ptr<std::vector<bool>> aops_in) :
               id_ranges(id_ranges_in), aops(aops_in) {}; 
     GammaInfo(std::shared_ptr<std::vector<bool>> full_aops_vec, std::shared_ptr<std::vector<std::string>> full_idx_ranges, std::shared_ptr<std::vector<int>> idxs_pos);


     GammaInfo(){};
     ~GammaInfo(){};

};

class GammaIntermediate {

   public :
     std::shared_ptr<std::vector<std::string>> full_id_ranges ;
     std::shared_ptr<std::vector<bool>> full_aops ;
     std::shared_ptr<std::vector<int>> ids_pos ;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos ;
     int my_sign ;

     GammaIntermediate(std::shared_ptr<std::vector<std::string>> full_id_ranges_in, std::shared_ptr<std::vector<bool>> full_aops_in, std::shared_ptr<std::vector<int>> ids_pos_in,
              std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in, int my_sign_in): 
     full_id_ranges(full_id_ranges_in), full_aops(full_aops_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), my_sign(my_sign_in) {}; 

     GammaIntermediate(std::shared_ptr<std::vector<std::string>> full_id_ranges_in, int my_sign_in) {};
     GammaIntermediate(){};
     ~GammaIntermediate(){};

};

class GammaGenerator{ 
  friend GammaInfo; 
  public : 

    // variables
    std::shared_ptr<std::vector<bool>> orig_aops ;
    std::shared_ptr<std::vector<std::string>> orig_ids ;
    
    std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< std::unordered_map<std::string, std::pair<int,int> > >>> G_to_A_map;
    std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< GammaInfo >>> Gamma_map;

    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec;
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec;

    std::shared_ptr<std::map< char, int>>   op_order ;
    std::shared_ptr<std::map< std::string, int>>  idx_order ;

    // functions
    GammaGenerator(std::shared_ptr<std::vector<bool>> aops_init, std::shared_ptr<std::vector<std::string>> ids_init,
                   std::shared_ptr<std::unordered_map<std::string, std::shared_ptr<GammaInfo>>> Gamma_map_in, 
                   std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< std::unordered_map<std::string, std::pair<int,int> > >>> G_to_A_map_init );
    ~GammaGenerator(){};

    void add_gamma(std::shared_ptr<std::vector<std::string>> full_id_ranges_in, int my_sign_in) ;
    
    void norm_order();
    
    void alt_order();
 
    void optimized_alt_order();

    void Contract_remaining_indexes(int kk);
   
    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec );
    
    std::shared_ptr<pint_vec> Standardize_delta_ordering(std::shared_ptr<pint_vec> deltas_pos ) ;
 
    bool ordering(std::pair<std::string,std::string> a, std::pair<std::string,std::string> b) {
                  return( (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false ); };
   
    std::string get_Aname(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges, 
                          std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos );
   
    static std::string get_gamma_name(std::shared_ptr<std::vector<std::string>> full_idx_ranges,  std::shared_ptr<std::vector<bool>> full_aops,
                               std::shared_ptr<std::vector<int>> idxs_pos );

    bool RangeCheck(std::shared_ptr<std::vector<std::string>> full_id_ranges) ;
    
    bool gamma_survives( std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<std::vector<std::string>> id_ranges) ;
     
    std::string get_gamma_name(std::shared_ptr<std::vector<bool>> aops_vec, std::shared_ptr<std::vector<std::string>> full_idx_ranges,
                               std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos );

    static std::shared_ptr<std::vector<std::pair<int,int>>>
           Standardize_delta_ordering_generic(std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos ) ;

    std::vector<int> sort_ranges(const std::vector<std::string> &rngs) ;

    std::vector<int> get_standard_order (const std::vector<std::string>& rngs ); 

    std::vector<int> get_standardized_alt_order ( const std::vector<std::string>& rngs ,const std::vector<bool>& aops ) ;

 };
 #endif
