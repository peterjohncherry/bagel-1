 #ifndef __SRC_SMITH_Gamma_Generator_H
 #define __SRC_SMITH_Gamma_Generator_H

 #include <src/smith/wicktool/WickUtils.h>
 #include <src/smith/wicktool/states_info.h> 
using namespace WickUtils;

class AContribInfo { 

   public : 
     std::vector<std::vector<int>> id_orders;
     std::vector<std::pair<int,int>> factors;

     AContribInfo(std::vector<int> init_order_in , std::pair<int,int> factor_in ){
                  id_orders.push_back(init_order_in);
                  factors.push_back(factor_in); 
     };
     ~AContribInfo(){};

     std::vector<int>   id_order(int qq) {return id_orders[qq]; }; 
     std::pair<int,int> factor(int qq) {return factors[qq]; }; 

};

class GammaInfo {

   public :
     std::shared_ptr<std::vector<std::string>> id_ranges ;
     std::shared_ptr<std::vector<bool>> aops ;
     std::shared_ptr<std::vector<bool>> spins ; 
     std::shared_ptr<std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>>> ssector_sequence;
     std::pair<int,int> factor ;
     std::string name;     
     std::shared_ptr<std::vector<std::string>> one_el_gammas;

     std::shared_ptr<CIVecInfo<double>> Bra_info;
     std::shared_ptr<CIVecInfo<double>> Ket_info;

     int Bra_num;
     int Bra_norb;
     int Bra_nalpha;
     int Bra_nbeta;

     int Ket_num;
     int Ket_norb;
     int Ket_nalpha;
     int Ket_nbeta;

     GammaInfo( std::shared_ptr<CIVecInfo<double>> Bra_info, std::shared_ptr<CIVecInfo<double>> Ket_info, 
                std::shared_ptr<std::vector<bool>> full_aops_vec, std::shared_ptr<std::vector<std::string>> full_idx_ranges,
                std::shared_ptr<std::vector<int>> idxs_pos );

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

     GammaIntermediate(std::shared_ptr<std::vector<std::string>> full_id_ranges_in,
                       std::shared_ptr<std::vector<bool>> full_aops_in,
                       std::shared_ptr<std::vector<int>> ids_pos_in,
                       std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                       int my_sign_in): 
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
    
    // key    : name of this gamma
    // result : map containing names of relevant A-tensors, list of reorderings, and factor for each reordering
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, AContribInfo>> >> G_to_A_map;

    // key    : name of this gamma
    // result : information used here and in compute routines
    std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> Gamma_map;

    // key    : name of A-tensor
    // result : list of reorderings which much be applied to this A-tensor before it is contracted with this gamma.
    //          second part of pair is factor associated with each reordering.
    std::shared_ptr<std::map<std::string, std::vector< std::pair<std::vector<int>, std::pair<int,int>> >> > Aid_orders_map;

    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec;
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec;

    std::shared_ptr<std::map< char, int>>         op_order ;
    std::shared_ptr<std::map< std::string, int>>  idx_order ;
    
    int Ket_num;
    int Bra_num; 

    std::shared_ptr<StatesInfo<double>> TargetStates;

    GammaGenerator(std::shared_ptr<StatesInfo<double>> TargetStates, int Ket_num, int Bra_num,
                   std::shared_ptr<std::vector<bool>> aops_init, std::shared_ptr<std::vector<std::string>> ids_init,
                   std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>> Gamma_map_in, 
                   std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, AContribInfo  >>>> G_to_A_map_in);
    ~GammaGenerator(){};

    void add_gamma(std::shared_ptr<std::vector<std::string>> full_id_ranges_in, int my_sign_in) ;
    
    void norm_order();
    
    void optimized_alt_order();

    void Contract_remaining_indexes(int kk);
   
    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec );
    
    std::shared_ptr<pint_vec> Standardize_delta_ordering(std::shared_ptr<pint_vec> deltas_pos ) ;
 
    bool ordering(std::pair<std::string,std::string> a, std::pair<std::string,std::string> b) {
                  return( (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false ); };
   
    std::string get_Aname(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges, 
                          std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos );
   
    bool RangeCheck(std::shared_ptr<std::vector<std::string>> full_id_ranges) ;
    
    bool Forbidden_Index(std::shared_ptr<std::vector<std::string>> full_id_ranges, int position );

    bool gamma_survives( std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<std::vector<std::string>> id_ranges) ;
     
    std::string get_gamma_name(std::shared_ptr<std::vector<bool>> aops_vec, std::shared_ptr<std::vector<std::string>> full_idx_ranges,
                               std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos );

    static std::shared_ptr<std::vector<std::pair<int,int>>>
           Standardize_delta_ordering_generic(std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos ) ;

    std::vector<int> get_standard_range_order(const std::vector<std::string> &rngs) ;

    std::vector<int> get_position_order(const std::vector<int> &positions) ;

    std::vector<int> get_Aid_order ( const std::vector<int>& id_pos );
 
    std::vector<int> get_standard_order (const std::vector<std::string>& rngs ); 

    std::vector<int> get_standard_idx_order(const std::vector<std::string>& idxs) ;

    std::vector<int> get_standardized_alt_order ( const std::vector<std::string>& rngs ,const std::vector<bool>& aops ) ;

 };
 #endif
