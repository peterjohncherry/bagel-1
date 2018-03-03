#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>

using namespace WickUtils;

class AContribInfo {

   public :
     std::string name_;
     std::vector<std::vector<int>> id_orders;
     std::vector<std::pair<double,double>> factors;
     int total_uses_;
     int remaining_uses_;

     AContribInfo( std::string name,  std::vector<int> init_order_in , std::pair<double,double> factor_in ):
                   name_(name), id_orders(std::vector<std::vector<int>>(1,init_order_in)),
                   factors(std::vector<std::pair<double,double>>(1,factor_in)),
                   total_uses_(1), remaining_uses_(1) {};
     ~AContribInfo(){};

     std::vector<int>   id_order(int qq) {return id_orders[qq]; };
     std::pair<double,double> factor(int qq) {return factors[qq]; };

     std::string name() {return name_ ;}

     int total_uses() const { return total_uses_ ; }
     int remaining_uses() const { return remaining_uses_ ; }

     void increase_total_uses() { total_uses_+=1; }
     void increase_remaining_uses() { remaining_uses_+=1; }
     void decrease_total_uses() { total_uses_-=1; }
     void decrease_remaining_uses() { remaining_uses_-=1; }

};

class GammaInfo {

   private :
     int order_ ;

     //Follows convention : <I | i j | K > < K | .........| J >
     std::shared_ptr<CIVecInfo<double>> Bra_info_;      // < I |
     std::shared_ptr<CIVecInfo<double>> Ket_info_;      // | J >
     std::shared_ptr<CIVecInfo<double>> prev_Bra_info_; // < K |

     std::shared_ptr<std::vector<std::string>> id_ranges_ ;
     std::shared_ptr<std::vector<std::string>> sigma_id_ranges_ ;

     std::shared_ptr<std::vector<bool>> aops_ ;
     std::pair<int,int> factor_ ;

     std::string name_;
     std::string sigma_name_;

     std::vector<std::string> prev_gammas_;
     std::vector<std::string> prev_sigmas_;

   public :
     GammaInfo( std::shared_ptr<CIVecInfo<double>> Bra_info, std::shared_ptr<CIVecInfo<double>> Ket_info,
                std::shared_ptr<const std::vector<bool>> full_aops_vec, std::shared_ptr<const std::vector<std::string>> full_idx_ranges,
                std::shared_ptr<std::vector<int>> idxs_pos, std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> Gamma_map_in );
     GammaInfo(){};
     ~GammaInfo(){};

     int order() { return order_ ; };

     std::shared_ptr<std::vector<std::string>> id_ranges(){ return id_ranges_ ;};
     std::shared_ptr<std::vector<std::string>> sigma_id_ranges(){ return sigma_id_ranges_ ;};

     std::string name() { return name_; };
     std::string sigma_name() { return sigma_name_; };

     std::vector<std::string> prev_sigmas() { return prev_sigmas_; };
     std::string prev_sigmas(int ii ) { return prev_sigmas_[ii]; };
     std::string prev_sigma_name() { return prev_sigmas_[0]; };

     std::vector<std::string> prev_gammas() { return prev_gammas_; };
     std::string prev_gammas(int ii ) { return prev_gammas_[ii] ; };
     std::string prev_gamma_name() { return prev_gammas_[0]; };

     std::string Bra_name() { return Bra_info_->name(); };
     std::string Ket_name() { return Ket_info_->name(); };
     std::string Prev_Bra_name() { return prev_Bra_info_->name(); };

     std::shared_ptr<CIVecInfo<double>> Bra_info() {return  Bra_info_; };
     std::shared_ptr<CIVecInfo<double>> prev_Bra_info(){ return prev_Bra_info_; };
     std::shared_ptr<CIVecInfo<double>> Ket_info() { return Ket_info_; };

     int Bra_nalpha() const {return  Bra_info_->nalpha(); };
     int prev_Bra_nalpha() const { return prev_Bra_info_->nalpha(); };
     int Ket_nalpha() const { return Ket_info_->nalpha(); };

     int Bra_nbeta() const {return  Bra_info_->nbeta(); };
     int prev_Bra_nbeta() const { return prev_Bra_info_->nbeta(); };
     int Ket_nbeta() const { return Ket_info_->nbeta(); };

};

class GammaIntermediate {

   public :
     std::shared_ptr<const std::vector<std::string>> full_id_ranges;
     std::shared_ptr<std::vector<int>> ids_pos;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos;
     std::shared_ptr<std::vector<int>> proj_id_order;
     int my_sign;

     GammaIntermediate( std::shared_ptr<const std::vector<std::string>> full_id_ranges_in,
                        std::shared_ptr<std::vector<int>> ids_pos_in,
                        std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                        int my_sign_in) :
     full_id_ranges(full_id_ranges_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), 
     proj_id_order(std::make_shared<std::vector<int>>(0)), my_sign(my_sign_in) {};


     GammaIntermediate( std::shared_ptr<const std::vector<std::string>> full_id_ranges_in,
                        std::shared_ptr<std::vector<int>> ids_pos_in,
                        std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                        int my_sign_in, std::shared_ptr<std::vector<int>> proj_id_order_in ) :
     full_id_ranges(full_id_ranges_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), 
     proj_id_order(proj_id_order_in), my_sign(my_sign_in) {};

     ~GammaIntermediate(){};

};

class GammaGenerator{
  friend GammaInfo;

  private : 
    // TODO Should replace this with the states infos  and pointer to CIVecInfo map,
    //      Do cycling over states outside of GammaGenerator.
    std::shared_ptr<StatesInfo<double>> target_states_;
    int Ket_num_;
    int Bra_num_;

    std::shared_ptr<const std::vector<bool>> orig_aops_ ;
    std::shared_ptr<const std::vector<std::string>> orig_ids_ ;

    std::shared_ptr<const std::vector<bool>> free_aops_;
    std::shared_ptr<const std::vector<std::string>> free_ids_;

    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec;

    std::shared_ptr<std::map< char, int>>         op_order ;
    std::shared_ptr<std::map< std::string, int>>  idx_order ;

    double bk_factor;

    // key    : name of this gamma
    // result : map containing names of relevant A-tensors, list of reorderings, and factor for each reordering
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, AContribInfo>> >> G_to_A_map;

    // key    : name of this gamma
    // result : information used here and in compute routines
    std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> Gamma_map;

    bool projected_bra_;
    bool projected_ket_; 

    std::shared_ptr<std::vector<std::string>> Bra_names_;
    std::shared_ptr<std::vector<std::string>> Ket_names_;

    // key    : name of A-tensor
    // result : list of reorderings which much be applied to this A-tensor before it is contracted with this gamma.
    //          second part of pair is factor associated with each reordering.
    std::shared_ptr<std::map<std::string, std::vector< std::pair<std::vector<int>, std::pair<int,int>> >> > Aid_orders_map;

  public :
    GammaGenerator( std::shared_ptr<StatesInfo<double>> target_states_, int Ket_num, int Bra_num,
                    std::shared_ptr<const std::vector<std::string>> orig_ids, std::shared_ptr< const std::vector<bool>> orig_aops,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>> Gamma_map_in,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, AContribInfo  >>>> G_to_A_map_in,
                    double bk_factor );

    ~GammaGenerator(){};


    //TODO make this private again when finished testing!!!
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec;

    void add_gamma(std::shared_ptr<range_block_info> block_info );

    bool norm_order();

    bool anti_norm_order();

    bool optimized_alt_order();

    bool generic_reorderer( std::string reordering_name, bool first_reordering, bool final_reordering );

    void normal_order( int kk );

    void anti_normal_order( int kk );

    void alternating_order( int kk );

    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name );

    bool check_if_same_sector( std::string bra_name, std::string ket_name );

    bool generic_reorderer_same_sector( std::string reordering_name, std::string bra_name,
                                        bool final_reordering );

    bool generic_reorderer_different_sector( std::string reordering_name, std::string bra_name,
                                             std::string ket_name, bool final_reordering );

    bool proj_onto_map( std::shared_ptr<GammaIntermediate> gint,
                        std::map<char, int> bra_hole_map, std::map<char, int> bra_elec_map,
                        std::map<char, int> ket_hole_map, std::map<char, int> ket_elec_map );

    void Contract_remaining_indexes(int kk);

    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec );

    std::shared_ptr<pint_vec> Standardize_delta_ordering(std::shared_ptr<pint_vec> deltas_pos );

    bool ordering( std::pair<std::string,std::string> a, std::pair<std::string,std::string> b ){
                   return( (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false ); }

    bool RangeCheck( std::shared_ptr<const std::vector<std::string>> full_id_ranges );

    bool Forbidden_Index( std::shared_ptr<const std::vector<std::string>> id_ranges, int position );

    bool all_active_ranges( std::shared_ptr<GammaIntermediate> gint);

    bool gamma_survives( std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<const std::vector<std::string>> id_ranges) ;

    static std::shared_ptr<std::vector<std::pair<int,int>>>
           Standardize_delta_ordering_generic(std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos ) ;

    std::vector<int> get_standard_range_order(const std::vector<std::string> &rngs) ;

    std::vector<int> get_position_order(const std::vector<int> &positions) ;

    std::vector<int> get_Aid_order ( const std::vector<int>& id_pos );

    std::vector<int> get_standard_order (const std::vector<std::string>& rngs );

    std::vector<int> get_standard_idx_order(const std::vector<std::string>& idxs) ;

    std::vector<int> get_standardized_alt_order( const std::vector<std::string>& rngs ,const std::vector<bool>& aops ) ;

    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec, std::string name );
    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec, std::string name,
                                    std::string bra_name, std::string ket_name );
};
#endif
