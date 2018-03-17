#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>

using namespace WickUtils;

class GammaIntermediate {

   public :
     std::shared_ptr<const std::vector<std::string>> full_id_ranges;
     std::shared_ptr<std::vector<int>> ids_pos;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos;
     int my_sign;
     long unsigned int plus_pnum_;
     long unsigned int kill_pnum_;

     GammaIntermediate( std::shared_ptr<const std::vector<std::string>> full_id_ranges_in,
                        std::shared_ptr<std::vector<int>> ids_pos_in,
                        std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                        int my_sign_in ) :
     full_id_ranges(full_id_ranges_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), 
     my_sign(my_sign_in) {};

     GammaIntermediate( std::shared_ptr<const std::vector<std::string>> full_id_ranges_in,
                        std::shared_ptr<std::vector<int>> ids_pos_in,
                        std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                        int my_sign_in, long unsigned int  plus_pnum, long unsigned int kill_pnum ) :
     full_id_ranges(full_id_ranges_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), 
     my_sign(my_sign_in), plus_pnum_(plus_pnum), kill_pnum_(kill_pnum) {};

     ~GammaIntermediate(){};

};

class GammaIntermediateUnranged {

   public :
     std::shared_ptr<std::vector<int>> ids_pos_;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_;
     int my_sign_;

     GammaIntermediateUnranged( std::shared_ptr<std::vector<int>> ids_pos,
                        std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos,
                        int my_sign ) :
     ids_pos_(ids_pos), deltas_pos_(deltas_pos), my_sign_(my_sign) {};

     ~GammaIntermediateUnranged(){};

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

    std::shared_ptr<std::vector<int>> free_pos_;
    std::shared_ptr<std::vector<bool>> free_aops_;
    std::shared_ptr<std::vector<std::string>> free_ids_;
    std::shared_ptr<std::vector<std::string>> free_ranges_;

    std::map<int,int> orig_to_free_pos_;

    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec;

    std::shared_ptr<std::map< char, int>>         op_order_ ;
    std::shared_ptr<std::map< std::string, int>>  idx_order_ ;
    std::vector<int> standardized_idx_order_;
    std::vector<int> order_map_;

    std::vector<int> standard_order_;

    double bk_factor;

    // key    : name of this gamma
    // result : map containing names of relevant A-tensors, list of reorderings, and factor for each reordering
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo>>> >> G_to_A_map;

    // key    : name of this gamma
    // result : information used here and in compute routines
    std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> Gamma_map;

    bool projected_bra_;
    bool projected_ket_; 
    bool orb_exc_deriv_; 
    char proj_op_name_;  


    std::shared_ptr<std::vector<std::string>> Bra_names_;
    std::shared_ptr<std::vector<std::string>> Ket_names_;

    // key    : name of A-tensor
    // result : list of reorderings which much be applied to this A-tensor before it is contracted with this gamma.
    //          second part of pair is factor associated with each reordering.
    std::shared_ptr<std::map<std::string, std::vector< std::pair<std::vector<int>, std::pair<int,int>> >> > Aid_orders_map;

  public :
    GammaGenerator( std::shared_ptr<StatesInfo<double>> target_states_, int Ket_num, int Bra_num,
                    std::shared_ptr<const std::vector<std::string>> orig_ids, std::shared_ptr< const std::vector<bool>> orig_aops,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>>& Gamma_map_in,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo>  >>>>& G_to_A_map_in,
                    double bk_factor );

    ~GammaGenerator(){};


    //TODO make this private again when finished testing!!!
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec;
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateUnranged>>> gamma_vec_unranged_;

    void add_gamma(std::shared_ptr<Range_Block_Info> block_info );

    bool norm_order();

    bool anti_norm_order();

    bool optimized_alt_order();

    bool generic_reorderer( std::string reordering_name, bool first_reordering, bool final_reordering );

    void normal_order( int kk );

    void anti_normal_order( int kk );

    void alternating_order( int kk );

    void normal_order_unranged( int kk );

    void anti_normal_order_unranged( int kk );

    void alternating_order_unranged( int kk );


    void set_standardized_alt_order_unranged ( int kk, std::vector<int>& alt_order_unranged );

    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name );

    void add_Acontrib_to_map_orb_deriv( int kk, std::string bra_name, std::string ket_name );
 
    bool check_if_same_sector( std::string bra_name, std::string ket_name );

    bool generic_reorderer_same_sector( std::string reordering_name, std::string bra_name,
                                        bool final_reordering );

    bool generic_reorderer_different_sector( std::string reordering_name, std::string bra_name,
                                             std::string ket_name, bool final_reordering );

    void swap_unranged( int ii, int jj, int kk );
    bool generic_reorderer_unranged( std::string reordering_name, bool first_reordering, bool final_reordering );
    bool generic_reorderer_different_sector_unranged( std::string reordering_name, std::string bra_name,
                                                      std::string ket_name, bool final_reordering         );

    bool proj_onto_map( std::shared_ptr<GammaIntermediate> gint,
                        std::map<char, int> bra_hole_map, std::map<char, int> bra_elec_map,
                        std::map<char, int> ket_hole_map, std::map<char, int> ket_elec_map );

//    bool braket_survival_check( std::shared_ptr<GammaIntermediate> gint, std::shared_ptr<CIVecInfo> bra_info, 
//                                std::shared_ptr<CIVecInfo> ket_info );

    void Contract_remaining_indexes(int kk);

    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> gamma_vec );

    std::shared_ptr<pint_vec> Standardize_delta_ordering(std::shared_ptr<pint_vec> deltas_pos );

    bool ordering( std::pair<std::string,std::string> a, std::pair<std::string,std::string> b ){
                   return( (idx_order_->at(a.first) < idx_order_->at(b.first) ) ? true : false ); }

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

    //reordering of initial idxs
    void get_standard_idx_order_init();

    //reorders with respect to the range obtained in the above 
    std::vector<int> get_standard_idx_order(const std::vector<std::string>& idxs) ;

    std::vector<int> get_standardized_alt_order( const std::vector<std::string>& rngs ,const std::vector<bool>& aops ) ;

    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec, std::string name );
    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediate>>> final_gamma_vec, std::string name,
                                    std::string bra_name, std::string ket_name );
};
#endif
