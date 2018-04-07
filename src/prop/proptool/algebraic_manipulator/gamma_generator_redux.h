#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_Redux_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_Redux_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>

using namespace WickUtils;

class GammaIntermediateRedux {

   public :
     std::shared_ptr<std::vector<int>> ids_pos;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos;
     int my_sign;

     GammaIntermediateRedux( std::shared_ptr<std::vector<int>> ids_pos_in,
                             std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in,
                             int my_sign_in ) :
     ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), my_sign(my_sign_in) {};

     ~GammaIntermediateRedux(){};

};

class GammaGeneratorRedux{
  friend GammaInfo;

  private : 
    std::shared_ptr<StatesInfo<double>> target_states_;
    std::shared_ptr<std::vector<std::string>> Bra_names_;
    std::shared_ptr<std::vector<std::string>> Ket_names_;

    std::shared_ptr<std::vector<bool>> orig_aops_ ;
    std::shared_ptr<std::vector<std::string>> orig_ids_ ;
    std::shared_ptr<std::vector<std::string>> orig_rngs_ ;

    std::shared_ptr<const std::vector<bool>> std_aops_ ;
    std::shared_ptr<const std::vector<std::string>> std_ids_ ;
    std::vector<std::string> std_rngs_ ;

    // key    : name of this gamma
    // result : map containing names of relevant A-tensors, list of reorderings, and factor for each reordering
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo>>> >> G_to_A_map;

    // key    : name of this gamma
    // result : information used here and in compute routines
    std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> Gamma_map;

    std::shared_ptr<TensOp_Base> total_op_;

    double bk_factor_;
    int orig_aops_half_size_;

    std::vector<int> standard_order_;
    std::vector<int> block_to_std_order_;
    std::vector<std::string> standardized_full_ids_;
    std::vector<std::string> standardized_full_id_ranges_;

    std::shared_ptr<std::vector<bool>> block_aops_;
    std::shared_ptr<std::vector<char>> block_aops_rngs_;
    std::vector<int> block_ids_pos_;
    std::vector<std::string> block_rngs_; 

    std::shared_ptr<std::vector<int>> aops_trans_;
    std::shared_ptr<std::vector<int>> rngs_trans_;
    std::shared_ptr<std::vector<int>> idxs_trans_;

    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateRedux>>> final_gamma_vec;
    
    // key    : name of A-tensor
    // result : list of reorderings which much be applied to this A-tensor before it is contracted with this gamma.
    //          second part of pair is factor associated with each reordering.
    std::shared_ptr<std::map<std::string, std::vector< std::pair<std::vector<int>, std::pair<int,int>> >> > Aid_orders_map;

  public :
    //TODO make this private again when finished testing!!!
    std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateRedux>>> gamma_vec;

    GammaGeneratorRedux( std::shared_ptr<StatesInfo<double>> target_states_, int Ket_num, int Bra_num,
                         std::shared_ptr<TensOp_Base> multitensop, 
                         std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>>& Gamma_map_in,
                         std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo>  >>>>& G_to_A_map_in,
                         double bk_factor );

    ~GammaGeneratorRedux(){};


    void add_gamma( const std::shared_ptr<Range_Block_Info> block_info, std::shared_ptr<std::vector<char>>  trans_aops_rngs,
                    std::shared_ptr<std::vector<bool>> trans_aops );

    void add_gamma( const std::shared_ptr<Range_Block_Info> block_info, const std::vector<std::string>& range_block );

    bool generic_reorderer( std::string reordering_name, bool first_reordering, bool final_reordering );

    bool generic_reorderer_different_sector( std::string reordering_name, std::string bra_name,
                                             std::string ket_name, bool final_reordering );
    bool anti_norm_order();

    void normal_order( int kk );

    void anti_normal_order( int kk );

    void alternating_order( int kk );

    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name );

    //TODO Replace this, but keep for now as very clear, if slow.
    bool proj_onto_map( std::shared_ptr<GammaIntermediateRedux> gint,
                        std::map<char, int> bra_hole_map, std::map<char, int> bra_elec_map,
                        std::map<char, int> ket_hole_map, std::map<char, int> ket_elec_map );

    void swap( int ii, int jj, int kk, std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateRedux>>> gamma_vec );

    std::shared_ptr<std::vector<std::pair<int,int>>>
    standardize_delta_ordering_generic(std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos );

    void set_standardized_alt_order_unranged ( int kk , std::vector<int>& standard_alt_order);

    //routines for reorderings
    std::vector<int> get_standard_order (const std::vector<std::string>& rngs );

    std::vector<int> get_standard_range_order(const std::vector<std::string> &rngs) ;

    std::vector<int> get_position_order(const std::vector<int> &positions) ;

    std::vector<int> get_Aid_order( const std::vector<int>& id_pos );

    //reordering of initial idxs
    void get_standard_idx_order_init();

    //reorders with respect to the range obtained in the above 
    std::vector<int> get_standard_idx_order(const std::vector<std::string>& idxs) ;

    std::vector<int> get_standardized_alt_order( const std::vector<std::string>& rngs ,const std::vector<bool>& aops ) ;

    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateRedux>>> final_gamma_vec, std::string name );
    void print_gamma_contributions( std::shared_ptr<std::vector<std::shared_ptr<GammaIntermediateRedux>>> final_gamma_vec, std::string name,
                                    std::string bra_name, std::string ket_name );

 
    void print_gamma_intermediate( std::shared_ptr<GammaIntermediateRedux> gint );

};
#endif
