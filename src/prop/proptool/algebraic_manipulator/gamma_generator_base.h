#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_Base_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_Base_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>

using namespace WickUtils;

class GammaIntermediate_Base {

   public :
     std::shared_ptr<std::vector<int>> ids_pos_;
     std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_;
     std::pair<double,double> factors_;

     GammaIntermediate_Base( std::shared_ptr<std::vector<int>> ids_pos,
                             std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos,
                             std::pair<double,double> factors ) :
     ids_pos_(ids_pos), deltas_pos_(deltas_pos), factors_(factors) {};

     ~GammaIntermediate_Base(){};


     virtual
     std::shared_ptr<std::vector<std::pair<int,int>>> A_A_deltas_pos() {
       throw std::logic_error( "Should not access (A, A) contractions from gammaintermediate_base");
       return deltas_pos_  ;
     } 

     virtual
     std::shared_ptr<std::vector<std::pair<int,int>>> target_A_deltas_pos() {
       throw std::logic_error( "Should not access (target, A) contractions from gammaintermediate_base");
       return deltas_pos_  ;
     } 

     virtual
     std::shared_ptr<std::vector<std::pair<int,int>>> target_target_deltas_pos() {
       throw std::logic_error( "Should not access (target, target) contractions from gammaintermediate_base");
       return deltas_pos_;
     } 
};

class GammaGenerator_Base{
  friend GammaInfo_Base;

  protected : 
    std::shared_ptr<StatesInfo_Base> target_states_;
    std::shared_ptr<std::vector<std::string>> Bra_names_;
    std::shared_ptr<std::vector<std::string>> Ket_names_;

    std::shared_ptr<std::vector<bool>> orig_aops_ ;
    std::shared_ptr<std::vector<std::string>> orig_ids_ ;
    std::shared_ptr<std::vector<std::string>> orig_rngs_ ;

    std::shared_ptr<const std::vector<bool>> std_aops_ ;
    std::shared_ptr<const std::vector<std::string>> std_ids_ ;
    std::vector<std::string> std_rngs_ ;

    std::string bra_name_;
    std::string ket_name_;

    std::shared_ptr<std::map<char,int>> bra_hole_map_;
    std::shared_ptr<std::map<char,int>> ket_hole_map_;
    std::shared_ptr<std::map<char,int>> bra_elec_map_;
    std::shared_ptr<std::map<char,int>> ket_elec_map_;

    // key    : name of this gamma
    // result : information used here and in compute routines
    std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo_Base >>> Gamma_map;

    std::shared_ptr<TensOp_Base> total_op_;

    std::pair<double,double> bk_factor_;
    int orig_aops_half_size_;

    std::vector<int> standard_order_;
    std::vector<int> block_to_std_order_;
    std::vector<std::string> standardized_full_ids_;
    std::vector<std::string> standardized_full_id_ranges_;

    std::shared_ptr<std::vector<bool>> block_aops_;
    std::shared_ptr<std::vector<char>> block_aops_rngs_;
    std::shared_ptr<std::vector<std::string>> block_rngs_; 
    std::vector<int> block_ids_pos_;
    
    std::vector<std::string> block_idxs_;

    std::shared_ptr<std::vector<int>> aops_trans_;
    std::shared_ptr<std::vector<int>> rngs_trans_;
    std::shared_ptr<std::vector<int>> idxs_trans_;
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_;

    std::shared_ptr<std::vector< std::shared_ptr<GammaIntermediate_Base>>> final_gamma_vec_;
    
    // key    : name of A-tensor
    // result : list of reorderings which much be applied to this A-tensor before it is contracted with this gamma.
    //          second part of pair is factor associated with each reordering.
    std::shared_ptr<std::map<std::string, std::vector< std::pair<std::vector<int>, std::pair<int,int>> >> > Aid_orders_map;

  public :
    std::shared_ptr<std::vector< std::shared_ptr<GammaIntermediate_Base> >> gamma_vec_;

    std::shared_ptr<Op_Info> op_info_;

    GammaGenerator_Base( std::shared_ptr<StatesInfo_Base> target_states_, int Ket_num, int Bra_num,
                                   std::shared_ptr<TensOp_Base> multitensop, 
                                   std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo_Base>>>& Gamma_map_in,
                                   std::pair<double,double> bk_factor );

    ~GammaGenerator_Base(){};

    std::shared_ptr<const std::vector<std::string>> std_ids() { return std_ids_; }

    std::shared_ptr<std::vector<int>> idxs_trans() { return idxs_trans_; }

    std::vector<std::string>& block_idxs() { return block_idxs_; }

    std::vector<int>& block_ids_pos() { return block_ids_pos_; }

    bool generic_reorderer( std::string reordering_name, bool first_reordering, bool final_reordering );

    bool generic_reorderer_different_sector( std::string reordering_name, bool final_reordering );

    void normal_order();

    void anti_normal_order();

    void alternating_order();

    //TODO Replace this, but keep for now as very clear, if slow.
    bool proj_onto_map( std::shared_ptr<GammaIntermediate_Base> gint,
                        std::map<char, int> bra_hole_map, std::map<char, int> bra_elec_map,
                        std::map<char, int> ket_hole_map, std::map<char, int> ket_elec_map );

    std::shared_ptr<std::vector<std::pair<int,int>>>
    standardize_delta_ordering_generic(std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos );

    void set_standardized_alt_order_unranged ( std::shared_ptr<GammaIntermediate_Base>& gint , std::vector<int>& standard_alt_order);

    std::vector<int> get_Aid_order( const std::vector<int>& id_pos );

    std::vector<int> get_position_order(const std::vector<int> &positions) ;
 
    void print_gamma_intermediate( std::shared_ptr<GammaIntermediate_Base> gint );

    virtual void add_gamma( const std::shared_ptr<Range_Block_Info> block_info, std::shared_ptr<std::vector<bool>> trans_aops );

    virtual void swap( int ii, int jj, int kk );

    virtual
    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name ) {
          throw std::logic_error("Add Acontrib to map is term specific; cannot call from gamma_generator_base!! Aborting!!"); };
};
#endif
