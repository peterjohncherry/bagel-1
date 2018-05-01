#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_orbderiv_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_orbderiv_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_base.h>

using namespace WickUtils;

template<typename DataType> 
class GammaIntermediate_OrbExcDeriv : public GammaIntermediate_Base {

   public :
     
     std::shared_ptr<std::vector<std::pair<int,int>>> target_A_deltas_pos_;
     std::shared_ptr<std::vector<std::pair<int,int>>> target_target_deltas_pos_;

     GammaIntermediate_OrbExcDeriv( std::shared_ptr<std::vector<int>> ids_pos,
                                 std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos,
                                 std::shared_ptr<std::vector<std::pair<int,int>>> target_A_deltas_pos,
                                 std::shared_ptr<std::vector<std::pair<int,int>>> target_target_deltas_pos,
                                 std::pair<double,double> factors ) :
                                 GammaIntermediate_Base( ids_pos, deltas_pos, factors),
                                 target_A_deltas_pos_(target_A_deltas_pos), target_target_deltas_pos_(target_target_deltas_pos) {} ;

     ~GammaIntermediate_OrbExcDeriv(){};

     std::shared_ptr<std::vector<std::pair<int,int>>> A_A_deltas_pos() { return deltas_pos_  ; } 
     std::shared_ptr<std::vector<std::pair<int,int>>> target_A_deltas_pos() { return target_A_deltas_pos_  ; } 
     std::shared_ptr<std::vector<std::pair<int,int>>> target_target_deltas_pos() { return target_target_deltas_pos_;} 

};

template class GammaIntermediate_OrbExcDeriv<double>;
template class GammaIntermediate_OrbExcDeriv<std::complex<double>>;

template<typename DataType> 
class GammaGenerator_OrbExcDeriv : public GammaGenerator_Base {
  friend GammaInfo_Base;

  public :
    std::shared_ptr<std::map<std::string,
                    std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> >> block_G_to_A_map_;

    std::string target_op_;

    std::string target_block_name_;

    std::vector<int> target_block_id_pos_;
    int target_block_end_;
    int target_block_start_;

    GammaGenerator_OrbExcDeriv( std::shared_ptr<StatesInfo_Base> target_states, int Ket_num, int Bra_num,
                                std::shared_ptr<TensOp_Base> multitensop, 
                                std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo_Base>>>& Gamma_map_in,
                                std::shared_ptr<std::map<std::string,
                                                std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> >>& block_G_to_A_map,
                                std::pair<double,double> bk_factor, std::string target_op ) :
                                GammaGenerator_Base( target_states, Ket_num, Bra_num, multitensop, Gamma_map_in, bk_factor ), 
                                block_G_to_A_map_(block_G_to_A_map), target_op_(target_op) {} ;

    ~GammaGenerator_OrbExcDeriv(){};

    void add_gamma( const std::shared_ptr<Range_Block_Info> block_info, std::shared_ptr<std::vector<bool>> trans_aops  );

    void swap( int ii, int jj, int kk );

    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name );

    void transformation_tester( int kk  );
};
#endif
