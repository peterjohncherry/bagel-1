#ifndef __SRC_PROP_PROPTOOL_Gamma_Generator_Redux_H
#define __SRC_PROP_PROPTOOL_Gamma_Generator_Redux_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_base.h>

using namespace WickUtils;

template<typename DataType> 
class GammaIntermediateRedux : public GammaIntermediate_Base {

   public :
     GammaIntermediateRedux( std::shared_ptr<std::vector<int>> ids_pos,
                             std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos,
                             std::pair<double,double> factors ) :
                             GammaIntermediate_Base( ids_pos, deltas_pos, factors ) {};

    ~GammaIntermediateRedux(){};

};

template class GammaIntermediateRedux<double>;
template class GammaIntermediateRedux<std::complex<double>>;

template<typename DataType> 
class GammaGeneratorRedux : public GammaGenerator_Base {
  friend GammaInfo_Base;

  public :

    // key    : name of this gamma
    // result : map containing names of relevant A-tensors, list of reorderings, and factor for each reordering
    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo_Base>>> >> G_to_A_map;

    GammaGeneratorRedux<DataType>( std::shared_ptr<StatesInfo_Base> target_states, int Ket_num, int Bra_num,
                                   std::shared_ptr<TensOp_Base> multitensop, 
                                   std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo_Base>>>& Gamma_map_in,
                                   std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, std::shared_ptr<AContribInfo_Base>  >>>>& G_to_A_map_in,
                                   std::pair<double,double> bk_factor ):
                                   GammaGenerator_Base( target_states, Ket_num, Bra_num, multitensop, Gamma_map_in, bk_factor ), 
                                   G_to_A_map(G_to_A_map_in) {} ;

    ~GammaGeneratorRedux<DataType>(){};

    void add_gamma( const std::shared_ptr<Range_Block_Info> block_info, std::shared_ptr<std::vector<bool>> trans_aops );

    void add_Acontrib_to_map( int kk, std::string bra_name, std::string ket_name );

    void block_trans_test( int kk );
};
#endif
