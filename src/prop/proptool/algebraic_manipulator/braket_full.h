#ifndef __SRC_PROP_PROPTOOL_BRAKET_FULL_H
#define __SRC_PROP_PROPTOOL_BRAKET_FULL_H
#include <src/prop/proptool/algebraic_manipulator/braket.h>
template<typename DataType> 
class BraKet_Full : public BraKet_Base {

  public :

    BraKet_Full( std::shared_ptr<Op_Info> multiop_info, std::pair<double, double> factor, int bra_num, int ket_num,  std::string type) :
                 BraKet_Base( multiop_info, factor, bra_num, ket_num, type) {} 
   ~BraKet_Full(){};

    void generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string,std::shared_ptr<TensOp_Base>>> MT_map,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> G_to_A_map,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                              std::shared_ptr<StatesInfo_Base> target_states,
                                              std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map );
    
};

#endif
