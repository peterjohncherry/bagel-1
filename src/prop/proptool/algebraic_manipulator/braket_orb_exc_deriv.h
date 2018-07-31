#ifndef __SRC_PROP_PROPTOOL_BRAKET_ORB_EXC_DERIV_H
#define __SRC_PROP_PROPTOOL_BRAKET_ORB_EXC_DERIV_H
#include <src/prop/proptool/algebraic_manipulator/braket.h>

template<typename DataType> 
class BraKet_OrbExcDeriv : public BraKet_Base {

  public :

    std::string target_op_;        

    BraKet_OrbExcDeriv( std::shared_ptr<Op_Info> multiop_info, std::pair<double,double> factor, int bra_num, int ket_num, std::string type) :
                        BraKet_Base( multiop_info, factor, bra_num, ket_num, type) {} 

   ~BraKet_OrbExcDeriv(){};

    void generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string, std::shared_ptr<TensOp_Base>>> MT_map,                
                                              std::shared_ptr<std::map<std::string,
                                                              std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo_Base> >>>> >> block_G_to_A_map,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo_Base >>> gamma_info_map,
                                              std::shared_ptr<StatesInfo_Base> target_states,
                                              std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map ); 
    

   std::string target_op(){ return target_op_;  }

};


#endif
