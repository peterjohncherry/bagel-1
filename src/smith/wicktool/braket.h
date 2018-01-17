#ifndef __SRC_SMITH_BraKet_H
#define __SRC_SMITH_BraKet_H
  #include <src/smith/wicktool/term.h>
  #include <src/smith/wicktool/wickutils.h>
  #include <src/smith/wicktool/gamma_generator.h>
  #include <src/smith/wicktool/tensop.h>
  #include <src/smith/wicktool/spin_manager.h>
  #include <src/smith/wicktool/states_info.h>

template<typename DataType> 
class BraKet{

      public:
        // Operator, only sub ranges are relevant        
        std::shared_ptr<MultiTensOp::MultiTensOp<DataType>> Total_Op_;

        // state numbers 
        int bra_num_;
        int ket_num_;       
 
        //factor; this is specific to the expression to which this BraKet object belongs
        DataType factor_;

        std::shared_ptr<StatesInfo<DataType>> target_states_;
                
        BraKet( BraKet_Init<DataType>& braket_info,  
                std::shared_ptr<std::map<std::string,std::shared_ptr<MultiTensOp::MultiTensOp<DataType>>>> MT_map,                
                std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map,
                std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> GammaMap,
                std::shared_ptr<StatesInfo<DataType>> target_states ); 
        
        ~BraKet(){};
       
};
#endif
