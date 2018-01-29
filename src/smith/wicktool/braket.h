#ifndef __SRC_SMITH_BRAKET_H
#define __SRC_SMITH_BRAKET_H
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/gamma_generator.h>
#include <src/smith/wicktool/tensop.h>
#include <src/smith/wicktool/states_info.h>

template<typename DataType> 
class BraKet{

     public :
       const std::vector<std::string> op_list_;
       const std::shared_ptr<std::vector<std::vector<int>>> op_state_ids_;
       const DataType factor_;
       const int bra_num_;
       const int ket_num_;
       const std::string type_ ; // should be "ci_deriv" or "full" 
       const std::string multiop_name_;

       BraKet( std::vector<std::string>& op_list, DataType factor, int bra_num, int ket_num, 
               std::shared_ptr<std::vector<std::vector<int>>> op_state_ids, std::string type) :
               op_list_(op_list), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
               op_state_ids_(op_state_ids), type_(type),
               multiop_name_(std::accumulate(op_list_.begin(), op_list_.end(), std::string(""))) {};
      ~BraKet(){};

       std::shared_ptr<MultiTensOp::MultiTensOp<DataType>> Total_Op_;

       void generate_gamma_Atensor_contractions( int bra_num, int ket_num, std::shared_ptr<std::map< std::string, std::vector<int>>> Op_state_ids, 
                                                 std::shared_ptr<std::map<std::string,std::shared_ptr<MultiTensOp::MultiTensOp<DataType>>>> MT_map,                
                                                 std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map,
                                                 std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> GammaMap,
                                                 std::shared_ptr<StatesInfo<DataType>> target_states ); 
        
       
};
#endif
