#ifndef __SRC_PROP_PROPTOOL_BRAKET_H
#define __SRC_PROP_PROPTOOL_BRAKET_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>

template<typename DataType> 
class BraKet{
  private : 
    std::string name_;

  public :
    std::vector<std::string> op_list_;
    std::vector<char> op_trans_list_;
    std::shared_ptr<std::vector<std::vector<int>>> op_state_ids_;
    DataType factor_;
    std::pair<DataType,DataType> ReIm_factors_;
    int bra_num_;
    int ket_num_;
    std::string type_ ; // should be "ci_deriv" or "full" 
    std::string multiop_name_;
    std::string proj_op_name_;        
    bool proj_op_;        
    bool projected_bra_;
    bool projected_ket_;

    std::vector<int> op_order_;
    std::shared_ptr<TensOp_Base> Total_Op_;

    std::shared_ptr<MultiOp_Info> multiop_info_;

    BraKet( std::shared_ptr<MultiOp_Info> multiop_info, std::pair<DataType,DataType> factor, int bra_num, int ket_num, std::string type);

    BraKet( std::vector<std::string>& op_list, std::vector<char>& op_trans_list,
            DataType factor, int bra_num, int ket_num, 
            std::shared_ptr<std::vector<std::vector<int>>> op_state_ids, std::string type); 

   ~BraKet(){};


    std::string bk_name() { return name_ ; } 
    std::string name() { return name_ ; } 
    DataType factor() const { return factor_ ; } 

    std::shared_ptr<std::vector<int>> op_order() { return multiop_info_->op_order_; };
    std::shared_ptr<std::vector<char>> op_trans_list(){ return multiop_info_->transformations_; };
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> op_state_ids() { return multiop_info_->state_ids_; };
   // void add_required_tens_block( std::string block_name ) { required_blocks.emplace( block_name ); } 
   // std::shared_ptr<std::set<std::string>> required_blocks( std::string block_name ) { required_blocks.emplace( block_name ); } 

    void generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string,std::shared_ptr<TensOp_Base>>> MT_map,                
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo<DataType>> >>>> G_to_A_map,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo<DataType> >>> gamma_info_map,
                                              std::shared_ptr<StatesInfo<DataType>> target_states,
                                              std::shared_ptr<std::set<std::string>> required_blocks,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map );         
    
   void print_gamma_Atensor_contractions(std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo<DataType>> >>>> G_to_A_map,
		                                        bool has_orb_exc );

};
#endif
