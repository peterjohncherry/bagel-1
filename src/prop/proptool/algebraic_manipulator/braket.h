#ifndef __SRC_PROP_PROPTOOL_BRAKET_H
#define __SRC_PROP_PROPTOOL_BRAKET_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/a_contrib_info.h>

class BraKet_Base{
  protected : 
    std::string name_;

  public :
    std::vector<std::string> op_list_;
    std::vector<char> op_trans_list_;
    std::shared_ptr<std::vector<std::vector<int>>> op_state_ids_;
    int bra_num_;
    int ket_num_;
    std::string type_ ; // should be "ci_deriv" or "full" 
    std::string multiop_name_;
    std::string target_op_;        
    bool proj_op_;        
    bool projected_bra_;
    bool projected_ket_;
    bool orb_exc_deriv_;

    std::pair<double,double> factor_;
    std::pair<double,double> ReIm_factors_;
   
    std::vector<int> op_order_;
    std::shared_ptr<TensOp_Base> Total_Op_;
   
    std::shared_ptr<MultiOp_Info> multiop_info_;
   
    BraKet_Base( std::shared_ptr<MultiOp_Info> multiop_info, std::pair<double, double> factor, int bra_num, int ket_num, std::string type);
  
   ~BraKet_Base(){};

    std::string bk_name() { return name_ ; } 
    std::string name() { return name_ ; } 
    std::pair<double,double> factor() const { return factor_ ; } 

    std::shared_ptr<std::vector<int>> op_order() { return multiop_info_->op_order_; };
    std::shared_ptr<std::vector<char>> op_trans_list(){ return multiop_info_->transformations_; };
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> op_state_ids() { return multiop_info_->state_ids_; };

};


template<typename DataType> 
class BraKet : public BraKet_Base {

  public :

    BraKet( std::shared_ptr<MultiOp_Info> multiop_info, std::pair<double,double> factor, int bra_num, int ket_num, std::string type);

   ~BraKet(){};

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
