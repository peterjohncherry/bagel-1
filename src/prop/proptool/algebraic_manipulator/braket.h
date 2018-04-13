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
    const std::vector<std::string> op_list_;
    const std::vector<char> op_trans_list_;
    const std::shared_ptr<std::vector<std::vector<int>>> op_state_ids_;
    const DataType factor_;
    const int bra_num_;
    const int ket_num_;
    const std::string type_ ; // should be "ci_deriv" or "full" 
    const std::string multiop_name_;
    const std::string proj_op_name_;        
    const bool proj_op_;        
    bool projected_bra_;
    bool projected_ket_;

    std::vector<int> op_order_;
    std::shared_ptr<TensOp_Base> Total_Op_;

    std::shared_ptr<MultiOp_Info> multiop_info_;

    BraKet( std::vector<std::string>& op_list, std::vector<char>& op_trans_list,
            DataType factor, int bra_num, int ket_num, 
            std::shared_ptr<std::vector<std::vector<int>>> op_state_ids, std::string type); 

    BraKet( std::vector<std::string>& op_list, std::vector<char>& op_trans_list, DataType factor, int bra_num, int ket_num, 
            std::shared_ptr<std::vector<std::vector<int>>> op_state_ids, std::string type,
            std::string proj_op_name) :
            op_list_(op_list), op_trans_list_(op_trans_list), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
            op_state_ids_(op_state_ids), type_(type),
            multiop_name_(std::accumulate(op_list_.begin(), op_list_.end(), std::string(""))),
            proj_op_(true), proj_op_name_(proj_op_name) {
                

            if (type_[0] == 'c' )// checking if derivative  
              name_ = "c_{I}"; 
 
            name_ = "<" + std::to_string(bra_num)+ "| ";
            
              for ( int ii = 0 ; ii != op_list_.size(); ii++ ) {
                WickUtils::print_vector(op_state_ids_->at(ii), "op_state_ids->at("+std::to_string(ii)+")" ); std::cout << std::endl;
                name_ += op_list_[ii] ;
                if (op_state_ids_->at(ii).size() > 0 ) {
                  name_ +=  "^{"; 
                  for( int jj = 0; jj != op_state_ids_->at(ii).size(); jj++ ) 
                    name_ += std::to_string(op_state_ids_->at(ii)[jj]); 
                  name_ += "}"; 
                }
              }
              name_ += " |"+ std::to_string(ket_num) + ">";
               
              projected_bra_ = false;
              projected_ket_ = false;
            }; 

   ~BraKet(){};


    std::string bk_name() { return name_ ; } 
    std::string name() { return name_ ; } 
    DataType factor() const { return factor_ ; } 

   // void add_required_tens_block( std::string block_name ) { required_blocks.emplace( block_name ); } 
   // std::shared_ptr<std::set<std::string>> required_blocks( std::string block_name ) { required_blocks.emplace( block_name ); } 

    void generate_gamma_Atensor_contractions( std::shared_ptr<std::map<std::string,std::shared_ptr<TensOp_Base>>> MT_map,                
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo> >>>> G_to_A_map,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> gamma_info_map,
                                              std::shared_ptr<StatesInfo<DataType>> target_states,
                                              std::shared_ptr<std::set<std::string>> required_blocks,
                                              std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> ctp_map );         
    
   void print_gamma_Atensor_contractions(std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, std::shared_ptr<AContribInfo> >>>> G_to_A_map,
		                                        bool has_orb_exc );

};
#endif
