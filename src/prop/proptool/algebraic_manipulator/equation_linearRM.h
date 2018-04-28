#ifndef __SRC_PROP_PROPTOOL_ALGMAN_EQUATION_linearRM_H
#define __SRC_PROP_PROPTOOL_AlGMAN_EQUATION_linearRM_H
#include <src/prop/proptool/algebraic_manipulator/equation.h>
 
// Generates an Equation object to evaluate all f_ij 
// if is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_LinearRM : public Equation_Base<DataType> {

    using Equation_Base<DataType>::name_;
    using Equation_Base<DataType>::type_;
    using Equation_Base<DataType>::states_info_;
    using Equation_Base<DataType>::term_braket_map_;
    using Equation_Base<DataType>::expression_term_map_;

    using Equation_Base<DataType>::expression_map_;

    using Equation_Base<DataType>::gamma_info_map_;
    using Equation_Base<DataType>::braket_map_;
    using Equation_Base<DataType>::ACompute_map_;

    using Equation_Base<DataType>::MT_map_ ;      

    using Equation_Base<DataType>::CTP_map_;      

    using Equation_Base<DataType>::term_braket_map_state_spec_; 
    using Equation_Base<DataType>::expression_term_map_state_spec_;
    using Equation_Base<DataType>::term_map_; 

  public :
                   
    Equation_LinearRM( std::string name, std::string type, std::shared_ptr<StatesInfo<DataType>> states_info, 
                       std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>>  term_braket_map,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map, 
                       std::shared_ptr<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>,
                                                      std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>> term_braket_map_state_spec,
                       std::shared_ptr<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>, 
                                                 std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_state_spec ) 
                    : Equation_Base<DataType>( name, type, states_info, term_braket_map, expression_term_map, term_braket_map_state_spec, expression_term_map_state_spec){ cout << " eqn_lrm" << endl;}  

   ~Equation_LinearRM(){};

    void generate_state_specific_terms(); 
  
    void add_term( std::pair<std::string, std::vector<std::pair<std::string,int>>>&  new_key, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>> expr_bk_list );
}; 
#endif
