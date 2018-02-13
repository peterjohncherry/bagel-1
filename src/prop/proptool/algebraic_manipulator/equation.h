#ifndef __SRC_PROP_PROPTOOL_ALGMAN_EQUATION_H
#define __SRC_PROP_PROPTOOL_AlGMAN_EQUATION_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>  
#include <src/prop/proptool/algebraic_manipulator/expression.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/initialization/tensop_info.h>  

template<typename DataType>
class Equation_Base {

   protected :
     std::string name_;
     std::string type_;

     std::shared_ptr<StatesInfo<DataType>> states_info_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_;

     std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map_;

     // Maps are not defined in constructor so that equations can be independent,  needed for building of expression_map_
     std::shared_ptr< std::map <std::string, std::shared_ptr< GammaInfo >>> gamma_info_map_;
     std::shared_ptr< std::map <std::string, std::shared_ptr< BraKet<DataType>>>> braket_map_;
     std::shared_ptr< std::map <std::string, std::shared_ptr< std::vector<std::shared_ptr<CtrOp_base>>>>> ACompute_map_;

     std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp::TensOp<DataType>>>> T_map_;      
     std::shared_ptr< std::map< std::string, std::shared_ptr< MultiTensOp::MultiTensOp<DataType>>>> MT_map_;      

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DataType>>>> CTP_map_;      
     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> >>> CMTP_map_;  

//     std::shared_ptr<std::vector<std::<std::tuple<std::pair<std::string,int>>>>> term_val_list_;

   public :
     Equation_Base( std::string name, std::string type, std::shared_ptr<StatesInfo<DataType>> states_info, 
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map )
                    : name_(name), type_(type), states_info_(states_info), term_braket_map_(term_braket_map),
                      expression_term_map_(expression_term_map) {}

    ~Equation_Base(){};

     void set_maps( std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map,
                    std::shared_ptr< std::map<std::string, std::shared_ptr< GammaInfo >>> gamma_info_map,
                    std::shared_ptr< std::map<std::string, std::shared_ptr< std::vector<std::shared_ptr<CtrOp_base>>>>> ACompute_map,
                    std::shared_ptr< std::map<std::string, std::shared_ptr< TensOp::TensOp<DataType>>>> T_map,
                    std::shared_ptr< std::map<std::string, std::shared_ptr< MultiTensOp::MultiTensOp<DataType>>>> MT_map,
                    std::shared_ptr< std::map<std::string, std::shared_ptr< CtrTensorPart<DataType>>>> CTP_map,     
                    std::shared_ptr< std::map<std::string, std::shared_ptr< CtrMultiTensorPart<DataType> >>> CMTP_map);


     void generate_all_expressions();  

     std::shared_ptr<Expression<DataType>> build_expression (std::string expression_name ) ;
     std::shared_ptr<Expression<DataType>> build_expression( std::shared_ptr<std::vector<BraKet<DataType>>> expr_bk_list );

     std::shared_ptr< std::map <std::string, std::shared_ptr< GammaInfo >>> gamma_info_map() { return gamma_info_map_; }

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DataType>>>> CTP_map() { return CTP_map_; }      
     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> >>> CMTP_map() { return CMTP_map_; }  
     std::shared_ptr< std::map< std::string, std::shared_ptr< MultiTensOp::MultiTensOp<DataType>>>> MT_map() { return MT_map_    ; }      
     std::shared_ptr< std::map <std::string, std::shared_ptr< std::vector<std::shared_ptr<CtrOp_base>>>>>  ACompute_map() { return ACompute_map_; }
     std::shared_ptr< std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map() const  { return expression_term_map_; };
     std::shared_ptr< std::map <std::string, std::shared_ptr< Expression<DataType>>>> expression_map() const { return expression_map_; }

     std::string name() { return name_ ; }; 
     std::string type() { return type_ ; }; 

}; 
 
// Generates an Equation object to evaluate all f_ij 
// f is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Value : public Equation_Base<DataType> {

     using Equation_Base<DataType>::name_;
     using Equation_Base<DataType>::type_;
     using Equation_Base<DataType>::states_info_;
     using Equation_Base<DataType>::term_braket_map_;
     using Equation_Base<DataType>::expression_term_map_;

     using Equation_Base<DataType>::expression_map_;

     using Equation_Base<DataType>::gamma_info_map_;
     using Equation_Base<DataType>::braket_map_;
     using Equation_Base<DataType>::ACompute_map_;

     using Equation_Base<DataType>::T_map_ ;      
     using Equation_Base<DataType>::MT_map_ ;      

     using Equation_Base<DataType>::CTP_map_;      
     using Equation_Base<DataType>::CMTP_map_;  


   public :

     Equation_Value( std::string name, std::string type, std::shared_ptr<StatesInfo<DataType>> states_info, 
                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map )
                     : Equation_Base<DataType>( name, type, states_info, term_braket_map, expression_term_map  ) {}  

    ~Equation_Value(){};

     void generate_all_expressions(); 
}; 

#endif
