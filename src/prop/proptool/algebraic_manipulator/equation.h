#ifndef __SRC_PROP_PROPTOOL_wICKTOOL_EQUATION_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_EQUATION_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/smith/tensor.h>

template<typename DataType>
class Equation_Base {

   public :

     std::string name_;
     std::string type_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_;

     Equation_Base( std::string name, std::string type, 
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                    std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map )
                    : name_(name), type_(type), term_braket_map_(term_braket_map), expression_term_map_(expression_term_map){} 

    ~Equation_Base(){};

     virtual void generate_expressions() = 0; 

}; 
 
// Generates an Equation object to evaluate all f_ij 
// f is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Value : public Equation_Base<DataType> {

   public :

     using Equation_Base<DataType>::name_;
     using Equation_Base<DataType>::type_;
     using Equation_Base<DataType>::term_braket_map_;
     using Equation_Base<DataType>::expression_term_map_;

     Equation_Value( std::string name, std::string type, 
                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map )
                     : Equation_Base<DataType>( name, type, term_braket_map, expression_term_map  ) {}  

    ~Equation_Value(){};

     void generate_expressions();

}; 
#endif

