#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQUATIONCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQUATIONCOMPUTER_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/equation.h>

template<typename DataType>
class Equation_Computer_Base {

   public :

     std::string name_;
     std::string type_;

     Equation_Computer_Base(); 
    ~Equation_Computer_Base(){};

}; 
 
// Generates an Equation object to evaluate all f_ij 
// f is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Computer_Value : public Equation_Computer_Base<DataType> {

   using Equation_Computer_Base<DataType>::name_;
   using Equation_Computer_Base<DataType>::type_;

   public :

     Equation_Computer_Value(); 
    ~Equation_Computer_Value(){};

}; 

template<typename DataType>
class Equation_Computer_LinearRM : public Equation_Computer_Base<DataType> {

   using Equation_Computer_Base<DataType>::name_;
   using Equation_Computer_Base<DataType>::type_;

   public :

     Equation_Computer_LinearRM(); 
    ~Equation_Computer_LinearRM(){};

}; 
#endif
