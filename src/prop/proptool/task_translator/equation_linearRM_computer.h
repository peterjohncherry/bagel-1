#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQN_LRM_COMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQN_LRM_COMPUTER_H
#include <src/prop/proptool/task_translator/equation_computer.h>

namespace bagel { 

template<typename DataType>
class Equation_Computer_LinearRM : public Equation_Computer_Base<DataType> {

   using Equation_Computer_Base<DataType>::name_;
   using Equation_Computer_Base<DataType>::type_;
   using Equation_Computer_Base<DataType>::equation_;
   using Equation_Computer_Base<DataType>::range_conversion_map_;
   using Equation_Computer_Base<DataType>::gamma_computer_;
   using Equation_Computer_Base<DataType>::tensop_data_map_;
   using Equation_Computer_Base<DataType>::expression_computer_;
 
   //TODO replace with MatType thing so can handle double and complex<double>
   std::shared_ptr<Matrix> h_eff;

   int ref_space_dim;

   public :

     Equation_Computer_LinearRM( std::shared_ptr<Equation_LinearRM<DataType>> equation,
                                 std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map ) : 
                                 Equation_Computer_Base<DataType>( equation, range_conversion_map ) {};

    ~Equation_Computer_LinearRM(){};

     void solve_equation(); 
}; 
}
#endif
