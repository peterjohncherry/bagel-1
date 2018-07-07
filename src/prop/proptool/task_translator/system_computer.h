#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_SYSTEMCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_SYSTEMCOMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/multitensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/integrals/moint_computer.h>
#include <src/prop/proptool/algebraic_manipulator/system_info.h>
#include <src/prop/proptool/task_translator/equation_computer.h>
#include <src/prop/proptool/task_translator/equation_computer_linearRM.h>

namespace bagel {
namespace  System_Computer { 
template<class DataType> 
class System_Computer {
  private: 
    // Tensor operators in MO basis and maps 
    std::shared_ptr<SMITH::Tensor_<DataType>> h1_; 
    std::shared_ptr<SMITH::Tensor_<DataType>> f1_; 
    std::shared_ptr<SMITH::Tensor_<DataType>> v2_; 

    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map_;

    std::shared_ptr<System_Info<DataType>> system_info_;
    std::shared_ptr<MOInt_Computer<DataType>> moint_computer_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_ ;
    std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer_;

  public:
    System_Computer( std::shared_ptr<System_Info<DataType>> system_info, std::shared_ptr<MOInt_Computer<DataType>> moint_computer,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map, 
                     std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_machine );
   ~System_Computer(){};

    void build_equation_computer(std::string equation_name );
    void build_tensop( std::string tensop_name ) ;

    std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer(){ return b_gamma_computer_; } 
};
}
}
#endif
