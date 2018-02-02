#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>

namespace bagel {

namespace Expression_Computer{

template<typename DataType>
class Expression_Computer {
  private :  
    std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression<DataType>>>> expression_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> Gamma_info_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;

  public:
   
    std::shared_ptr<std::map< std::string, DataType>> scalar_results_map;
    std::shared_ptr<StatesInfo<DataType>> TargetsInfo;

    Expression_Computer( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer,
                         std::shared_ptr<std::map< std::string , std::shared_ptr<Expression<DataType>>>> expression_map, 
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map );
    ~Expression_Computer() {};
    
    void Evaluate_Expression(std::string expression_name );

    void print_AContraction_list(std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base>>> ACompute_list, std::string A_Contrib_name );

    bool check_AContrib_factors(AContribInfo& AC_info );
   
    void set_gamma_maps( std::string expression_name,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map );

    struct compare_string_length {
      bool operator()(const std::string& first, const std::string& second) {
          return first.size() > second.size();
      }
    };


};
}//end Expression_Computer namespace

}//end bagel namespace
#endif

