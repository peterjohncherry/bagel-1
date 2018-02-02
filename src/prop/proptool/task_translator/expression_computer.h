#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/gamma_computer.h>

namespace bagel {

namespace Expression_Computer{

template<class DataType>
class Expression_Computer {
  
  public:
   
    std::shared_ptr<const Dvec> civectors; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression<double>>>> Expression_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> Gamma_info_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> TensOp_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Gamma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Sigma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> CIvec_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map;

    std::shared_ptr<std::map< std::string, double>> scalar_results_map;
    std::shared_ptr<Gamma_Computer::Gamma_Computer> Gamma_Machine;
  
    std::shared_ptr<StatesInfo<double>> TargetsInfo;

    Expression_Computer( std::shared_ptr<const Dvec > civectors, 
                         std::shared_ptr<std::map< std::string , std::shared_ptr<Expression<double>>>> Expression_map, 
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> TensOp_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Gamma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Sigma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> CIvec_data_map );
    ~Expression_Computer() {};
    
    void Evaluate_Expression(std::string expression_name );

    void print_AContraction_list(std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base>>> ACompute_list, std::string A_Contrib_name );

    bool check_AContrib_factors(AContribInfo& AC_info );

    void Known_TensOp_Initializer( std::shared_ptr<TensOp::TensOp<DataType>> TensOp_info  );

    struct compare_string_length {
      bool operator()(const std::string& first, const std::string& second) {
          return first.size() > second.size();
      }
    };


};
}//end Expression_Computer namespace

}//end bagel namespace
#endif

