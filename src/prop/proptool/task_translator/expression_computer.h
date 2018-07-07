#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EXPRESSIONCOMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/algebraic_manipulator/expression_full.h>
#include <src/prop/proptool/algebraic_manipulator/expression_orb_exc_deriv.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>
#include <src/prop/proptool/integrals/moint_computer.h>


namespace bagel {

namespace Expression_Computer{

template<typename DataType>
class Expression_Computer {
  private :  
    std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression<DataType>>>> expression_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> Gamma_info_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;
    std::shared_ptr<MOInt_Computer<DataType>> moint_computer_;

    std::shared_ptr<TensOp_Computer::TensOp_Computer<DataType>> tensop_machine_;
                                                                                 
  public:
   
    std::shared_ptr<std::map< std::string, DataType>> scalar_results_map;

    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensor_results_map_;

    std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<std::string, DataType>>>> gamma_contrib_maps_; 

    std::shared_ptr<StatesInfo<DataType>> TargetsInfo;

    Expression_Computer( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer,
                         std::shared_ptr<std::map< std::string , std::shared_ptr<Expression<DataType>>>> expression_map, 
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map,
                         std::shared_ptr<MOInt_Computer<DataType>> moint_computer  );
    ~Expression_Computer() {};
    
    void evaluate_expression( std::shared_ptr<Expression<DataType>> expression );

    void evaluate_expression( std::string expression_name );

    void evaluate_expression_orb_exc_deriv( std::shared_ptr<Expression<DataType>> expression );

    void evaluate_expression_full( std::shared_ptr<Expression<DataType>> expression );

    bool check_acontrib_factors(AContribInfo_Base& AC_info );
   
    void set_gamma_maps( std::string expression_name,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map,
                         std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map );

    void print_scalar_result( std::string expression_name, bool print_gamma_contribs = true );

    //TESTING
    void test_trace_subtraction();
    void test_sum_reordered_tensor_list();
    void check_rdms();
};

}//end Expression_Computer namespace
}//end bagel namespace
#endif

