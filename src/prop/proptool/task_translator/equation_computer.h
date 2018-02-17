#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQUATIONCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_EQUATIONCOMPUTER_H

#include <cstdarg>
#include <src/util/math/matrix.h>
#include <src/prop/proptool/proputils.h>
#include <src/smith/tensor.h>
#include <src/smith/multitensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/task_translator/expression_computer.h>
#include <src/prop/proptool/algebraic_manipulator/system_info.h>

namespace bagel { 
template<typename DataType>
class Equation_Computer_Base {
 
   protected :
     const std::string name_;
     const std::string type_;
     const std::shared_ptr<Equation_Base<DataType>> equation_;
     const std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;

     std::shared_ptr<Expression_Computer::Expression_Computer<DataType>> expression_computer_;
     std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer_;

     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map_; 
     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map_;

     DataType
     get_scalar_result( std::string result_name, std::vector<std::pair<std::string, int>>& fixed_idxs  );

     DataType
     get_scalar_result( std::string result_name, std::vector<std::pair<std::string, int>>& fixed_idxs,
                                                 std::vector<std::pair<std::string, int>>& summed_idxs );
 
     std::shared_ptr<SMITH::Tensor_<DataType>>
     get_tensop( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idxs  );
   
     std::shared_ptr<SMITH::Tensor_<DataType>>
     get_tensop( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idxs,
                                          std::vector<std::pair<std::string, int>>& summed_idxs );
  
     std::shared_ptr<SMITH::MultiTensor_<DataType>>
     get_tensop_vector( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idxs  );

     std::shared_ptr<SMITH::MultiTensor_<DataType>>
     get_tensop_vector( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idxs,
                                          std::vector<std::pair<std::string, int>>& summed_idxs );

     void 
     evaluate_term( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idxs );

     void 
     evaluate_term( std::string tensop_name, std::vector<std::pair<std::string, int>>& fixed_idx,
                                             std::vector<std::pair<std::string, int>>& summed_idxs );

   public :

     Equation_Computer_Base( std::shared_ptr<Equation_Base<DataType>> equation,
                             std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map ) : 
                             name_(equation->name()+"computer"), type_("Base"), equation_(equation),
                             range_conversion_map_(range_conversion_map){} ; 

    ~Equation_Computer_Base(){};

     std::string const name() { return name_; } 
     std::string const type() { return type_; } 

     void set_computers( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer );

     void set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,  
                    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map );

     void build_expression_computer(); 
     virtual void solve_equation() = 0; 

}; 
 
// Generates an Equation object to evaluate all f_ij 
// f is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Computer_Value : public Equation_Computer_Base<DataType> {

   using Equation_Computer_Base<DataType>::name_;
   using Equation_Computer_Base<DataType>::type_;
   using Equation_Computer_Base<DataType>::equation_;
   using Equation_Computer_Base<DataType>::range_conversion_map_;
   using Equation_Computer_Base<DataType>::gamma_computer_;
   using Equation_Computer_Base<DataType>::tensop_data_map_;
   using Equation_Computer_Base<DataType>::expression_computer_;

   public :

     Equation_Computer_Value( std::shared_ptr<Equation_Base<DataType>> equation,
                              std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map ) : 
                              Equation_Computer_Base<DataType>( equation, range_conversion_map ) {  
                              assert(std::dynamic_pointer_cast<Equation_Value<DataType>>(equation)); }; 


    ~Equation_Computer_Value(){};

     void solve_equation(); 
}; 
 
}
#endif
