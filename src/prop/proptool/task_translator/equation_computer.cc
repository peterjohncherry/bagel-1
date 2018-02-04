#include <src/prop/proptool/task_translator/equation_computer.h>

using namespace std;
using namespace bagel;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_computers( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_computers " << endl;

  gamma_computer_ = gamma_computer;

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,  
                                                 std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_maps() " << endl;
 
  gamma_data_map_ = gamma_data_map;
  tensop_data_map_ = tensop_data_map;
  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::build_expression_computer(){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::build_expression_computer " << endl;

  cout << "equation_->expression_map()->size() = "<< equation_->expression_map()->size() << endl; 
  expression_computer_ =
  make_shared<Expression_Computer::Expression_Computer<DataType>>( gamma_computer_, equation_->expression_map(), range_conversion_map_, tensop_data_map_ );
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Value<DataType>::solve_equation(){ /*not really solving anything here*/ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::solve_equation() " << endl;
  for ( auto& expr_map_elem : *equation_->expression_map() )
    expression_computer_->evaluate_expression( expr_map_elem.first ) ; 

  return; 
}
////////////////////////////////////////////////////////////////////
template class Equation_Computer_Base<double>;
template class Equation_Computer_Value<double>;
template class Equation_Computer_LinearRM<double>;
////////////////////////////////////////////////////////////////////
